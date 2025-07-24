import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import yaml
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


from Dataloader import Dataloader
from NetUtils import TrainableSiLU
from LossUtils import QuantileLoss, MultiheadLoss
from LayerUtils import EnergyLayer, QuantileOrderingLayer


from MiscUtils import *
from PlotUtils import PlotMetric, PlotHist, PlotCompare2D


def main():
    # save all params used for training in a dict and after training dupm to yaml
    params = {}

    files = []
    with open('dl_train_files.txt', 'r') as file_cfg:
        files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    # dataloader.Load(files)
    dataloader.Load('../train_data/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M_800/nano_0.root')
    X_train, input_names, y_train, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)

    params['input_names'] = input_names
    params['target_names'] = target_names

    standardize = False
    params['standardize'] = standardize
    params['target_train_means'] = None
    params['target_train_scales'] = None
    params['input_train_means'] = None
    params['input_train_scales'] = None
    
    input_scaler = None
    target_scaler = None
    if standardize:
        input_scaler = StandardScaler()
        X_train = input_scaler.fit_transform(X_train)
        target_scaler = StandardScaler()
        y_train = target_scaler.fit_transform(y_train)
        
        params['target_train_means'] = target_scaler.mean_.tolist()
        params['target_train_scales'] = target_scaler.scale_.tolist()
        params['input_train_means'] = input_scaler.mean_.tolist()
        params['input_train_scales'] = input_scaler.scale_.tolist()

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
    os.environ['TF_DETERMINISTIC_OPS'] = '1'
    tf.random.set_seed(42)

    model_name = 'model'
    params['model_name'] = model_name
    model_dir = 'miultihead_loss_test'
    params['model_dir'] = model_dir
    os.makedirs(model_dir, exist_ok=True)

    learning_rate_scheduler = tf.keras.callbacks.LearningRateScheduler(Scheduler)

    input_shape = X_train.shape[1:]
    num_units = 384
    num_layers = 12
    l2_strength = 0.01
    quantiles = [0.16, 0.5, 0.84]
    # quantiles = None
    num_quantiles = len(quantiles) if quantiles else 1
    epochs = 50
    batch_size = 4096
    batch_norm = False
    momentum = 0.01
    bn_eps = 0.01
    learning_rate = 1e-4
    use_energy_layer = False
    use_quantile_ordering = True
    use_quantile_width_penalty = True
    num_output_units = num_quantiles
    add_mass_loss = True

    params['input_shape'] = list(input_shape)
    params['num_units'] = num_units
    params['num_layers'] = num_layers
    params['l2_strength'] = l2_strength
    params['quantiles'] = quantiles
    params['num_quantiles'] = num_quantiles
    params['epochs'] = epochs
    params['batch_size'] = batch_size
    params['batch_norm'] = batch_norm
    params['momentum'] = momentum
    params['bn_eps'] = bn_eps
    params['learning_rate'] = learning_rate
    params['use_energy_layer'] = use_energy_layer
    params['use_quantile_ordering'] = use_quantile_ordering
    params['use_quantile_width_penalty'] = use_quantile_width_penalty
    params['num_output_units'] = num_output_units
    params['add_mass_loss'] = add_mass_loss
    
    # individual width penalty rates
    # some errors are overestimated, others are underestimated
    width_penalty_rates = None
    if use_quantile_width_penalty:
        width_penalty_rates = {}
        width_penalty_rates['genHVV_E'] = 0.05
        width_penalty_rates['genHVV_px'] = 0.05
        width_penalty_rates['genHVV_py'] = 0.05
        width_penalty_rates['genHVV_pz'] = 0.05

        width_penalty_rates['genHbb_E'] = 0.065
        width_penalty_rates['genHbb_px'] = 0.05
        width_penalty_rates['genHbb_py'] = 0.05
        width_penalty_rates['genHbb_pz'] = 0.05

    params['width_penalty_rates'] = width_penalty_rates

    # prepare base part of the model
    inputs = tf.keras.layers.Input(shape=input_shape)
    if batch_norm:
        x = tf.keras.layers.Dense(num_units, kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(inputs)
        x = tf.keras.layers.BatchNormalization(momentum=momentum, epsilon=bn_eps)(x)
        x = TrainableSiLU(units=num_units)(x)
        for _ in range(num_layers - 1):
            x = tf.keras.layers.Dense(num_units, kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(x)
            x = tf.keras.layers.BatchNormalization(momentum=momentum, epsilon=bn_eps)(x)
            x = TrainableSiLU(units=num_units)(x)
    else:
        x = tf.keras.layers.Dense(num_units, 
                                  activation=TrainableSiLU(units=num_units), 
                                  kernel_initializer='random_normal', 
                                  bias_initializer='random_normal',
                                  kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(inputs)
        for _ in range(num_layers - 1):
            x = tf.keras.layers.Dense(num_units, 
                                      activation=TrainableSiLU(units=num_units), 
                                      kernel_initializer='random_normal', 
                                      bias_initializer='random_normal',
                                      kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(x)

    # prepeare output heads
    num_targets = len(target_names)
    losses = {}
    loss_weights = {}
    loss_names = {}
    outputs = []
    ys_train = []
    for idx, target_name in enumerate(target_names):
        target_array = y_train[:, idx]
        target_array = np.reshape(target_array, (-1, 1))

        # losses[target_name] = QuantileLoss(quantiles=quantiles, 
        #                                    width_penalty_rate=width_penalty_rates[target_name], 
        #                                    name=f'quantile_loss_{target_name}')
        losses[target_name] = tf.keras.losses.LogCosh()
        loss_names[target_name] = losses[target_name].name

        # repeat each target array num_quantiles times so that each output node has ground truth value
        target_array = np.repeat(target_array, repeats=num_quantiles, axis=-1)
        ys_train.append(target_array)

        # if custom energy calculation is used, skip energy outputs, they'll be inserted separately
        if use_energy_layer and 'E' in target_name:
            continue

        output = None
        if use_quantile_ordering:
            output = tf.keras.layers.Dense(num_output_units)(x)
            output = QuantileOrderingLayer(name=target_name)(output)
        else:
            output = tf.keras.layers.Dense(num_output_units, name=target_name)(x)
        outputs.append(output)

        loss_weights[target_name] = 1.0

    if use_energy_layer:
        # energy layer of H->bb takes outputs of heads producing p3 of H->bb, but not the rest of the model
        # energy layer of H->VV takes outputs of heads producing p3 of H->VV, but not the rest of the model
        # energy layers have to be inserted back into list of outputs for compatibility with order of targets
        hvv_en_idx = target_names.index('genHVV_E')
        hbb_en_idx = target_names.index('genHbb_E')

        hbb_means = None
        hbb_scales = None
        hvv_means = None
        hvv_scales = None
        if standardize:
            hbb_means = target_scaler.mean_[4:]
            hbb_scales = target_scaler.scale_[4:]

            hvv_means = target_scaler.mean_[:4]
            hvv_scales = target_scaler.scale_[:4]

        hbb_energy_layer = EnergyLayer(num_quantiles=num_quantiles, 
                                       normalize=batch_norm, 
                                       means=hbb_means,
                                       scales=hbb_scales,
                                       name='genHbb_E')(tf.concat(outputs[3:], axis=-1))
        hvv_energy_layer = EnergyLayer(num_quantiles=num_quantiles, 
                                       normalize=batch_norm, 
                                       means=hvv_means,
                                       scales=hvv_scales,
                                       name='genHVV_E')(tf.concat(outputs[:3], axis=-1))
        outputs.insert(hvv_en_idx, hvv_energy_layer)
        outputs.insert(hbb_en_idx, hbb_energy_layer)

        loss_weights['genHbb_E'] = 1.0
        loss_weights['genHVV_E'] = 1.0

    if add_mass_loss:
        # mass loss needs outputs from all heads => concatenate heads into list of tensors and append to outputs
        ys_train.append(y_train)
        combined_outputs = tf.keras.layers.Concatenate(name='combined')(outputs)
        outputs.append(combined_outputs)
        losses['combined'] = MultiheadLoss()
        loss_weights['combined'] = 1.0

    params['loss_names'] = loss_names
    params['loss_weights'] = loss_weights

    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    model.compile(loss=losses,
                  loss_weights=loss_weights, 
                  optimizer=tf.keras.optimizers.Adam(learning_rate))

    history = model.fit(X_train,
                        ys_train,
                        validation_split=0.2,
                        verbose=1,
                        batch_size=batch_size,
                        epochs=epochs,
                        shuffle=True,
                        callbacks=[learning_rate_scheduler])

    tf.keras.utils.plot_model(model, os.path.join(model_dir, f"summary_{model.name}.pdf"), show_shapes=True)
    model.save(os.path.join(model_dir, f"{model_name}.keras"))

    # save training parameters
    with open(os.path.join(model_dir, 'params.yaml'), 'w') as outfile:
        yaml.dump(params, outfile, sort_keys=False, default_flow_style=False)

    plotting_dir = os.path.join(model_dir, 'plots')
    os.makedirs(plotting_dir, exist_ok=True)
    history_keys = list(history.history.keys())
    drawable_metrics = [key for key in history_keys if 'loss' in key and 'val' not in key]
    for metric in drawable_metrics:
        PlotMetric(history, model_name, metric, plotting_dir=plotting_dir)    

    # form testing dataloader
    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)

    X_test, _, y_test, _ = dataloader.Get(TestSelection, 2, 1, 800, 1)
    if standardize:
        X_test = input_scaler.transform(X_test)
    ys_pred = model.predict(X_test)

    # transform predicted data back to normal
    if standardize:
        # each item of ys_pred is a array of shape (n, 3) for 3 quantiles
        # each quantile for given variable should be rescaled back by the same scale
            for pred_idx, arr in enumerate(ys_pred):
                if quantiles:
                    for q_idx in range(num_quantiles):
                        arr[:, q_idx] *= target_scaler.scale_[pred_idx]
                        arr[:, q_idx] += target_scaler.mean_[pred_idx]
                else:
                    arr *= target_scaler.scale_[pred_idx]
                    arr += target_scaler.mean_[pred_idx]

    assert len(ys_pred) == len(target_names), f"mismatch between number of predicted values and ground truth values: ({len(ys_pred)} vs {len(target_names)})"

    # ys_pred is a list of np arrays of shape (num_events, num_quantiles) if quantile regression
    # or (num_events, 1) if normal regression
    pred_dict = {}
    for i, name in enumerate(target_names):
        pred_array = ys_pred[i]
        if quantiles:
            for q_idx, quantile in enumerate(quantiles):
                q_pred_array = pred_array[:, q_idx]
                q = int(100*quantile)
                pred_dict[f'{name}_q{q}'] = q_pred_array
        else:
            pred_dict[name] = pred_array.reshape(-1)

    pred_df = pd.DataFrame.from_dict(pred_dict)

    if quantiles:
        for quantile in quantiles:
            q = int(100*quantile)

            for i, col in enumerate(target_names):
                ground_truth = y_test[:, i]
                output = pred_df[f'{col}_q{q}']
                bin_left = np.min([np.quantile(ground_truth, 0.05), np.quantile(output, 0.05)]) - 1.0
                bin_right = np.max([np.quantile(ground_truth, 0.95), np.quantile(output, 0.95)]) + 1.0
                bins = np.linspace(bin_left, bin_right, 100)
                PlotCompare2D(ground_truth, output, col, bins=bins, quantile=q, plotting_dir=plotting_dir)

            X_E_pred = pred_df[f"genHbb_E_q{q}"] + pred_df[f"genHVV_E_q{q}"]
            X_px_pred = pred_df[f"genHbb_px_q{q}"] + pred_df[f"genHVV_px_q{q}"]
            X_py_pred = pred_df[f"genHbb_py_q{q}"] + pred_df[f"genHVV_py_q{q}"]
            X_pz_pred = pred_df[f"genHbb_pz_q{q}"] + pred_df[f"genHVV_pz_q{q}"]

            X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
            X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
            print(f"quantile {q}: positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

            X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
            PlotHist(data=X_mass_pred, 
                     bins=np.linspace(0, 2000, 100),
                     title="Predicted X->HH mass",
                     ylabel='Count',
                     xlabel='Mass, [GeV]',
                     plotting_dir=plotting_dir,
                     peak=True,
                     width=True,
                     count=True,
                     file_name=f'pred_mass_q{q}')
    else:
        for i, col in enumerate(target_names):
                ground_truth = y_test[:, i]
                bin_left = np.min([np.quantile(ground_truth, 0.05), np.quantile(pred_df[col], 0.05)]) - 1.0
                bin_right = np.max([np.quantile(ground_truth, 0.95), np.quantile(pred_df[col], 0.95)]) + 1.0
                bins = np.linspace(bin_left, bin_right, 100)
                PlotCompare2D(ground_truth, pred_df[col], col, bins=bins, plotting_dir=plotting_dir)

        X_E_pred = pred_df["genHbb_E"] + pred_df["genHVV_E"]
        X_px_pred = pred_df["genHbb_px"] + pred_df["genHVV_px"]
        X_py_pred = pred_df["genHbb_py"] + pred_df["genHVV_py"]
        X_pz_pred = pred_df["genHbb_pz"] + pred_df["genHVV_pz"]

        X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
        X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
        print(f"positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

        X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
        PlotHist(data=X_mass_pred, 
                bins=np.linspace(0, 2000, 100),
                title="Predicted X->HH mass",
                ylabel='Count',
                xlabel='Mass, [GeV]',
                plotting_dir=plotting_dir,
                peak=True,
                width=True,
                count=True,
                file_name='pred_mass')

    # plot errors
    if num_quantiles > 1:
        for i, name in enumerate(target_names):
            up = pred_df[f'{name}_q84']
            down = pred_df[f'{name}_q16']
            ground_truth = y_test[:, i]
            err = up - down
            pred = pred_df[f'{name}_q50']
            correct_predictions = np.abs(pred - ground_truth) < err
            proba = len(ground_truth[correct_predictions])/len(ground_truth)
            print(f'Target {name}:')
            print(f'\tprediction coverage probability: {proba:.2f}')
            print(f'\tquantile crossing fraction: {len(err[err < 0.0])/len(err):.2f}')

            bin_left = np.quantile(err, 0.05) - 1.0
            bin_right = np.quantile(err, 0.95) + 1.0
            bins = np.linspace(bin_left, bin_right, 50)
            PlotHist(data=err, 
                    bins=bins,
                    title=f'Predicted {name} error',
                    ylabel='Count',
                    xlabel='Error, [GeV]',
                    plotting_dir=plotting_dir,
                    file_name=f'pred_error_{name}')



if __name__ == '__main__':
    main()