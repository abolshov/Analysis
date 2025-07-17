import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


from Dataloader import Dataloader
from NetUtils import TrainableSiLU, EpochLossUpdater
from LossUtils import QuantileLoss
from LayerUtils import EnergyLayer, QuantileOrderingLayer


from MiscUtils import *
from PlotUtils import PlotMetric, PlotHist


def PlotCompare2D(target, output, quantity, quantile, plotting_dir=None):
    plt.grid(False)
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist2d(target, output, bins=bins)
    var = quantity.split('_')[1]
    plt.title(f'{var} comparison')
    plt.ylabel(f'predicted {ground_truth_map[quantity]}')
    plt.xlabel(f'true {ground_truth_map[quantity]}')
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f"cmp2d_{quantity}_{quantile}.pdf"), bbox_inches='tight')
    else:
        plt.savefig(f"cmp2d_{quantity}_{quantile}.pdf", bbox_inches='tight')
    plt.clf()


def main():
    files = []
    with open('dl_train_files.txt', 'r') as file_cfg:
        files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(files)
    X_train, input_names, y_train, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)

    standardize = False
    input_scaler = None
    target_scaler = None
    if standardize:
        input_scaler = StandardScaler()
        X_train = input_scaler.fit_transform(X_train)
        target_scaler = StandardScaler()
        y_train = target_scaler.fit_transform(y_train)

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
    os.environ['TF_DETERMINISTIC_OPS'] = '1'
    tf.random.set_seed(42)

    model_name = "model_hh_dl_train_silu_qloss_width_pen"
    model_dir = model_name
    os.makedirs(model_dir, exist_ok=True)

    learning_rate_scheduler = tf.keras.callbacks.LearningRateScheduler(Scheduler)

    input_shape = X_train.shape[1:]
    num_units = 384
    num_layers = 12
    l2_strength = 0.01
    quantiles = [0.16, 0.5, 0.84]
    num_quantiles = len(quantiles)
    epochs = 50
    batch_size = 512
    batch_norm = False
    momentum = 0.01
    use_energy_layer = False
    use_quantile_ordering = True
    width_penalty_rate = 0.02

    # prepare base part of the model
    inputs = tf.keras.layers.Input(shape=input_shape)
    if batch_norm:
        x = tf.keras.layers.Dense(num_units, kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(inputs)
        x = tf.keras.layers.BatchNormalization(momentum=momentum)(x)
        x = TrainableSiLU(units=num_units)(x)
        for _ in range(num_layers - 1):
            x = tf.keras.layers.Dense(num_units, kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(x)
            x = tf.keras.layers.BatchNormalization(momentum=momentum)(x)
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
    outputs = []
    ys_train = []
    for idx, target_name in enumerate(target_names):
        target_array = y_train[:, idx]
        target_array = np.reshape(target_array, (-1, 1))

        losses[target_name] = QuantileLoss(quantiles=quantiles, 
                                           width_penalty_rate=width_penalty_rate, 
                                           name=f'loss_{target_name}')

        # repeat each target array num_quantiles times so that each output node has ground truth value
        target_array = np.repeat(target_array, repeats=num_quantiles, axis=-1)
        ys_train.append(target_array)

        # if custom energy calculation is used, skip energy outputs, they'll be inserted separately
        if use_energy_layer and 'E' in target_name:
            continue

        output = None
        if use_quantile_ordering:
            output = tf.keras.layers.Dense(num_quantiles)(x)
            output = QuantileOrderingLayer(name=target_name)(output)
        else:
            output = tf.keras.layers.Dense(num_quantiles, name=target_name)(x)
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

    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    model.compile(loss=losses,
                  loss_weights=loss_weights, 
                  optimizer=tf.keras.optimizers.Adam(3e-4))

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

    history_keys = list(history.history.keys())
    drawable_metrics = [key for key in history_keys if 'loss' in key and 'val' not in key]
    for metric in drawable_metrics:
        PlotMetric(history, model_name, metric, plotting_dir=model_dir)    

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
        # eachitem of ys_pred is a array of shape (n, 3) for 3 quantiles
        # each quantile for given variable should be rescaled back by the same scale
        for pred_idx, arr in enumerate(ys_pred):
            for q_idx in range(num_quantiles):
                arr[:, q_idx] *= target_scaler.scale_[pred_idx]
                arr[:, q_idx] += target_scaler.mean_[pred_idx]

    assert len(ys_pred) == len(target_names), f"mismatch between number of predicted values and ground truth values: ({len(ys_pred)} vs {len(target_names)})"

    # ys_pred is a list of np arrays of shape (num_events, num_quantiles)
    pred_dict = {}
    for i, name in enumerate(target_names):
        pred_array = ys_pred[i]
        for q_idx, quantile in enumerate(quantiles):
            q_pred_array = pred_array[:, q_idx]
            q = int(100*quantile)
            pred_dict[f'{name}_q{q}'] = q_pred_array

    pred_df = pd.DataFrame.from_dict(pred_dict)

    bins = np.linspace(0, 2000, 100)
    for quantile in quantiles:
        q = int(100*quantile)

        for i, col in enumerate(target_names):
            ground_truth = y_test[:, i]
            PlotCompare2D(ground_truth, pred_df[f'{col}_q{q}'], col, q, plotting_dir=model_dir)

        X_E_pred = pred_df[f"genHbb_E_q{q}"] + pred_df[f"genHVV_E_q{q}"]
        X_px_pred = pred_df[f"genHbb_px_q{q}"] + pred_df[f"genHVV_px_q{q}"]
        X_py_pred = pred_df[f"genHbb_py_q{q}"] + pred_df[f"genHVV_py_q{q}"]
        X_pz_pred = pred_df[f"genHbb_pz_q{q}"] + pred_df[f"genHVV_pz_q{q}"]

        X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
        X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
        print(f"quantile {q}: positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

        X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
        PlotHist(data=X_mass_pred, 
                 bins=bins,
                 title="Predicted X->HH mass",
                 ylabel='Count',
                 xlabel='Mass, [GeV]',
                 plotting_dir=model_dir,
                 peak=True,
                 width=True,
                 count=True,
                 file_name=f'pred_mass_q{q}')

    # plot errors
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

        bins = np.linspace(np.min(err) - 10, np.max(err) + 10, 50)
        PlotHist(data=err, 
                 bins=bins,
                 title=f'Predicted {name} error',
                 ylabel='Count',
                 xlabel='Error, [GeV]',
                 plotting_dir=model_dir,
                 file_name=f'pred_error_{name}')

        actual_error = pred - ground_truth
        bins = np.linspace(np.min(actual_error) - 10, np.max(actual_error) + 10, 50)
        PlotHist(data=actual_error, 
                 bins=bins,
                 title=f'Actual {name} error',
                 ylabel='Count',
                 xlabel='Error, [GeV]',
                 plotting_dir=model_dir,
                 file_name=f'error_{name}')


if __name__ == '__main__':
    main()