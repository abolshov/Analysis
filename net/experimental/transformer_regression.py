import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import yaml
import matplotlib.pyplot as plt

from Dataloader import Dataloader
from LossUtils import QuantileLoss, MultiheadLoss
from LayerUtils import Encoder, EnergyLayer

from MiscUtils import ground_truth_map, Scheduler, pretty_vars, objects
from PlotUtils import PlotMetric, PlotHist, PlotCompare2D, PlotCovarMtrx

def main():
    params = {}

    model_dir = 'transformer_v1'
    model_name = 'transformer_v1'
    os.makedirs(model_dir, exist_ok=True)

    quantiles = None
    num_quantiles = len(quantiles) if quantiles else 1
    epochs = 30
    batch_size = 512
    learning_rate = 1e-3
    use_energy_layer = True
    use_quantile_ordering = False
    use_quantile_width_penalty = False
    num_output_units = num_quantiles
    add_mass_loss = False
    standardize = False
    batch_norm = False

    width_penalty_rates = {}
    width_penalty_rates['genHVV_E'] = 0.035
    width_penalty_rates['genHVV_px'] = 0.035
    width_penalty_rates['genHVV_py'] = 0.035
    width_penalty_rates['genHVV_pz'] =  0.02

    width_penalty_rates['genHbb_E'] = 0.035
    width_penalty_rates['genHbb_px'] = 0.035
    width_penalty_rates['genHbb_py'] = 0.035
    width_penalty_rates['genHbb_pz'] =  0.025

    files = []
    input_files = 'files_Run3_2022.txt'
    with open(input_files, 'r') as file_cfg:
        files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(files)
    # dataloader.Load('../train_data/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M_800/nano_0.root')
    input_objects = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                          as_df=False,
                                          unravel=True,
                                          selector=lambda df: df['event'] % 2 == 0)

    target_objects = dataloader.GetObjData(objects=['genHbb', 'genHVV'], 
                                           as_df=False,
                                           unravel=True,
                                           selector=lambda df: df['event'] % 2 == 0)

    # prepare ground truth values
    target_names = [] 
    for obj, (names, arr) in target_objects.items():
        target_names.extend(names)
    targets = np.hstack([target_objects['genHbb'][1], target_objects['genHVV'][1]])

    # prepare model inputs
    lep_inputs = [input_objects['lep1'][1], input_objects['lep2'][1]]
    met_inputs = [input_objects['met'][1]]
    jet_inputs = [input_objects[name][1] for name in input_objects.keys() if 'centralJet' in name]
    fatjet_inputs = [input_objects[name][1] for name in input_objects.keys() if 'SelectedFatJet' in name]
    X_train = lep_inputs + met_inputs + jet_inputs + fatjet_inputs

    # prepare input layer for the model
    # each particle will have its own input layer
    # then all layers will be concatenated to take list of particles as model input
    lep_input_layer = tf.keras.layers.Input(shape=lep_inputs[0].shape[1:], name='leptons')
    met_input_layer = tf.keras.layers.Input(shape=met_inputs[0].shape[1:], name='MET')
    jet_input_layer = tf.keras.layers.Input(shape=jet_inputs[0].shape[1:], name='AK4Jets')
    fatjet_input_layer = tf.keras.layers.Input(shape=fatjet_inputs[0].shape[1:], name='AK8Jets')

    model_input_layers = [lep_input_layer]*2 + [met_input_layer] + [jet_input_layer]*10 + [fatjet_input_layer]*3
    model_input_layer = tf.keras.layers.Concatenate(axis=1)(model_input_layers)

    input_mapping = [1]*2 + [2] + [3]*10 + [4]*3
    encoder = Encoder(num_encoder_layers=6, d_model=32, num_heads=8, dff=128, input_mapping=input_mapping)
    
    # base transformer mode: inputs -> concatenate -> encoder -> flatten
    x = encoder(model_input_layers)
    x = tf.keras.layers.Flatten()(x)
    x = tf.keras.layers.Dense(units=256, activation='silu')(x)
    x = tf.keras.layers.Dense(units=256, activation='silu')(x)

    # attach heads to transformer
    num_targets = len(target_names)
    losses = {}
    loss_weights = {}
    loss_names = {}
    outputs = []
    ys_train = []
    y_train = targets
    for idx, target_name in enumerate(target_names):
        target_array = y_train[:, idx]
        target_array = np.reshape(target_array, (-1, 1))

        # losses[target_name] = QuantileLoss(quantiles=quantiles, 
        #                                    width_penalty_rate=width_penalty_rates[target_name] if use_quantile_width_penalty else None, 
        #                                    name=f'quantile_loss_{target_name}')
        losses[target_name] = tf.keras.losses.LogCosh()
        loss_names[target_name] = losses[target_name].name

        # repeat each target array num_quantiles times so that each output node has ground truth value
        if num_quantiles > 1:
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

        # loss_weights[target_name] = 0.675
        loss_weights[target_name] = 0.75
        # loss_weights[target_name] = 1.0

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

        hbb_energy_input = tf.keras.layers.Concatenate(name='hbb_energy_input')(outputs[:3])
        hbb_energy_layer = EnergyLayer(num_quantiles=num_quantiles, 
                                       normalize=batch_norm, 
                                       means=hbb_means,
                                       scales=hbb_scales,
                                       name='genHbb_E')(hbb_energy_input)
        
        hvv_energy_input = tf.keras.layers.Concatenate(name='hvv_energy_input')(outputs[3:])
        hvv_energy_layer = EnergyLayer(num_quantiles=num_quantiles, 
                                       normalize=batch_norm, 
                                       means=hvv_means,
                                       scales=hvv_scales,
                                       name='genHVV_E')(hvv_energy_input)
        
        outputs.insert(hvv_en_idx, hvv_energy_layer)
        outputs.insert(hbb_en_idx, hbb_energy_layer)

        loss_weights['genHbb_E'] = 0.25
        loss_weights['genHVV_E'] = 0.25

    if add_mass_loss:
        # mass loss needs outputs from all heads => concatenate heads into list of tensors and append to outputs
        ys_train.append(y_train)
        combined_outputs = tf.keras.layers.Concatenate(name='combined')(outputs)
        outputs.append(combined_outputs)
        
        means = None
        scales = None
        if batch_norm and standardize:
            means = target_scaler.mean_
            scales = target_scaler.scale_

        losses['combined'] = MultiheadLoss(num_quantiles=num_quantiles, means=means, scales=scales)
        loss_weights['combined'] = 0.001

    model = tf.keras.Model(inputs=model_input_layers, outputs=outputs, name='transformer')
    tf.keras.utils.plot_model(model, os.path.join(model_dir, f"summary_{model.name}.pdf"), show_shapes=True)

    model.compile(loss=losses,
                  loss_weights=loss_weights, 
                  optimizer=tf.keras.optimizers.Adam(learning_rate))

    history = model.fit(X_train,
                        ys_train,
                        validation_split=0.2,
                        verbose=1,
                        batch_size=batch_size,
                        epochs=epochs,
                        shuffle=True)

    model.save(os.path.join(model_dir, f"{model_name}.keras"))

    plotting_dir = os.path.join(model_dir, 'plots')
    os.makedirs(plotting_dir, exist_ok=True)
    history_keys = list(history.history.keys())
    drawable_metrics = [key for key in history_keys if 'loss' in key and 'val' not in key]
    for metric in drawable_metrics:
        PlotMetric(history, model_name, metric, plotting_dir=plotting_dir)

    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)


    input_objects = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                          as_df=False,
                                          unravel=True,
                                          selector=lambda df: np.logical_and(df['event'] % 2 == 1, df['X_mass'] == 800))

    lep_inputs = [input_objects['lep1'][1], input_objects['lep2'][1]]
    met_inputs = [input_objects['met'][1]]
    jet_inputs = [input_objects[name][1] for name in input_objects.keys() if 'centralJet' in name]
    fatjet_inputs = [input_objects[name][1] for name in input_objects.keys() if 'SelectedFatJet' in name]
    X_test = lep_inputs + met_inputs + jet_inputs + fatjet_inputs
    ys_pred = model.predict(X_test)

    target_objects = dataloader.GetObjData(objects=['genHbb', 'genHVV'], 
                                           as_df=False,
                                           unravel=True,
                                           selector=lambda df: np.logical_and(df['event'] % 2 == 1, df['X_mass'] == 800))
    y_test = np.hstack([target_objects['genHbb'][1], target_objects['genHVV'][1]])

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

    is_quantile = quantiles is not None
    has_central_quantile = is_quantile and 0.5 in quantiles
    if is_quantile and has_central_quantile:
        q = int(100*0.5)

        for i, col in enumerate(target_names):
            ground_truth = y_test[:, i]
            output = pred_df[f'{col}_q{q}']
            bin_left = np.min([np.quantile(ground_truth, 0.01), np.quantile(output, 0.01)]) - 1.0
            bin_right = np.max([np.quantile(ground_truth, 0.99), np.quantile(output, 0.99)]) + 1.0
            bins = np.linspace(bin_left, bin_right, 100)
            
            var = col.split('_')[-1]
            obj = col.split('_')[0]
            PlotCompare2D(ground_truth, 
                          output, 
                          col.split('_')[-1], 
                          col.split('_')[0],
                          bins=bins, 
                          title=f'{objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} comparison',
                          xlabel=f'True {ground_truth_map[col]}',
                          ylabel=f'Predicted {ground_truth_map[col]}',
                          quantile=q, 
                          plotting_dir=plotting_dir)

        X_E_pred = pred_df[f"genHbb_E_q{q}"] + pred_df[f"genHVV_E_q{q}"]
        X_px_pred = pred_df[f"genHbb_px_q{q}"] + pred_df[f"genHVV_px_q{q}"]
        X_py_pred = pred_df[f"genHbb_py_q{q}"] + pred_df[f"genHVV_py_q{q}"]
        X_pz_pred = pred_df[f"genHbb_pz_q{q}"] + pred_df[f"genHVV_pz_q{q}"]

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
    elif not is_quantile:
        # not quantiles <=> model contains single output for central value without quantile loss
        for i, col in enumerate(target_names):
                ground_truth = y_test[:, i]
                bin_left = np.min([np.quantile(ground_truth, 0.01), np.quantile(pred_df[col], 0.01)]) - 1.0
                bin_right = np.max([np.quantile(ground_truth, 0.99), np.quantile(pred_df[col], 0.99)]) + 1.0
                bins = np.linspace(bin_left, bin_right, 100)
                
                var = col.split('_')[-1]
                obj = col.split('_')[0]
                PlotCompare2D(ground_truth, 
                              pred_df[col], 
                              var, 
                              obj,
                              bins=bins, 
                              title=f'{objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} comparison',
                              xlabel=f'True {ground_truth_map[col]}',
                              ylabel=f'Predicted {ground_truth_map[col]}',
                              plotting_dir=plotting_dir)

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

            bin_left = np.min(err) - 1.0
            bin_right = np.max(err) + 1.0
            bins = np.linspace(bin_left, bin_right, 50)
            PlotHist(data=err, 
                     bins=bins,
                     title=f'Predicted {name} error',
                     ylabel='Count',
                     xlabel='Error, [GeV]',
                     plotting_dir=plotting_dir,
                     file_name=f'pred_error_{name}')

            # plot up vs down errors to see if they are symmetric
            var = name.split('_')[-1]
            obj = name.split('_')[0]
            xlabel = f'Up {objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} error'
            ylabel = f'Down {objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} error'
            title = f'{objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} error comparison'
            
            bin_left = np.min([up, down]) - 1.0
            bin_right = np.max([up, down]) + 1.0
            bins = np.linspace(bin_left, bin_right, 100)
            PlotCompare2D(up, 
                          down, 
                          f'{var}_error', 
                          obj,
                          title=title,
                          bins=bins, 
                          xlabel=xlabel,
                          ylabel=ylabel,
                          plotting_dir=plotting_dir)

if __name__ == '__main__':
    main()