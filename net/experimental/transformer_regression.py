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

    quantiles = [0.16, 0.5, 0.84]
    num_quantiles = len(quantiles) if quantiles else 1
    epochs = 20
    batch_size = 256
    learning_rate = 3e-4
    use_energy_layer = True
    use_quantile_ordering = False
    use_quantile_width_penalty = True
    num_output_units = num_quantiles
    add_mass_loss = True
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

    # files = []
    # input_files = 'files_Run3_2022.txt'
    # with open(input_files, 'r') as file_cfg:
    #     files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    # dataloader.Load(files)
    dataloader.Load('../train_data/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M_800/nano_0.root')
    input_objects = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                          as_df=False,
                                          unravel=True,
                                          selector=lambda df: df['event'] % 2 == 0)

    print('Inputs:')
    for obj, (names, arr) in input_objects.items():
        print(f'\t{obj}: {arr.shape}')

    target_objects = dataloader.GetObjData(objects=['genHbb', 'genHVV'], 
                                           as_df=False,
                                           unravel=True,
                                           selector=lambda df: df['event'] % 2 == 0)

    # prepare ground truth values
    print('Targets:')
    target_names = [] 
    for obj, (names, arr) in target_objects.items():
        target_names.extend(names)
        print(f'\t{obj}: {arr.shape}')
    print(target_names)

    targets = np.hstack([target_objects['genHbb'][1], target_objects['genHVV'][1]])
    print(targets.shape)

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
    encoder = Encoder(num_encoder_layers=8, d_model=64, num_heads=8, dff=256, input_mapping=input_mapping)
    
    # base transformer mode: inputs -> concatenate -> encoder -> flatten
    x = encoder(model_input_layers)
    x = tf.keras.layers.Flatten()(x)

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

        losses[target_name] = QuantileLoss(quantiles=quantiles, 
                                           width_penalty_rate=width_penalty_rates[target_name] if use_quantile_width_penalty else None, 
                                           name=f'quantile_loss_{target_name}')
        # losses[target_name] = tf.keras.losses.LogCosh()
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
    print(model.summary())
    tf.keras.utils.plot_model(model, 'summary_transformer.pdf', show_shapes=True)

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

if __name__ == '__main__':
    main()