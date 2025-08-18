import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import yaml
import matplotlib.pyplot as plt

from Dataloader import Dataloader
from LossUtils import QuantileLoss
from LayerUtils import Encoder

from MiscUtils import ground_truth_map, Scheduler, pretty_vars, objects
from PlotUtils import PlotMetric, PlotHist, PlotCompare2D, PlotCovarMtrx

def main():
    params = {}

    # files = []
    # input_files = 'files_Run3_2022.txt'
    # with open(input_files, 'r') as file_cfg:
    #     files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    # dataloader.Load(files)
    dataloader.Load('../nano_0.root')
    input_objects = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                          as_df=False,
                                          unravel=True,
                                          selector=lambda df: df['event'] % 2 == 0)

    print('Inputs:')
    for obj, (names, arr) in inputs.items():
        print(f'\t{obj}: {arr.shape}')

    target_objects = dataloader.GetObjData(objects=['genHbb', 'genHVV'], 
                                           as_df=False,
                                           unravel=True,
                                           selector=lambda df: df['event'] % 2 == 0)

    print('Targets:')
    for obj, (names, arr) in targets.items():
        print(f'\t{obj}: {arr.shape}')

    # prepare model inputs
    lep_inputs = [inputs['lep1'][1], inputs['lep2'][1]]
    met_inputs = [inputs['met'][1]]
    jet_inputs = [inputs[name][1] for name in inputs.keys() if 'centralJet' in name]
    fatjet_inputs = [inputs[name][1] for name in inputs.keys() if 'SelectedFatJet' in name]
    inputs = lep_inputs + met_inputs + jet_inputs + fatjet_inputs

    # prepare input layer for the model
    # each particle will have its own input layer
    # then all layers will be concatenated to take list of particles as model input
    lep_input_layer = tf.keras.layers.Input(shape=lep_inputs[0].shape[1:])
    met_input_layer = tf.keras.layers.Input(shape=met_inputs[0].shape[1:])
    jet_input_layer = tf.keras.layers.Input(shape=jet_inputs[0].shape[1:])
    fatjet_input_layer = tf.keras.layers.Input(shape=fatjet_inputs[0].shape[1:])
    model_input_layer = tf.layers.layer.Concatenate([lep_input_layer]*2 + [met_input_layer] + [jet_input_layer]*10 + [fatjet_input_layer]*3)

    input_mapping = [1]*2 + [2] + [3]*10 + [4]*3
    encoder = Encoder(num_encoder_layers=6, d_model=64, num_heads=8, dff=256, input_mapping=input_mapping)
    

if __name__ == '__main__':
    main()