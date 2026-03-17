import tensorflow as tf

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print(f"Enabled memory growth for {len(gpus)} GPU(s).")
    except RuntimeError as e:
        print(e)

import uproot
import pathlib
import awkward as ak
import numpy as np
import os
import re
import yaml

from MiscUtils import MemoryMonitor, load_file, map_input_files
from AnomalyDetection import Autoencoder

def main():
    mm = MemoryMonitor()
    mm.print_memory_usage(msg="Training autoencoder")

    weight_file_pattern = re.compile(r"nParity(\d)_Merged_weight.root")
    event_file_pattern = re.compile(r"nParity(\d)_Merged.root")
    train_file_dir = "/home/artem/Desktop/CMS/data/DNN/SL/resolved/Run3_2022/Dataset"
    parity_file_map = map_input_files(directory=train_file_dir,
                                      weight_file_pattern=weight_file_pattern,
                                      event_file_pattern=event_file_pattern)

    train_parity = 0
    train_event_file_path, train_weight_file_path = parity_file_map[train_parity]
    use_extra_variables = True

    with open ('classifier_cfg.yaml', 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    branches_to_load = cfg['base_branches']
    if use_extra_variables:
        branches_to_load.extend(cfg['extra_branches'])

    X = load_file(tree_name="Events", 
                  file_path=train_event_file_path, 
                  list_of_branches=branches_to_load,
                  convert_to_numpy=True)
    
    y = load_file(tree_name="weight_tree", 
                  file_path=train_weight_file_path, 
                  list_of_branches=["class_target"],
                  convert_to_numpy=True)
    
    y = 1 - y
    sig_mask = y == 1
    bkg_mask = ~sig_mask
    
    mm.print_memory_usage(msg="After loading train set")

    X_bkg = X[bkg_mask]
    X_sig = X[sig_mask]

    input_dim = X_bkg.shape[-1]
    anomaly_detector = Autoencoder(encoder_dims=[input_dim, 64, 32, 16, 8], 
                                   decoder_dims=[16, 32, 64, input_dim],
                                   hidden_activation="relu",
                                   output_activation=None,
                                   name="autoencoder_anomaly_detector")
    
    anomaly_detector.compile(
        loss=tf.keras.losses.MeanSquaredError(),
        optimizer=tf.keras.optimizers.Adam(3e-4)    
    )

    train_history = anomaly_detector.fit(
        X_bkg,
        X_bkg,
        validation_data=(X_sig, X_sig),
        epochs=50,
        batch_size=2048
    )



if __name__ == "__main__":
    main()