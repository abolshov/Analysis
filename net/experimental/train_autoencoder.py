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

from MiscUtils import MemoryMonitor, load_file, map_input_files, clean_extreme_values
from AnomalyDetection import Autoencoder
from PlotUtils import PlotMetric

tf.random.set_seed(42)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'

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

    X_train = load_file(tree_name="Events", 
                        file_path=train_event_file_path, 
                        list_of_branches=branches_to_load,
                        convert_to_numpy=True)
    train_mask, X_train = clean_extreme_values(data=X_train)
    
    y_train = load_file(tree_name="weight_tree", 
                        file_path=train_weight_file_path, 
                        list_of_branches=["class_target"],
                        convert_to_numpy=True)
    y_train = y_train[train_mask]
    
    y_train = 1 - y_train
    sig_mask_train = y_train == 1
    bkg_mask_train = ~sig_mask_train
    
    sig_mask_train = sig_mask_train.reshape(-1)
    bkg_mask_train = bkg_mask_train.reshape(-1)

    mm.print_memory_usage(msg="After loading train set")

    X_bkg_train = X_train[bkg_mask_train]
    # X_sig_train = X_train[sig_mask]

    input_dim = X_bkg_train.shape[-1]
    input_shape = X_bkg_train.shape[1:]

    # precompute mean and variance
    # train_mean = np.mean(X_bkg_train, axis=0)
    # train_variance = np.var(X_bkg_train, axis=0)
    train_mins = np.min(X_bkg_train, axis=0)
    train_maxs = np.max(X_bkg_train, axis=0)

    is_bad_feature = np.isclose(train_maxs, train_mins)
    num_bad_features = np.sum(is_bad_feature)
    if num_bad_features > 0:
        bad_feature_idxs = np.flatnonzero(is_bad_feature)
        bad_feature_names = [branches_to_load[idx] for idx in bad_feature_idxs]
        raise RuntimeError(f"Have {num_bad_features} features with coinciding min and max: {bad_feature_names}")

    X_bkg_train = (X_bkg_train - train_mins)/(train_maxs - train_mins)

    inputs = tf.keras.Input(shape=(X_bkg_train.shape[1],))
    x = tf.keras.layers.Identity()(inputs)
    for dim in [64, 32, 16, 8]:
        x = tf.keras.layers.Dense(dim, activation='relu')(x)

    for dim in [16, 32, 64]:
        x = tf.keras.layers.Dense(dim, activation='relu')(x)

    outputs = tf.keras.layers.Dense(input_dim, activation='sigmoid')(x)

    anomaly_detector = tf.keras.models.Model(inputs=inputs, 
                                             outputs=outputs, 
                                             name="autoencoder")
    
    anomaly_detector.compile(
        loss=tf.keras.losses.MeanSquaredError(),
        optimizer=tf.keras.optimizers.Adam(3e-4)    
    )

    print(anomaly_detector.summary())

    val_parity = 1
    val_event_file_path, val_weight_file_path = parity_file_map[val_parity]
    X_val = load_file(tree_name="Events", 
                      file_path=val_event_file_path, 
                      list_of_branches=branches_to_load,
                      convert_to_numpy=True)
    val_mask, X_val = clean_extreme_values(data=X_val)

    y_val = load_file(tree_name="weight_tree", 
                      file_path=val_weight_file_path, 
                      list_of_branches=["class_target"],
                      convert_to_numpy=True)
    y_val = y_val[val_mask]
    mm.print_memory_usage(msg="After loading validation set")

    bkg_mask_val = y_val == 0
    bkg_mask_val = bkg_mask_val.reshape(-1)
    X_val = X_val[bkg_mask_val]

    val_mins = np.min(X_val, axis=0)
    val_maxs = np.max(X_val, axis=0)
    X_val = (X_val - val_mins)/(val_maxs - val_mins)

    callbacks = [
        tf.keras.callbacks.EarlyStopping(patience=10, restore_best_weights=True)
    ]

    train_history = anomaly_detector.fit(
        X_bkg_train,
        X_bkg_train,
        validation_data=(X_val, X_val),
        epochs=200,
        batch_size=2048,
        callbacks=callbacks,
        verbose=0,
        shuffle=True
    )

    model_dir = os.path.join("/home/artem/Desktop/CMS/models/anomaly_detection/", anomaly_detector.name)
    PlotMetric(train_history, anomaly_detector.name, "loss", plotting_dir=model_dir)
    anomaly_detector.save(os.path.join(model_dir, f"{anomaly_detector.name}.keras"))

    # reconstructions = anomaly_detector.predict(X_val, batch_size=128)
    # val_loss = tf.keras.losses.MSE(reconstructions, X_val)

    # plt.hist(val_loss[None, :], bins=50)
    # plt.xlabel("Validation loss")
    # plt.ylabel("No of examples")
    # plt.savefig(os.path.join(model_dir, "anomaly_score.pdf"), format="pdf")

if __name__ == "__main__":
    main()