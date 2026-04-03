import tensorflow as tf

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print(f"Enabled memory growth for {len(gpus)} GPU(s).")
    except RuntimeError as e:
        print(e)

import numpy as np
import os
import re
import yaml
import argparse

from typing import List
from PlotUtils import PlotMetric
from MiscUtils import MemoryMonitor, load_file, nearest_pow2, map_input_files, threshold_cleaning, make_dataset, quantile_cleaning
from sklearn.utils import class_weight
from LayerUtils import ResidualBlock

tf.random.set_seed(42)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'
os.environ['TF_GPU_ALLOCATOR'] = 'cuda_malloc_async'
np.random.seed(42)

def main():
    mm = MemoryMonitor()
    mm.print_memory_usage(msg="Training DNN for event classification per mass point")

    parser = argparse.ArgumentParser(
        description="Train event classification model according to provided cofiguration",
    )
    parser.add_argument(
        '--config',
        type=str,
        help="Model configuration file",
        required=True
    )
    args = parser.parse_args()
    cfg_path = args.config
    print(f"Will use configuration: {cfg_path}")

    with open (cfg_path, 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    weight_file_pattern = re.compile(r"nParity(\d)_Merged_weight.root")
    event_file_pattern = re.compile(r"nParity(\d)_Merged.root")
    train_file_dir = cfg["input_directory"]
    parity_file_map = map_input_files(directory=train_file_dir,
                                      weight_file_pattern=weight_file_pattern,
                                      event_file_pattern=event_file_pattern)

    hyperparameters = cfg['hyperparameters']

    train_parity = hyperparameters['training_parity']
    train_event_file_path, train_weight_file_path = parity_file_map[train_parity]
    print(f"Data file for trainig: {train_event_file_path}")
    print(f"Weight file for trainig: {train_weight_file_path}")

    use_extra_variables = hyperparameters['use_extra_branches']

    branches_to_load = cfg['base_branches']
    if use_extra_variables:
        branches_to_load.extend(cfg['extra_branches'])

    X_train = load_file(tree_name="Events", 
                        file_path=train_event_file_path, 
                        list_of_branches=branches_to_load,
                        convert_to_numpy=True)
    
    cleaning_cfg = cfg["data_cleaning"]
    mm.print_memory_usage(msg=f"Before cleaning: {X_train.shape}")
    cleaning_type = cleaning_cfg["type"]
    high = cleaning_cfg["high"]
    low = cleaning_cfg["low"]
    if cleaning_type == "threshold":
        clean_mask_train, X_train = threshold_cleaning(
            data=X_train,
            pos_thrsh=low,
            neg_thrsh=high
        )
    elif cleaning_type == "quantile":
        clean_mask_train, X_train = quantile_cleaning(
            data=X_train,
            q_low=low,
            q_high=high
        )
    else:
        raise RuntimeError(f"Unsupported data cleaning type {cleaning_type}.")
    mm.print_memory_usage(msg=f"After cleaning: {X_train.shape}")

    y_train = load_file(tree_name="weight_tree", 
                        file_path=train_weight_file_path, 
                        list_of_branches=["class_target"],
                        convert_to_numpy=True)
    y_train = y_train[clean_mask_train]
    y_train = y_train.reshape(-1)

    multiclass = cfg['multiclass']
    if not multiclass:
        # swap classes: 1=>signal, 0=>background
        y_train = 1 - y_train
    
    apply_class_weights = cfg['class_weights']
    if apply_class_weights:
        class_weights = class_weight.compute_class_weight(
            'balanced',
            classes=np.unique(y_train),
            y=y_train
        )
        class_weight_dict = dict(enumerate(class_weights))
        print(f"Applying weights: {class_weight_dict}")
    else:
        class_weight_dict = None

    mm.print_memory_usage(msg=f"Loaded {X_train.shape[1]} variables for {X_train.shape[0]} events for training set")

    # shuffle train set
    indices = np.random.permutation(len(X_train))
    X_train = X_train[indices]
    y_train = y_train[indices]

    val_parity = hyperparameters['validation_parity']
    val_event_file_path, val_weight_file_path = parity_file_map[val_parity]
    X_val = load_file(tree_name="Events", 
                      file_path=val_event_file_path, 
                      list_of_branches=branches_to_load,
                      convert_to_numpy=True)
    
    y_val = load_file(tree_name="weight_tree", 
                     file_path=val_weight_file_path, 
                     list_of_branches=["class_target"],
                     convert_to_numpy=True)
    y_val = y_val.reshape(-1)

    if not multiclass:
        y_val = 1 - y_val

    mm.print_memory_usage(msg=f"Loaded {X_val.shape[1]} variables for {X_val.shape[0]} events for validation set")

    # precompute mean and variance
    train_mean = np.mean(X_train, axis=0)
    train_variance = np.var(X_train, axis=0)

    assert train_mean.shape == train_variance.shape and train_mean.shape[0] == X_train.shape[1]
    
    # load X_mass parameter
    X_mass = load_file(
        tree_name="Events", 
        file_path=train_event_file_path, 
        list_of_branches=["X_mass"],
        convert_to_numpy=True
    )
    masspoints = np.unique(X_mass[X_mass > 0.0])
    
    # create full dataset: contains all mass of signal + background
    full_dataset = tf.data.Dataset.from_tensor_slices((X_train, y_train, X_mass))

    batch_size = hyperparameters["batch_size"]

    # loop over all masspoints and train model for each mass
    for mp in masspoints:
        filtered_dataset = full_dataset.filter(lambda x, y, mx: mx == mp)
    
        # Remove the parameter column for training
        mp_ds = filtered_dataset.map(lambda x, y, mp: (x, y))
        # mp_ds = mp_ds.shuffle(10000).batch(batch_size).prefetch(tf.data.AUTOTUNE)

        sig_mp_ds = mp_ds.filter(lambda x, y: y == 1).shuffle(10000)
        bkg_mp_ds = mp_ds.filter(lambda x, y: y == 0).shuffle(100000)

        tf.keras.backend.clear_session()

if __name__ == "__main__":
    main()