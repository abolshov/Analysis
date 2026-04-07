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

from MiscUtils import MemoryMonitor, load_file, map_input_files, threshold_cleaning, quantile_cleaning
from PlotUtils import PlotMetric

tf.random.set_seed(42)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'

def main():
    mm = MemoryMonitor()
    mm.print_memory_usage(msg="Training autoencoder")

    parser = argparse.ArgumentParser(
        description="Train autoencoder according to provided cofiguration",
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
    parity_file_map = map_input_files(
        directory=train_file_dir,
        weight_file_pattern=weight_file_pattern,
        event_file_pattern=event_file_pattern
    )

    hyperparameters = cfg['hyperparameters']

    train_parity = hyperparameters['training_parity']
    train_event_file_path, train_weight_file_path = parity_file_map[train_parity]
    print(f"Data file for trainig: {train_event_file_path}")
    print(f"Weight file for trainig: {train_weight_file_path}")

    use_extra_variables = hyperparameters['use_extra_branches']
    branches_to_load = cfg['base_branches']
    if use_extra_variables:
        branches_to_load.extend(cfg['extra_branches'])

    X_train = load_file(
        tree_name="Events", 
        file_path=train_event_file_path, 
        list_of_branches=branches_to_load,
        convert_to_numpy=True
    )
    
    cleaning_cfg = cfg["data_cleaning"]
    cleaning_type = cleaning_cfg["type"]
    high = cleaning_cfg["high"]
    low = cleaning_cfg["low"]
    if cleaning_type == "threshold":
        train_mask, X_train = threshold_cleaning(
            data=X_train,
            pos_thrsh=high,
            neg_thrsh=low
        )
    elif cleaning_type == "quantile":
        train_mask, X_train = quantile_cleaning(
            data=X_train,
            q_low=low,
            q_high=high
        )
    else:
        raise RuntimeError(f"Unsupported data cleaning type {cleaning_type}.")
    
    mm.print_memory_usage(msg=f"Cleaned dataset contains {len(X_train)} events")

    y_train = load_file(
        tree_name="weight_tree", 
        file_path=train_weight_file_path, 
        list_of_branches=["class_target"],
        convert_to_numpy=True
    )
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

    # precompute mean and variance
    train_mean = np.mean(X_bkg_train, axis=0)
    train_variance = np.var(X_bkg_train, axis=0)

    # train_mins = np.min(X_bkg_train, axis=0)
    # train_maxs = np.max(X_bkg_train, axis=0)

    # is_bad_feature = np.isclose(train_maxs, train_mins)
    # num_bad_features = np.sum(is_bad_feature)
    # if num_bad_features > 0:
    #     bad_feature_idxs = np.flatnonzero(is_bad_feature)
    #     bad_feature_names = [branches_to_load[idx] for idx in bad_feature_idxs]
    #     raise RuntimeError(f"Have {num_bad_features} features with coinciding min and max: {bad_feature_names}")

    # X_bkg_train = (X_bkg_train - train_mins)/(train_maxs - train_mins)

    encoder_activation = hyperparameters["encoder_activation"]
    decoder_activation = hyperparameters["decoder_activation"]
    input_noise_std = hyperparameters["input_noise_std"]
    encoder_dropout = hyperparameters["encoder_dropout"]
    decoder_dropout = hyperparameters["decoder_dropout"]
    bottleneck_activation = hyperparameters["bottleneck_activation"]
    bottleneck_noise_std = hyperparameters["bottleneck_noise_std"]

    # encoder part
    inputs = tf.keras.Input(shape=(X_bkg_train.shape[1],))
    x = tf.keras.layers.Identity()(inputs)
    if input_noise_std > 0.0:
        x = tf.keras.layers.GaussianNoise(input_noise_std)(x)
    x = tf.keras.layers.Normalization(mean=train_mean, variance=train_variance)(x)
    # [64, 48, 24, 12] ? 
    for dim in [64, 32, 16]:
        x = tf.keras.layers.BatchNormalization()(x)
        x = tf.keras.layers.Dense(dim, activation=encoder_activation)(x)
        if encoder_dropout > 0.0:
            x = tf.keras.layers.Dropout(encoder_dropout)(x)

    # bottleneck part
    x = tf.keras.layers.Dense(8, activation=bottleneck_activation)(x)
    x = tf.keras.layers.BatchNormalization()(x)
    if bottleneck_noise_std > 0.0:
        x = tf.keras.layers.GaussianNoise(bottleneck_noise_std)(x)

    # decoder part
    for dim in [16, 32, 64]:
        x = tf.keras.layers.BatchNormalization()(x)
        x = tf.keras.layers.Dense(dim, activation=decoder_activation)(x)
        if decoder_dropout > 0.0:
            x = tf.keras.layers.Dropout(decoder_dropout)(x)

    outputs = tf.keras.layers.Dense(input_dim)(x)
    output_activation_name = hyperparameters["output_activation"]
    if output_activation_name:
        output_activation_fn = tf.keras.activations.get(output_activation_name)
        outputs = output_activation_fn(outputs)

    model_name = cfg['model_name']
    anomaly_detector = tf.keras.models.Model(
        inputs=inputs, 
        outputs=outputs, 
        name=model_name
    )
    
    loss_cfg = cfg["loss"]
    loss_instance = tf.keras.losses.get(loss_cfg)
    optimizer_cfg = cfg['optimizer']
    optimizer_instance = tf.keras.optimizers.get(optimizer_cfg)
    anomaly_detector.compile(
        loss=loss_instance,
        optimizer=optimizer_instance,
    )

    print(anomaly_detector.summary())

    val_parity = hyperparameters['validation_parity']
    val_event_file_path, val_weight_file_path = parity_file_map[val_parity]
    X_val = load_file(
        tree_name="Events", 
        file_path=val_event_file_path, 
        list_of_branches=branches_to_load,
        convert_to_numpy=True
    )
    
    if cleaning_type == "threshold":
        val_mask, X_val = threshold_cleaning(
            data=X_val,
            pos_thrsh=low,
            neg_thrsh=high
        )
    elif cleaning_type == "quantile":
        val_mask, X_val = quantile_cleaning(
            data=X_val,
            q_low=low,
            q_high=high
        )
    else:
        raise RuntimeError(f"Unsupported data cleaning type {cleaning_type}.")

    y_val = load_file(
        tree_name="weight_tree", 
        file_path=val_weight_file_path, 
        list_of_branches=["class_target"],
        convert_to_numpy=True
    )
    y_val = y_val[val_mask]
    mm.print_memory_usage(msg="After loading validation set")

    bkg_mask_val = y_val == 0
    bkg_mask_val = bkg_mask_val.reshape(-1)
    X_val = X_val[bkg_mask_val]

    # val_mins = np.min(X_val, axis=0)
    # val_maxs = np.max(X_val, axis=0)
    # X_val = (X_val - val_mins)/(val_maxs - val_mins)

    X_val = (X_val - train_mean)/np.sqrt(train_variance)

    callbacks = [
        tf.keras.callbacks.EarlyStopping(patience=20, restore_best_weights=True),
        # tf.keras.callbacks.ReduceLROnPlateau(patience=20)
    ]

    epochs = hyperparameters["epochs"]
    batch_size = hyperparameters["batch_size"]
    train_history = anomaly_detector.fit(
        X_bkg_train,
        X_bkg_train,
        validation_data=(X_val, X_val),
        epochs=epochs,
        batch_size=batch_size,
        callbacks=callbacks,
        verbose=0,
        shuffle=True
    )

    model_dir = os.path.join(cfg["training_directory"], anomaly_detector.name)
    PlotMetric(train_history, anomaly_detector.name, "loss", plotting_dir=model_dir)
    anomaly_detector.save(os.path.join(model_dir, f"{anomaly_detector.name}.keras"))

if __name__ == "__main__":
    main()