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

from typing import List
from PlotUtils import PlotMetric, PlotConfMatrix, PlotROC, PlotPRC
from MiscUtils import MemoryMonitor, load_file, nearest_pow2, map_input_files, clean_extreme_values
from sklearn.utils import class_weight
from LayerUtils import DeepResidualBlock
import gc

tf.random.set_seed(42)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'

def main():
    mm = MemoryMonitor()
    mm.print_memory_usage(msg="Training DNN for event classification")

    weight_file_pattern = re.compile(r"nParity(\d)_Merged_weight.root")
    event_file_pattern = re.compile(r"nParity(\d)_Merged.root")
    train_file_dir = "/home/artem/Desktop/CMS/data/DNN/SL/resolved/Run3_2022/Dataset"
    parity_file_map = map_input_files(directory=train_file_dir,
                                      weight_file_pattern=weight_file_pattern,
                                      event_file_pattern=event_file_pattern)
    for p, (event_file, weight_file) in parity_file_map.items():
        print(f"Parity {p}: event_file={event_file}, weight_file={weight_file}")


    with open ('classifier_cfg.yaml', 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    hyperparameters = cfg['hyperparameters']

    train_parity = hyperparameters['training_parity']
    train_event_file_path, train_weight_file_path = parity_file_map[train_parity]

    use_extra_variables = hyperparameters['use_extra_branches']

    branches_to_load = cfg['base_branches']
    if use_extra_variables:
        branches_to_load.extend(cfg['extra_branches'])

    X_train = load_file(tree_name="Events", 
                        file_path=train_event_file_path, 
                        list_of_branches=branches_to_load,
                        convert_to_numpy=True)
    
    mm.print_memory_usage(msg=f"Before cleaning: {X_train.shape}")
    clean_mask_train, X_train = clean_extreme_values(data=X_train)
    mm.print_memory_usage(msg=f"After cleaning: {X_train.shape}")

    weights_train = load_file(tree_name="weight_tree", 
                              file_path=train_weight_file_path, 
                              list_of_branches=["class_weight", "class_target"],
                              convert_to_numpy=False)

    train_sample_weights = weights_train["class_weight"].to_numpy() # may add later for sample_weights
    y_train = weights_train["class_target"].to_numpy()
    y_train = y_train[clean_mask_train]

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

    val_parity = hyperparameters['validation_parity']
    val_event_file_path, val_weight_file_path = parity_file_map[val_parity]
    X_val = load_file(tree_name="Events", 
                      file_path=val_event_file_path, 
                      list_of_branches=branches_to_load,
                      convert_to_numpy=True)
    
    weights_val = load_file(tree_name="weight_tree", 
                            file_path=val_weight_file_path, 
                            list_of_branches=["class_weight", "class_target"],
                            convert_to_numpy=False)

    val_sample_weights = weights_val["class_weight"].to_numpy() # may add later for sample_weights
    y_val = weights_val["class_target"].to_numpy()
    if not multiclass:
        y_val = 1 - y_val

    mm.print_memory_usage(msg=f"Loaded {X_val.shape[1]} variables for {X_val.shape[0]} events for validation set")

    # precompute mean and variance
    train_mean = np.mean(X_train, axis=0)
    train_variance = np.var(X_train, axis=0)

    assert train_mean.shape == train_variance.shape and train_mean.shape[0] == X_train.shape[1]

    # define model
    num_hidden_layers = hyperparameters['num_hidden_layers']
    num_units = 2 * nearest_pow2(X_train.shape[1])
    dropout = hyperparameters['dropout']
    print(f"Using {num_hidden_layers} hidden layers with {num_units} units in each layer")

    inputs = tf.keras.Input(shape=(X_train.shape[1],))
    x = tf.keras.layers.Normalization(mean=train_mean, variance=train_variance)(inputs)
    x = tf.keras.layers.Dense(num_units, use_bias=False)(x)
    for _ in range(num_hidden_layers):
        # x = tf.keras.layers.Dense(num_units)(x)
        # x = tf.keras.layers.BatchNormalization()(x)
        # x = ResidualBlock(units=num_units, activation="swish")(x)
        x = DeepResidualBlock(units=num_units, activation="swish")(x)
        # x = tf.keras.activations.swish(x)
        x = tf.keras.layers.Dropout(dropout)(x)

    # output layer
    output_bias = cfg['output_bias']
    if output_bias is None:
        bias_initializer = "zeros"
    elif isinstance(output_bias, str) and output_bias == 'exact':
        num_pos_train = np.sum(y_train)
        tot_train = len(y_train)
        num_neg_train = tot_train - num_pos_train
        proba = num_pos_train/num_neg_train
        initial_bias = np.log(proba)
        bias_initializer = tf.keras.initializers.Constant(initial_bias)
    elif isinstance(output_bias, float):
        bias_initializer = tf.keras.initializers.Constant(output_bias)
    else:
        raise RuntimeError(f"Illegal option `{output_bias}` for `output_bias`. Allowed `'exact'`, `None` or specifica floating point value.")

    x = tf.keras.layers.Dense(1, bias_initializer=bias_initializer)(x)
    
    output_activation = tf.keras.activations.get(cfg['output_activation'])
    outputs = output_activation(x)

    model_name = cfg['model_name']
    model = tf.keras.models.Model(inputs=inputs, 
                                  outputs=outputs, 
                                  name=model_name)
    
    model_dir = os.path.join(cfg['training_directory'], model.name)
    os.makedirs(model_dir, exist_ok=True)
    print(model.summary())
    tf.keras.utils.plot_model(model, os.path.join(model_dir, f"summary_{model.name}.pdf"), show_shapes=True)

    metrics = [
        tf.keras.metrics.TruePositives(name='TruePositives'),
        tf.keras.metrics.FalsePositives(name='FalsePositives'),
        tf.keras.metrics.TrueNegatives(name='TrueNegatives'),
        tf.keras.metrics.FalseNegatives(name='FalseNegatives'),
        tf.keras.metrics.BinaryAccuracy(name='Accuracy'),
        tf.keras.metrics.Precision(name='Precision'),
        tf.keras.metrics.Recall(name='Recall'),
        tf.keras.metrics.AUC(name='AUC'),
        tf.keras.metrics.AUC(name='PRC', curve='PR'),
        # tf.keras.metrics.F1Score(name='F1'),
    ]

    # learning_rate = hyperparameters['learning_rate']
    loss_cfg = cfg["loss"]
    loss_instance = tf.keras.losses.get(loss_cfg)
    optimizer_cfg = cfg['optimizer']
    optimizer_instance = tf.keras.optimizers.get(optimizer_cfg)
    model.compile(
        loss=loss_instance,
        optimizer=optimizer_instance,
        metrics=metrics
    )

    # fit model
    batch_size = hyperparameters['batch_size']
    epochs = hyperparameters['epochs']
    print("Training the model ...")
    history = model.fit(X_train, 
                        y_train, 
                        shuffle=True,
                        validation_data=(X_val, y_val),
                        class_weight=class_weight_dict,
                        verbose=0,
                        batch_size=batch_size,
                        epochs=epochs)
    print("... Done!")
    
    PlotMetric(history, model.name, "loss", plotting_dir=model_dir)
    for m in metrics:
        PlotMetric(history, model.name, m.name, plotting_dir=model_dir)

    model.save(os.path.join(model_dir, f"{model.name}.keras"))

    del X_train
    del y_train
    del X_val
    del y_val
    coll = gc.collect()
    mm.print_memory_usage(msg=f"After invoking garbage collector, collected {coll}.")

    print("Loading test set")
    test_parity = hyperparameters['test_parity']
    test_event_file_path, test_weight_file_path = parity_file_map[test_parity]
    X_test = load_file(tree_name="Events", 
                       file_path=test_event_file_path, 
                       list_of_branches=branches_to_load,
                       convert_to_numpy=True)
    
    weights_test = load_file(tree_name="weight_tree", 
                             file_path=test_weight_file_path, 
                             list_of_branches=["class_weight", "class_target"],
                             convert_to_numpy=False)

    y_test = weights_test["class_target"].to_numpy()
    if not multiclass:
        y_test = 1 - y_test
    
    mm.print_memory_usage(msg=f"Loaded {X_test.shape[1]} variables for {X_test.shape[0]} events for test set")

    test_pred = model.predict(X_test, batch_size=batch_size)
    
    threshold = cfg['classification_threshold']
    PlotConfMatrix(labels=y_test,
                   predictions=test_pred,
                   threshold=threshold,
                   plotdir=model_dir)

    PlotROC(labels=y_test,
            predictions=test_pred,
            plotdir=model_dir)
    
    PlotPRC(labels=y_test,
            predictions=test_pred,
            plotdir=model_dir)

    print("Evaluate model on test set ...")
    test_history = model.evaluate(X_test,
                                  y_test, 
                                  batch_size=batch_size,
                                  verbose=0)
    print("Test set results:")
    loss = test_history[0]
    print(f"\tLoss: {loss:.4f}")
    for i, m in enumerate(metrics):
        print(f"\t{m.name}: {test_history[i + 1]:.4f}")

    # free resources (just in case)
    tf.keras.backend.clear_session()

if __name__ == "__main__":
    main()