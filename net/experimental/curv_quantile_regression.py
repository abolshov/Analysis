import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


from Dataset import Dataset
from NetUtils import TrainableSiLU, EpochLossUpdater
from LayerUtils import PtLayer, EtaLayer, PhiLayer, EnErrLayer

from MiscUtils import *
from PlotUtils import PlotMetric, PlotHist


def main():
    dataset = Dataset('dataset_config.yaml')
    dataset.Load()
    X_train, input_names, y_train, target_names = dataset.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
    os.environ['TF_DETERMINISTIC_OPS'] = '1'
    tf.random.set_seed(42)

    model_name = "model_test"
    model_dir = model_name
    os.makedirs(model_dir, exist_ok=True)

    learning_rate_scheduler = tf.keras.callbacks.LearningRateScheduler(Scheduler)

    input_shape = X_train.shape[1:]
    num_units = 384
    num_layers = 12
    l2_strength = 0.01
    epochs = 50
    batch_size = 256

    # prepare base part of the model
    inputs = tf.keras.layers.Input(shape=input_shape)
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
    # each output head takes all inputs of the base model
    # and does with them whatever it needs to do 
    

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

    # form testing dataset
    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)

    X_test, _, y_test, _ = dataset.Get(TestSelection, 2, 1, 800, 1)
    if standardize:
        X_test = input_scaler.transform(X_test)
    ys_pred = model.predict(X_test)

    assert len(ys_pred) == len(target_names), f"mismatch between number of predicted values and ground truth values: ({len(ys_pred)} vs {len(target_names)})"

    


if __name__ == '__main__':
    main()