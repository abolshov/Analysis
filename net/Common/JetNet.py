import tensorflow as tf
import pandas as pd
import numpy as np
import uproot
import pandas as pd
import os

from Common.JetNet_utils import * 


class JetNet():
    def __init__(self, cfg):
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
        os.environ['TF_DETERMINISTIC_OPS'] = '1'
        tf.random.set_seed(42)

        # training parameters
        self.lr = cfg['learning_rate']
        self.n_epochs = cfg['n_epochs']
        self.batch_size = cfg['batch_size']
        self.verbosity = cfg['verbosity']
        self.valid_split = cfg['valid_split']
        self.name = cfg['name']
        self.topology = cfg['topology']

        self.model = None


    def ConfigureModel(self, dataset_shape):
        # self.model = tf.keras.Sequential([tf.keras.layers.Dense(layer_size, 
        #                                                         activation='relu',
        #                                                         kernel_initializer='random_normal',
        #                                                         bias_initializer='random_normal') for layer_size in self.topology])

        self.model = tf.keras.Sequential()
        self.model.add(tf.keras.layers.BatchNormalization())
        for layer_size in self.topology:
            AddLayerWithBatchNorm(self.model, layer_size)

        # to predict px, py, pz of H->bb (first three) and px, py, pz of H->WW (last three)
        self.model.add(tf.keras.layers.Dense(6)) 

        self.model.compile(loss=MXLossFunc, optimizer=tf.keras.optimizers.Adam(self.lr))
        self.model.build(dataset_shape)


    def Fit(self, train_features, train_labels):
        if not self.model:
            raise RuntimeError("Model has not been configured before fitting")
        history = self.model.fit(train_features,
                                 train_labels,
                                 validation_split=self.valid_split,
                                 verbose=self.verbosity,
                                 batch_size=self.batch_size,
                                 epochs=self.n_epochs)
        return history


    def Predict(self, test_features):
        output = self.model.predict(test_features)
        pred_df = pd.DataFrame({"X_mass_pred": GetMXPred(output).numpy().ravel()})
        return pred_df


    def SaveModel(self, path):
        if path[-1] == '/':
            self.model.save(f"{path}{self.name}.keras")
        self.model.save(f"{path}/{self.name}.keras")

    
    def LoadModel(self, path_to_model):
        self.model = tf.keras.models.load_model(path_to_model, compile=False)
        print(self.model.summary())
