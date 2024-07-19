import tensorflow as tf
import pandas as pd
import numpy as np
import uproot

from JetNet_utils import MXLossFunc

class JetNet():
    def __init__(self, features, labels):
        self.features = features
        self.labels = labels

        # training parameters
        self.lr = 0.001
        self.n_epochs = 50
        self.batch_size = 64
        self.verbosity = 0
        self.valid_split = 0.33

        self.model = None


    def ConfigureModel(self, dataset_shape):
        self.model = tf.keras.Sequential([tf.keras.layers.Dense(64, activation='relu'),
                                          tf.keras.layers.Dense(32, activation='relu'),
                                          tf.keras.layers.Dense(16, activation='relu'),
                                          tf.keras.layers.Dense(3)])

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
        if not np.all(test_features.columns == self.features):
            raise RuntimeError(f"Features pased for prediction do not match expected features: passed {test_features.columns}, while expected {self.features}")
        return self.model.predict(test_features)


    def SaveModel(self, path):
        self.model.save(f"{path}JetNet_v1.keras")

    
    def LoadModel(self, model_name):
        # to load with compilation custom_objects needs to be passed
        # couldn't make it work
        # if model is loaded only for prediction, compilation is not needed
        self.model = tf.keras.models.load_model(model_name, compile=False)
        # self.model = tf.keras.models.load_model(model_name, custom_objects={'loss' : MXLossFunc(target, output)})
        print(f"Loaded model {model_name}")
        print(self.model.summary())
