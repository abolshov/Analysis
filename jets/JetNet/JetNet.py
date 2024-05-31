import tensorflow as tf
import pandas as pd
import numpy as np
import uproot

from JetNet_utils import MXLossFunc

class JetNet():
    def __init__(self):
        self.features = ['genbjet1_px', 'genbjet1_py', 'genbjet1_pz', 'genbjet1_E', 'genbjet2_px', 'genbjet2_py', 'genbjet2_pz', 'genbjet2_E']
        self.labels = ['H_WW_px', 'H_WW_py', 'H_WW_pz', 'H_WW_E', 'X_px', 'X_py', 'X_pz', 'X_E']

        # dataset for training
        self.train_features = pd.DataFrame()
        self.test_features = pd.DataFrame()

        self.train_labels = pd.DataFrame()
        self.test_labels = pd.DataFrame()

        # training parameters
        self.lr = 0.001
        self.n_epochs = 50
        self.batch_size = 64
        self.verbosity = 0
        self.valid_split = 0.33
        self.train_frac = 0.8

        self.model = None
    

    def ReadFile(self, file_name):
        file = uproot.open(file_name)
        tree = file['JetNetTree']
        branches = tree.arrays()

        data_dict = {}
        for name in tree.keys():
            data_dict[name] = np.array(branches[name])
            
        df = pd.DataFrame.from_dict(data_dict)

        train_df = df.sample(frac=self.train_frac, random_state=0)
        test_df = df.drop(train_df.index)

        self.train_features = train_df[self.features].copy()
        self.test_features = test_df[self.features].copy()
        self.train_labels = train_df[self.labels]
        self.test_labels = test_df[self.labels]


    def ConfigureModel(self):
        self.model = tf.keras.Sequential([tf.keras.layers.Dense(32, activation='relu'),
                                          tf.keras.layers.Dense(32, activation='relu'),
                                          tf.keras.layers.Dense(32, activation='relu'),
                                          tf.keras.layers.Dense(3)])

        self.model.compile(loss=MXLossFunc, optimizer=tf.keras.optimizers.Adam(self.lr))
        self.model.build(self.train_features.shape)


    def Fit(self):
        self.ConfigureModel()
        history = self.model.fit(self.train_features,
                                 self.train_labels,
                                 validation_split=self.valid_split,
                                 verbose=self.verbosity,
                                 batch_size=self.batch_size,
                                 epochs=self.n_epochs)

        return history


    def Predict(self):
        return self.model.predict(self.test_features)
