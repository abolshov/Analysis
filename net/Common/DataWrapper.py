import pandas as pd
import numpy as np
import awkward as ak
import uproot
import vector
from sklearn.utils import shuffle
from Common.DataWrapper_utils import *


class DataWrapper():
    def __init__(self, cfg):
        self.n_jets = cfg['n_jets']
        self.n_lep = cfg['n_lep']
        self.jet_obs = cfg['jet_observables']

        if self.n_lep > 2 or self.n_lep == 0:
            raise RuntimeError("Wrong number of leptons provided in configuration: n_lep can be 1 (SL) or 2 (DL)")

        self.tree_name = cfg['tree_name']
        self.labels = cfg['labels']
        self.extra_data = cfg['extra_data']

        self.modulo = cfg['modulo']
        self.train_val = cfg['train_val']
        self.test_val = cfg['test_val']

        self.bb_resolved = cfg['bb_resolved']
        self.qq_resolved = cfg['qq_resolved']

        self.bb_boosted = cfg['bb_boosted']
        self.qq_boosted = cfg['qq_boosted']

        self.apply_fiducial_cut = cfg['apply_fiducial_cut']

        if not self.bb_resolved and not self.bb_boosted:
            raise RuntimeError("Wrong topology configuration provided: check bb_resolved and bb_boosted parameters")

        if not self.qq_resolved and not self.qq_boosted:
            raise RuntimeError("Wrong topology configuration provided: check qq_resolved and qq_boosted parameters")

        # pd.DataFrame containing full dataset
        self.data = None


    def ReadFile(self, file_name):
        print(f"Reading from {file_name}:")
        file = uproot.open(file_name)
        tree = file[self.tree_name]
        branches = tree.arrays()

        auxilliary_columns = self.labels + self.extra_data
        df = pd.DataFrame({s: branches[s].to_numpy() for s in auxilliary_columns})

        AddKinematicFeatures(df, branches, self.n_lep, self.n_jets)
        AddJetFeatures(df, branches, self.jet_obs, self.n_jets)

        self.features = [name for name in df.columns if name not in auxilliary_columns]

        if self.apply_fiducial_cut:
            ApplyFiducialSelection(df, branches, self.n_lep)

        to_drop = [name for name in df.columns if name not in self.features + auxilliary_columns]
        print("columns to drop:")
        print(to_drop)
        df = df.drop(columns=to_drop)

        if self.data:
            self.data = pd.concat([self.data, df]) 
        else:
            self.data = df


    def ReadFiles(self, input_files):
        print("============START READING FILES============")
        for file_name in input_files:
            self.ReadFile(file_name)
        print("=============END READING FILES=============")


    def Shuffle(self):
        self.data = shuffle(self.data)


    def TrainTestSplit(self, is_test):
        if not is_test:
            self.Shuffle()
        train_df = self.SelectEvents(self.train_val, self.modulo)
        test_df = self.SelectEvents(self.test_val, self.modulo)

        self.test_labels = test_df[self.labels] # true X_mass (and nothing else!!!)
        test_df = test_df.drop(self.labels, axis=1)
        
        self.train_events = train_df['event'] 
        self.test_events = test_df['event']

        test_df = test_df.drop(['event'], axis=1)
        train_df = train_df.drop(['event'], axis=1)
       
        self.train_features = train_df[self.features] # contain all possible centralJet variables for all central jets (for train)
        self.test_features = test_df[self.features] # contain all possible centralJet variables for all central jets (for test)
        self.train_labels = train_df[self.labels] # contain X_mass

        if is_test:
            print(f"Test dataset contains {len(test_df)} events")
        else:
            print(f"Training dataset contains {len(train_df)} events")


    def SelectEvents(self, value, modulo):
        return self.data[self.data['event'] % modulo == value]


    def FormTrainSet(self):
        self.Shuffle()
        self.train_data = SelectEvents(self.train_val, self.modulo)

        self.train_labels = self.train_data[self.labels]
        self.train_features = self.train_data[self.features]


    def FormTestSet(self):
        self.test_data = SelectEvents(self.test_val, self.modulo)
        self.test_features = self.test_data[self.features]


    def Print(self):
        print(self.data.describe())
        for col in self.data.columns:
            print(col)