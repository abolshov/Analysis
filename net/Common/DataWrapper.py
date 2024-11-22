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
        self.apply_gen_reco_match = cfg['apply_gen_reco_match']
        self.is_bkg = cfg['is_bkg']

        if not self.bb_resolved and not self.bb_boosted:
            raise RuntimeError("Wrong topology configuration provided: check bb_resolved and bb_boosted parameters")

        if not self.qq_resolved and not self.qq_boosted:
            raise RuntimeError("Wrong topology configuration provided: check qq_resolved and qq_boosted parameters")

        # pd.DataFrame containing full dataset
        self.data = None


    def ReadFile(self, file_name):
        file = uproot.open(file_name)
        tree = file[self.tree_name]
        branches = tree.arrays()

        auxilliary_columns = self.labels + self.extra_data
        df = pd.DataFrame({s: branches[s].to_numpy().astype(float) for s in auxilliary_columns})

        # onshell_mass = np.maximum(df['genV1_mass'], df['genV2_mass'])
        # offshell_mass = np.minimum(df['genV1_mass'], df['genV2_mass'])

        # df['genV1_mass'] = onshell_mass
        # df['genV2_mass'] = offshell_mass

        df = AddKinematicFeatures(df, branches, self.n_lep, self.n_jets)
        df = AddJetFeatures(df, branches, self.jet_obs, self.n_jets)

        self.features = [name for name in df.columns if name not in auxilliary_columns]

        if self.apply_fiducial_cut:
            if not self.is_bkg:
                df = ApplyFiducialSelection(df, branches, self.n_lep, self.n_jets, self.apply_gen_reco_match)
            else:
                raise RuntimeError("Invalid attempt to apply fiducial selections for background sample")

        to_drop = [name for name in df.columns if name not in self.features + auxilliary_columns]
        df = df.drop(columns=to_drop)

        if self.data is not None:
            self.data = pd.concat([self.data, df]) 
        else:
            self.data = df

        print(f"Reading from {file_name}: {df.shape[0]} events")


    def ReadFiles(self, input_files):
        print("============START READING FILES============")
        for file_name in input_files:
            self.ReadFile(file_name)
        print("=============END READING FILES=============")


    def Shuffle(self):
        self.data = shuffle(self.data)


    def SelectEvents(self, value, modulo):
        return self.data[self.data['event'] % modulo == value]


    def FormTrainSet(self):
        self.Shuffle()
        self.train_data = self.SelectEvents(self.train_val, self.modulo)

        self.train_labels = self.train_data[self.labels]
        self.train_features = self.train_data[self.features]

        print(f"Training set contains {self.train_features.shape[1]} features for {self.train_features.shape[0]} events")


    def FormTestSet(self):
        self.test_data = self.SelectEvents(self.test_val, self.modulo)
        self.test_features = self.test_data[self.features]
        self.test_labels = self.test_data[self.labels]
        self.test_events = self.test_data['event']


    def Print(self):
        print(self.data.describe())
        print(self.train_data.head())
        for col in self.data.columns:
            print(col)