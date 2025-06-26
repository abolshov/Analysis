import pandas as pd
import numpy as np
import awkward as ak
import uproot
import vector
import itertools
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from DataWrapper_utils import *


class DataWrapper():
    def __init__(self):
        self.n_jets = 10
        self.n_lep = 2
        self.n_fatjets = 3
        self.jet_obs = ['btagPNetB',
                        'PNetRegPtRawCorr',
                        'PNetRegPtRawCorrNeutrino',
                        'PNetRegPtRawRes',
                        'btagPNetCvB',
                        'btagPNetCvL',
                        'btagPNetCvNotB',
                        'btagPNetQvG']

        self.fatjet_obs = ['particleNet_QCD',
                           'particleNetWithMass_QCD',
                           'particleNetWithMass_HbbvsQCD',
                           'particleNet_XbbVsQCD',
                           'particleNet_massCorr']

        self.feature_list = [f"jet{i + 1}_{var}" for i in range(self.n_jets) for var in self.jet_obs + ['px', 'py', 'pz', 'E']]
        self.feature_list.extend([f"fatjet{i + 1}_{var}" for i in range(self.n_fatjets) for var in self.fatjet_obs + ['px', 'py', 'pz', 'E']])
        self.feature_list.extend([f"lep{i + 1}_{var}" for i in range(2) for var in ['px', 'py', 'pz', 'E']])
        self.feature_list.extend(["met_px", "met_py", "met_E"])

        self.tree_name = "Events"
        self.extra_data = ["event"]

        self.labels = ["genHbb_px", "genHbb_py", "genHbb_pz", "genHbb_E",
                       "genHVV_px", "genHVV_py", "genHVV_pz", "genHVV_E"] 

        # self.labels = ["offshellV_px", 
        #                "offshellV_py", 
        #                "offshellV_pz", 
        #                "offshellV_E"] 

        # self.labels = ["offshellV_px", 
        #                "offshellV_py", 
        #                "offshellV_pz", 
        #                "offshellV_mass"] 

        self.modulo = 2
        self.train_val = 0
        self.test_val = 1

        self.apply_fiducial_cut = True

        self.data = None
        self.feature_scaler = None
        self.label_scaler = None


    def ReadFile(self, file_name):
        file = uproot.open(file_name)
        tree = file[self.tree_name]
        branches = tree.arrays()

        auxilliary_columns = self.extra_data
        df = pd.DataFrame({s: branches[s].to_numpy().astype(float) for s in auxilliary_columns})

        df = AddKinematicFeatures(df, branches, self.n_lep, self.n_jets)
        df = AddJetFeatures(df, branches, self.jet_obs, self.n_jets)
        df = AddFatJetFeatures(df, branches, self.fatjet_obs, self.n_fatjets)
        df = AddKinematicLabels(df, branches)

        self.features = [name for name in df.columns if name not in auxilliary_columns]

        if self.apply_fiducial_cut:
            df = ApplyFiducialSelection(df, branches, self.n_lep, self.n_jets, self.n_fatjets)

        to_keep = self.features + auxilliary_columns + self.labels
        to_drop = [name for name in df.columns if name not in to_keep]
        df = df.drop(columns=to_drop)

        if self.n_lep == 2:
            df = df.drop(columns=["lep1_pt", "lep2_pt"])
            self.features.remove("lep2_pt")
        else:
            df = df.drop(columns=["lep1_pt"])

        self.features.remove("lep1_pt")

        if self.data is not None:
            self.data = pd.concat([self.data, df]) 
        else:
            self.data = df

        print(f"Reading from {file_name}: {df.shape[0]} events")


    def ReadFiles(self, input_files):
        print("============START READING FILES============")
        for file_name in input_files:
            self.ReadFile(file_name)
            # break
        print("=============END READING FILES=============")


    def Shuffle(self):
        self.data = shuffle(self.data)


    def SelectEvents(self, value, modulo):
        return self.data[self.data['event'] % modulo == value]


    def FormTrainSet(self):
        self.Shuffle()
        self.train_data = self.SelectEvents(self.train_val, self.modulo)
        self.train_data = self.train_data.drop(columns=self.extra_data)

        self.train_labels = self.train_data[[col for col in self.train_data.columns if col not in self.feature_list]]
        self.train_features = self.train_data[[col for col in self.train_data.columns if col in self.feature_list]]

        self.train_features_scaled = self.StandardizeFeatures(self.train_features)
        self.train_labels_scaled = self.StandardizeLabels(self.train_labels)

        # print(self.train_features.shape)
        # print(self.train_labels.shape)

        print(f"Training set contains {self.train_features.shape[1]} features for {self.train_features.shape[0]} events")


    def FormTestSet(self):
        self.test_data = self.SelectEvents(self.test_val, self.modulo)
        self.test_events = self.test_data["event"]
        self.test_data = self.test_data.drop(columns=self.extra_data)

        self.test_labels = self.test_data[[col for col in self.test_data.columns if col not in self.feature_list]]
        self.test_features = self.test_data[[col for col in self.test_data.columns if col in self.feature_list]]

        self.test_features_scaled = self.StandardizeFeatures(self.test_features)


    def AugmentDataset(self):
        original = self.data.copy()
        directions = ['px', 'py', 'pz']

        # reflection against individual directions
        for direction in directions:
            copy = original.copy()
            for col in copy.columns:
                if direction in col:
                    copy[col] = -1.0*copy[col]
            self.data = pd.concat([self.data, copy])

        # reflection of pairs of directions
        for pair in list(itertools.combinations(directions, 2)):
            copy = original.copy()
            d1, d2 = pair
            for col in copy.columns:
                if d1 in col or d2 in col:
                    copy[col] = -1.0*copy[col]
            self.data = pd.concat([self.data, copy])

        # reflect all 3 directions simultaneously        
        copy = original.copy()   
        for col in copy.columns:
            if 'px' in col or 'py' in col or 'pz' in col:
                copy[col] = -1.0*copy[col]
        self.data = pd.concat([self.data, copy])


    def StandardizeFeatures(self, df):
        if self.feature_scaler is None:
            self.feature_scaler = StandardScaler()
            self.feature_scaler.fit(df)
        scaled_values = self.feature_scaler.fit_transform(df)
        return pd.DataFrame(scaled_values, columns=df.columns)


    def StandardizeLabels(self, df):
        if self.label_scaler is None:
            self.label_scaler = StandardScaler()
            self.label_scaler.fit(df)
        scaled_values = self.label_scaler.fit_transform(df)
        return pd.DataFrame(scaled_values, columns=df.columns)


    def Clear(self):
        self.data = None