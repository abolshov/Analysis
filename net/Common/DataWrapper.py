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
        self.lep_obs = cfg['lep_observables']
        self.met_obs = cfg['met_observables']

        if self.n_lep > 2 or self.n_lep == 0:
            raise RuntimeError("Wrong number of leptons provided in configuration: n_lep can be 1 (SL) or 2 (DL)")

        jet_featrues = [f"centralJet{i + 1}_{obs}" for i in range(self.n_jets) for obs in self.jet_obs]
        lep_features = [f"lep{i + 1}_{var}" for i in range(self.n_lep) for var in self.lep_obs]
        met_features = [f"met_{var}" for var in self.met_obs]
        features = jet_featrues + lep_features + met_features
        print(lep_features)
        self.features = features

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

        self.apply_acceptance = cfg['apply_acceptance']

        if not self.bb_resolved and not self.bb_boosted:
            raise RuntimeError("Wrong topology configuration provided: check bb_resolved and bb_boosted parameters")

        if not self.qq_resolved and not self.qq_boosted:
            raise RuntimeError("Wrong topology configuration provided: check qq_resolved and qq_boosted parameters")

        # pd.DataFrame containing full dataset
        self.data = pd.DataFrame(columns=[*self.features, *self.labels, *self.extra_data])


    def ReadFile(self, file_name):
        print(f"Reading from {file_name}:")
        file = uproot.open(file_name)
        tree = file[self.tree_name]
        branches = tree.arrays()

        top_sel = EventTopology(branches, self.bb_resolved, self.bb_boosted, self.qq_resolved, self.qq_boosted)
        selection = top_sel
        if self.apply_acceptance:
            accept_sel = Acceptance(branches, self.n_lep)
            selection = top_sel & accept_sel
        print(f"\tinitial number of events: {len(top_sel)}")

        d1 = {name: np.array(branches[name][selection]) for name in self.extra_data}
        d1["X_mass"] = np.array(branches['X_mass'][selection], dtype=float)

        centralJet_p4 = vector.zip({'pt': branches['centralJet_pt'], 
                                    'eta': branches['centralJet_eta'], 
                                    'phi': branches['centralJet_phi'], 
                                    'mass': branches['centralJet_mass']})

        lep1_p4 = vector.zip({'pt': branches['lep1_pt'], 
                             'eta': branches['lep1_eta'], 
                             'phi': branches['lep1_phi'], 
                             'mass': branches['lep1_mass']})

        zeros = ak.Array(np.zeros(len(branches)))
        lep2_p4 = vector.zip({'pt': zeros, 
                              'eta': zeros, 
                              'phi': zeros, 
                              'mass': zeros})
        if self.n_lep == 2:
            lep2_p4 = vector.zip({'pt': branches['lep2_pt'], 
                                  'eta': branches['lep2_eta'], 
                                  'phi': branches['lep2_phi'], 
                                  'mass': branches['lep2_mass']})

        met_p4 = vector.zip({'pt': branches['PuppiMET_pt'], 
                            'eta': 0, 
                            'phi': branches['PuppiMET_phi'], 
                            'mass': 0})

        centralJet_p4 = centralJet_p4[selection]
        lep1_p4 = lep1_p4[selection]
        lep2_p4 = lep2_p4[selection]
        met_p4 = met_p4[selection]

        PxPyPzE = ['px', 'py', 'pz', 'E']
        func_map = {'px': Px, 'py': Py, 'pz': Pz, 'E': E}

        d2 = {}
        for i in range(self.n_jets):
            for var in self.jet_obs:
                if var in PxPyPzE:
                    func = func_map[var]
                    var_awkward_array = func(centralJet_p4)
                else:
                    branch_name = f"centralJet_{var}"
                    var_awkward_array = branches[branch_name]
                    var_awkward_array = var_awkward_array[selection]
                d2[f"centralJet{i + 1}_{var}"] = GetNumPyArray(var_awkward_array, self.n_jets, i)

        for var in PxPyPzE:
            func = func_map[var]
            var_awkward_array = func(lep1_p4)
            d2[f"lep1_{var}"] = ak.to_numpy(var_awkward_array)
            var_awkward_array = func(lep2_p4)
            d2[f"lep2_{var}"] = ak.to_numpy(var_awkward_array)

        d2["met_px"] = ak.to_numpy(met_p4.px)
        d2["met_py"] = ak.to_numpy(met_p4.py)

        data_dict = d1 | d2

        df = pd.DataFrame.from_dict(data_dict)
        print(f"\tnumber of events selected: {df.shape[0]}")
        self.data = pd.concat([self.data, df]) 


    def ReadFiles(self, input_files):
        print("============START READING FILES============")
        for file_name in input_files:
            self.ReadFile(file_name)
        print("=============END READING FILES=============")


    def Shuffle(self):
        self.data = shuffle(self.data)


    def TrainTestSplit(self):
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

        print(f"Training dataset contains {len(train_df)} events")


    def SelectEvents(self, value, modulo):
        return self.data[self.data['event'] % modulo == value]


    def Print(self):
        for name in self.data.columns:
            if name not in self.features:
                print(name)
        print(self.train_features["lep2_px"].head())
        print(f"Train feature columns ({len(self.train_features.columns)}):")
        for col in self.train_features.columns:
            print(f"\t{col}")
        print(self.train_features.head(5))
        print(self.data.columns)
        print("\n")
        print(f"Data columns ({len(self.data.columns)}):")
        for col in self.data.columns:
            print(f"\t{col}")
