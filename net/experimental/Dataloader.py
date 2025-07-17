import numpy as np
import awkward as ak
import pandas as pd
import uproot
import yaml
import vector
import matplotlib.pyplot as plt
import os
import time


sample_type_map = {1: "GluGluToRadion",
                   2: "GluGluToBulkGraviton", 
                   3 : "VBFToRadion", 
                   4 : "VBFToBulkGraviton"}


def IsKinematic(var):
    if var in ['pt', 'eta', 'phi', 'mass', 'px', 'py', 'pz', 'E']:
        return True
    return False


class Dataloader:
    def __init__(self, branch_loading_cfg):
        self.loader_cfg = None
        self.df = pd.DataFrame()
        self.p4_cache = {}
        with open (branch_loading_cfg, 'r') as branch_loading_cfg_file:
            self.loader_cfg = yaml.safe_load(branch_loading_cfg_file)
        assert self.loader_cfg, "Must contain mapping of branches to load for each object"
        
        self.input_names = self.MakeNameCollection('input_objects')
        self.target_names = self.MakeNameCollection('target_objects')
        
        if self.loader_cfg['plot_cutflow']:
            os.makedirs(self.loader_cfg['plot_dir'], exist_ok=True)

    def MakeNameCollection(self, category):
        res = []
        obj_cfg = self.loader_cfg['objects'][category]
        for obj_name, obj_params in obj_cfg.items():
            dim = obj_params['target_shape']
            features = [f for f in obj_params['branches_to_load'] if not IsKinematic(f)]
            features.extend(obj_params['kinematics_output_format'])
            if dim > 1:
                res.extend([f'{obj_name}_{i + 1}_{f}' for f in features for i in range(dim)])
            else:
                res.extend([f'{obj_name}_{f}' for f in features])
        return res

    def ClearCache(self):
        self.p4_cache.clear()

    def LoadObject(self, obj_name, obj_cfg, tree):
        
        object_present = len([bn for bn in tree.keys() if obj_name in bn]) > 1
        if not object_present:
            # fill dataframe with nans and return
            n_events = tree.num_entries
            nan_dict = {}
            shape = obj_cfg['target_shape']
            for var in obj_cfg['branches_to_load']:
                if IsKinematic(var):
                    continue
                
                if shape > 1:
                    for i in range(shape):
                        nan_dict[f'{obj_name}_{i + 1}_{var}'] = np.full(tree.num_entries, np.nan)
                else:
                    nan_dict[f'{obj_name}_{var}'] = np.full(tree.num_entries, np.nan)
            
            for var in obj_cfg['kinematics_output_format']:
                if shape > 1:
                    for i in range(shape):
                        nan_dict[f'{obj_name}_{i + 1}_{var}'] = np.full(tree.num_entries, np.nan)
                else:
                    nan_dict[f'{obj_name}_{var}'] = np.full(tree.num_entries, np.nan)

            return pd.DataFrame.from_dict(nan_dict)

        # only load in memory branches specified in the config
        arrays = tree.arrays([f'{obj_name}_{var}' for var in obj_cfg['branches_to_load']])
        target_shape = obj_cfg['target_shape']
        kinematics_format = obj_cfg['kinematics_output_format']
        # needs_kinematic_conversion = kinematics_format not None

        # create 2 dicts: one with kinematic variables if conversion is needed
        # another for everything else
        # dict maps variable name (without object prefix) to array of values (possibly > one-dimensional) 

        momentum_branches = {}
        other_branches = {}
        for var in obj_cfg['branches_to_load']:
            branch_name = f'{obj_name}_{var}'
            min_depth, max_depth = arrays[branch_name].layout.minmax_depth
            array = arrays[branch_name]
            
            # add inner dimension to 1D arrays to treat them the same way as branches with nested arrays
            if max_depth == 1:
                array = ak.unflatten(array, 1)

            if IsKinematic(var):
                momentum_branches[var] = ak.fill_none(ak.pad_none(array, target_shape), 0)
            else:
                other_branches[var] = ak.fill_none(ak.pad_none(array, target_shape), 0)                
            
        p4 = vector.zip(momentum_branches)
        p4 = p4[:, :target_shape]
        self.p4_cache[obj_name] = p4 if target_shape > 1 else ak.flatten(p4)

        tmp_dict = {}
        # flattening arrays for each variable if needed
        for i in range(target_shape):
            prefix = f'{obj_name}_{i + 1}' if target_shape > 1 else obj_name
            # first deal with non-kinematic quantities
            for var, arr in other_branches.items():
                tmp_dict[f'{prefix}_{var}'] = arr[:, i].to_numpy()
            
            # now deal with momentum branches
            for var in kinematics_format:
                tmp_dict[f'{prefix}_{var}'] = (lambda x: getattr(x, var))(p4[:, i]).to_numpy()

        return pd.DataFrame.from_dict(tmp_dict)

    def LoadObjects(self, tree):
        # concatenate data from all objects horizontally
        df = pd.DataFrame()
        objects_to_load = self.loader_cfg['objects']
        for obj_type, obj_type_settings in objects_to_load.items():
            for obj_name, cfg in obj_type_settings.items():
                obj_df = self.LoadObject(obj_name, cfg, tree)
                df = pd.concat([df, obj_df], axis=1)
        return df     

    def Load(self, loadable, append=True, apply_selections=None):
        print(f'Loading data...')
        start = time.perf_counter()

        if apply_selections is None:
            apply_selections = self.loader_cfg['apply_selections']

        if isinstance(loadable, list):
            for file_name in loadable:
                self.LoadFile(file_name, append=append, apply_selections=apply_selections)
        elif isinstance(loadable, str):
                self.LoadFile(loadable, append=append, apply_selections=apply_selections)
        else:
            raise TypeError(f'Illegal type {type(loadable)} passed. Only str or list are allowed')

        # shuffle rows
        self.df = self.df.sample(frac=1, random_state=42).reset_index(drop=True)
        end = time.perf_counter()
        elapsed = end - start
        print(f'Done in {elapsed:.2f}s')
        print(f'Shape: {self.df.shape}')

    def LoadFile(self, file_name, append=True, apply_selections=True):
        tree_name = self.loader_cfg['tree_name']
        tree = uproot.open(f'{file_name}:{tree_name}')
        df = pd.concat([self.LoadObjects(tree), self.LoadBranches(tree)], axis=1)

        if apply_selections:
            df = self.ComputeCutflow(df)
            df = self.ApplySelections(df, self.loader_cfg['plot_cutflow'])

        if append:
            self.df = pd.concat([self.df, df], axis=0)
        else:
            self.df = df

        self.ClearCache()
    
    def LoadBranches(self, tree):
        aux_branches_dict = {}
        branch_names = self.loader_cfg['auxilliary_branches']
        for branch in branch_names:
            if branch in tree.keys():
                aux_branches_dict[branch] = tree[branch].array().to_numpy()
            else:
                aux_branches_dict[branch] = np.full(tree.num_entries, np.nan)
        return pd.DataFrame.from_dict(aux_branches_dict)

    def DropUnused(self):
        cols_to_drop = [col for col in self.df.columns if col not in self.input_names + self.target_names]
        self.df.drop(cols_to_drop, axis=1, inplace=True)

    def ComputeCutflow(self, df):
        # compute high-level kinematics needed for selections using p4_cache
        cut_dict = {}
        cut_dict['b1_accept_pt'] = self.p4_cache['genb1'].pt > 20
        cut_dict['b2_accept_pt'] = self.p4_cache['genb2'].pt > 20
        cut_dict['b1_accept_eta'] = np.abs(self.p4_cache['genb1'].eta) < 2.5
        cut_dict['b2_accept_eta'] = np.abs(self.p4_cache['genb2'].eta) < 2.5
        cut_dict['genb1_accept'] = np.logical_and(cut_dict['b1_accept_pt'], cut_dict['b1_accept_eta'])
        cut_dict['genb2_accept'] = np.logical_and(cut_dict['b2_accept_pt'], cut_dict['b2_accept_eta'])
        cut_dict['genb_accept'] = np.logical_and(cut_dict['genb1_accept'], cut_dict['genb2_accept'])
        cut_dict['bb_dr'] = self.p4_cache['genb1'].deltaR(self.p4_cache['genb2'])
        cut_dict['mindR_b1_Jet'] = ak.min(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        cut_dict['mindR_b2_Jet'] = ak.min(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        cut_dict['b1_match_idx'] = ak.argmin(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        cut_dict['b2_match_idx'] = ak.argmin(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        cut_dict['has_two_matches'] = np.logical_and(cut_dict['mindR_b1_Jet'] < 0.4, cut_dict['mindR_b2_Jet'] < 0.4)
        cut_dict['mindR_Hbb_FatJet'] = ak.min(self.p4_cache['genHbb'].deltaR(self.p4_cache['SelectedFatJet']), axis=1)
        cut_dict['is_boosted'] = cut_dict['bb_dr'] < 0.8
        cut_dict['is_resolved'] = ~cut_dict['is_boosted']
        cut_dict['reco_as_boosted'] = np.logical_and(cut_dict['mindR_Hbb_FatJet'] < 0.8, df['nSelectedFatJet'] >= 1)
        cut_dict['reco_2b'] = np.logical_and(cut_dict['has_two_matches'], cut_dict['b1_match_idx'] != cut_dict['b2_match_idx'])
        cut_dict['reco_as_resolved'] = np.logical_and(cut_dict['reco_2b'], df['ncentralJet'] >= 2)
        cut_dict['successful_jet_reco'] = np.logical_or(cut_dict['reco_as_resolved'], cut_dict['reco_as_boosted'])
        cut_dict['reco_lep1'] = np.logical_and(self.p4_cache['lep1'].pt > 0.0, df['lep1_legType'] == df['lep1_gen_kind'])
        cut_dict['reco_lep2'] = np.logical_and(self.p4_cache['lep2'].pt > 0.0, df['lep2_legType'] == df['lep2_gen_kind'])
        cut_dict['successful_lep_reco'] = np.logical_and(cut_dict['reco_lep1'], cut_dict['reco_lep2'])
        cut_dict['successful_reco'] = np.logical_and(cut_dict['successful_lep_reco'], cut_dict['successful_jet_reco'])

        cut_df = pd.DataFrame.from_dict(cut_dict)
        df = pd.concat([df, cut_df], axis=1)
        return df

    def ApplySelections(self, df, plot_cutflow):
        if not plot_cutflow:
            df = df[df['successful_reco']]
            return df
        else:
            cutflow = {}
            cutflow['total'] = len(df)
            df = df[df['b1_accept_pt']]
            cutflow['b1_accept_pt'] = len(df)
            df = df[df['b1_accept_eta']]
            cutflow['b1_accept_eta'] = len(df)
            df = df[df['b2_accept_pt']]
            cutflow['b2_accept_pt'] = len(df)
            df = df[df['b2_accept_eta']]
            cutflow['b2_accept_eta'] = len(df)

            if self.loader_cfg['use_topology'] == 'boosted':
                df = df[df['is_boosted']]
                cutflow['is_boosted'] = len(df)
                df = df[df['reco_lep1']]
                cutflow['reco_lep1'] = len(df)
                df = df[df['reco_lep2']]
                cutflow['reco_lep2'] = len(df)
                df = df[df['reco_as_boosted']]
                cutflow['reco_as_boosted'] = len(df)
            elif self.loader_cfg['use_topology'] == 'resolved':
                df = df[df['is_resolved']]
                cutflow['is_resolved'] = len(df)
                df = df[df['reco_lep1']]
                cutflow['reco_lep1'] = len(df)
                df = df[df['reco_lep2']]
                cutflow['reco_lep2'] = len(df)
                df = df[df['reco_as_resolved']]
                cutflow['reco_as_resolved'] = len(df)
            elif self.loader_cfg['use_topology'] == 'mixed':
                df = df[df['successful_lep_reco']]
                cutflow['successful_lep_reco'] = len(df)
                df = df[df['successful_jet_reco']]
                cutflow['successful_jet_reco'] = len(df)
            else:
                raise RuntimeError(f"Got unknown topology {self.loader_cfg['use_topology']}")

            cut_names = cutflow.keys()
            passed_evt_info = list(cutflow.values())
            if self.loader_cfg['use_percentages']:
                passed_evt_info = [val/passed_evt_info[0]*100 for val in passed_evt_info]
            fig, ax = plt.subplots()
            bar_container = ax.bar(cut_names, passed_evt_info)

            # print(df.describe())
            sample_type = sample_type_map[df['sample_type'].iloc[0]]
            mp = df['X_mass'].iloc[0]
            mass = int(mp) if not np.isnan(mp) else None
            title = f"Cutflow {sample_type} M={mass}" if mp else f"Cutflow {sample_type}"
            ax.set(ylabel='Number of events passed', title=title)
            if self.loader_cfg['use_percentages']:
                ax.bar_label(bar_container, fmt=lambda x: f"{x:.1f}%", fontsize=8)
            else:
                ax.bar_label(bar_container, fmt=lambda x: int(x), fontsize=8)
            plt.xticks(rotation=45, ha='right')

            cutlfow_name = f"cutflow_{sample_type}_M{mass}"
            plt.savefig(os.path.join(self.loader_cfg['plot_dir'], f"{cutlfow_name}.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            return df

    def Get(self, selection=None, *args, **kwargs):
        """
        returns tuple (X, [names of features], y, [names of targets]) or underlying dataframe
        args:
            selection: callable with signature df, *args
        """
        if selection:
            df = self.df[selection(self.df, *args)]
            if 'df' in kwargs and kwargs['df']:
                return df
            return df[self.input_names].values, self.input_names, df[self.target_names].values, self.target_names
        if 'df' in kwargs and kwargs['df']:
            return self.df
        return self.df[self.input_names].values, self.input_names, self.df[self.target_names].values, self.target_names
