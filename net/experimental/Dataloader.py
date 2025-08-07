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
                   3: "VBFToRadion", 
                   4: "VBFToBulkGraviton",
                   5: "DY",
                   8: "TTbar"}


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
        obj_max_depth = None
        for var in obj_cfg['branches_to_load']:
            branch_name = f'{obj_name}_{var}'
            min_depth, max_depth = arrays[branch_name].layout.minmax_depth
            array = arrays[branch_name]
            obj_max_depth = max_depth
            
            # add inner dimension to 1D arrays to treat them the same way as branches with nested arrays
            if max_depth == 1:
                array = ak.unflatten(array, 1)

            if IsKinematic(var):
                momentum_branches[var] = ak.fill_none(ak.pad_none(array, target_shape), 0)
            else:
                other_branches[var] = ak.fill_none(ak.pad_none(array, target_shape), 0)                
            
        p4 = vector.zip(momentum_branches)
        p4 = p4[:, :target_shape]
        self.p4_cache[obj_name] = p4 if target_shape > 1 or obj_max_depth > 1 else ak.flatten(p4)

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

    def Load(self, loadable, append=True, selection_cfg=None):
        print(f'Loading data...')
        start = time.perf_counter()

        if selection_cfg is None:
            selection_cfg = self.loader_cfg['selection_config']

        if isinstance(loadable, list):
            for file_name in loadable:
                self.LoadFile(file_name, append=append, selection_cfg=selection_cfg)
        elif isinstance(loadable, str):
                self.LoadFile(loadable, append=append, selection_cfg=selection_cfg)
        else:
            raise TypeError(f'Illegal type {type(loadable)} passed. Only str or list are allowed')

        # shuffle rows
        self.df = self.df.sample(frac=1, random_state=42).reset_index(drop=True)
        end = time.perf_counter()
        elapsed = end - start
        print(f'Done in {elapsed:.2f}s')
        print(f'Shape: {self.df.shape}')

    def LoadFile(self, file_name, append=True, selection_cfg=None):
        tree_name = self.loader_cfg['tree_name']
        tree = uproot.open(f'{file_name}:{tree_name}')
        df = pd.concat([self.LoadObjects(tree), self.LoadBranches(tree)], axis=1)

        if selection_cfg:
            is_bkg = np.all(df['sample_type'] > 4)
            cuts_df = self.ComputeCutflow(df, is_bkg, selection_cfg)
            df = self.ApplySelections(df, cuts_df, self.loader_cfg['plot_cutflow'])

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

    def ComputeGenLevelCuts(self, df, cuts):
        # compute signal specific cuts on gen level objects in order they will be applied
        b1_accept_pt = self.p4_cache['genb1'].pt > 20
        b1_accept_eta = np.abs(self.p4_cache['genb1'].eta) < 2.5
        cuts['b1_accept'] = np.logical_and(b1_accept_pt, b1_accept_eta)
        
        b2_accept_pt = self.p4_cache['genb2'].pt > 20
        b2_accept_eta = np.abs(self.p4_cache['genb2'].eta) < 2.5
        cuts['b2_accept'] = np.logical_and(b2_accept_pt, b2_accept_eta)

        cuts['lep_accept'] = np.logical_and(df['genV1prod1_pt'] > 5.0, df['genV2prod1_pt'] > 5.0)
        return cuts

    def ComputeRecoLevelCuts(self, df, cuts):
        # compute signal specific cuts on gen level objects in order they will be applied
        ak4_selection = np.logical_and(self.p4_cache['centralJet'].pt > 20.0, np.abs(self.p4_cache['centralJet'].eta) < 2.5)
        has_ak4_jets = ak.count_nonzero(ak4_selection, axis=1) >= 2
        ak8_selection = np.logical_and(self.p4_cache['SelectedFatJet'].pt > 200.0, np.abs(self.p4_cache['SelectedFatJet'].eta) < 2.5)
        has_ak8_jets = ak.count_nonzero(ak8_selection, axis=1) >= 1
        cuts['has_reco_jets'] = np.logical_or(has_ak8_jets, has_ak4_jets)
        cuts['has_reco_leptons'] = np.logical_and(self.p4_cache['lep1'].pt > 0.0, self.p4_cache['lep2'].pt > 0.0)
        return cuts

    def ComputeGenRecoMatchingCuts(self, df, cuts):
        topology = self.loader_cfg['use_topology']
        bb_dr = self.p4_cache['genb1'].deltaR(self.p4_cache['genb2'])
        is_boosted = bb_dr < 0.8

        mindR_b1_Jet = ak.min(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        mindR_b2_Jet = ak.min(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        b1_match_idx = ak.argmin(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        b2_match_idx = ak.argmin(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        has_two_matches = np.logical_and(mindR_b1_Jet < 0.4, mindR_b2_Jet < 0.4)
        reco_2b = np.logical_and(has_two_matches, b1_match_idx != b2_match_idx)

        mindR_Hbb_FatJet = ak.min(self.p4_cache['genHbb'].deltaR(self.p4_cache['SelectedFatJet']), axis=1)

        match topology:
            case 'resolved':
                cuts['is_resolved'] = ~is_boosted
                cuts['reco_as_resolved'] = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
            case 'boosted':
                cuts['is_boosted'] = is_boosted
                cuts['reco_as_boosted'] = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
            case 'mixed':
                as_fat = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
                as_slim = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
                cuts['successful_jet_reco'] = np.logical_or(as_fat, as_slim)
            case _:
                raise RuntimeError(f'Illegal use_topology value: {topology}')

        # lepton matching
        reco_lep1 = df['lep1_legType'] == df['lep1_gen_kind']
        reco_lep2 = df['lep2_legType'] == df['lep2_gen_kind']
        cuts['successful_lep_reco'] = np.logical_and(reco_lep1, reco_lep2)
        return cuts

    def ComputeCutflow(self, df, is_bkg, selection_cfg):
        # I want it to return list of cuts it computed in certain order for plotting cutflow
        assert selection_cfg is not None, 'Selection config must be provided to compute cutflow'
        
        cut_dict = {}
        if not is_bkg and selection_cfg['gen_level']:
            cut_dict = self.ComputeGenLevelCuts(df, cut_dict)
        
        if selection_cfg['reco_level']:
            cut_dict.update(self.ComputeRecoLevelCuts(df, cut_dict))

        if not is_bkg and selection_cfg['gen_reco_matching']:
            cut_dict.update(self.ComputeGenRecoMatchingCuts(df, cut_dict))
    
        cut_df = pd.DataFrame.from_dict(cut_dict)
        return cut_df

    def ApplySelections(self, df, cut_df, plot_cutflow):
        cutflow = {}
        valid = np.full(len(cut_df), True)
        cutflow['total'] = len(cut_df)
        for cut_name in cut_df.columns:
            cut = cut_df[cut_name]
            valid = np.logical_and(cut, valid)
            cutflow[cut_name] = np.sum(valid)

        df = df[valid]

        if plot_cutflow:        
            cut_names = cutflow.keys()
            passed_evt_info = list(cutflow.values())
            if self.loader_cfg['use_percentages']:
                passed_evt_info = [val/passed_evt_info[0]*100 for val in passed_evt_info]
            fig, ax = plt.subplots()
            bar_container = ax.bar(cut_names, passed_evt_info)

            sample_type = sample_type_map[df['sample_type'].iloc[0]]
            mp = df['X_mass'].iloc[0]
            mass = int(mp) if not np.isnan(mp) else None
            title = f"Cutflow {sample_type} M={mass}" if mp and mp > 0 else f"Cutflow {sample_type}"
            ax.set(ylabel='Number of events passed', title=title)
            if self.loader_cfg['use_percentages']:
                ax.bar_label(bar_container, fmt=lambda x: f"{x:.1f}%", fontsize=8)
            else:
                ax.bar_label(bar_container, fmt=lambda x: int(x), fontsize=8)
            plt.xticks(rotation=45, ha='right')

            cutlfow_name = f"cutflow_{sample_type}"
            if mp > 0:
                cutlfow_name += f'_{mp}'
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
