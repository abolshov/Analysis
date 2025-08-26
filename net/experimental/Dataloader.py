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
        self.object_feature_names = {}
        self.obj_shapes = {}
        # self.objects = {} # maps object name -> df with features of htat object
        with open (branch_loading_cfg, 'r') as branch_loading_cfg_file:
            self.loader_cfg = yaml.safe_load(branch_loading_cfg_file)
        assert self.loader_cfg, "Must contain mapping of branches to load for each object"
        
        self.input_names = self.MakeNameCollection('input_objects')
        self.target_names = self.MakeNameCollection('target_objects')
        
        if self.loader_cfg['plot_cutflow']:
            os.makedirs(self.loader_cfg['plot_dir'], exist_ok=True)

        # should it be provided in cfg or calculated based on how many lep objects are specified?
        self.channel = self.loader_cfg['channel']
        self.bb_topology = self.loader_cfg['bb_topology']
        self.qq_topology = self.loader_cfg['qq_topology']

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

        # fill dict with object name -> object feature names
        self.obj_shapes[obj_name] = obj_cfg['target_shape']

        if obj_cfg['target_shape'] > 1:
            self.object_feature_names[obj_name] = [f'{obj_name}_{idx + 1}_{branch}' for branch in obj_cfg['branches_to_load'] if not IsKinematic(branch) for idx in range(obj_cfg['target_shape'])]
            self.object_feature_names[obj_name].extend([f'{obj_name}_{idx + 1}_{branch}' for branch in obj_cfg['kinematics_output_format'] for idx in range(obj_cfg['target_shape'])])
        else:
            self.object_feature_names[obj_name] = [f'{obj_name}_{branch}' for branch in obj_cfg['branches_to_load'] if not IsKinematic(branch)]
            self.object_feature_names[obj_name].extend([f'{obj_name}_{branch}' for branch in obj_cfg['kinematics_output_format']])

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

        match self.channel:
            case 'DL':
                cuts['lep_accept'] = np.logical_and(df['genV1prod1_pt'] > 5.0, df['genV2prod1_pt'] > 5.0)
            case 'SL':
                q1_accept_pt = self.p4_cache['genV2prod1'].pt > 20
                q1_accept_eta = np.abs(self.p4_cache['genV2prod1'].eta) < 5.0
                cuts['q1_accept'] = np.logical_and(q1_accept_pt, q1_accept_eta)
                
                q2_accept_pt = self.p4_cache['genV2prod2'].pt > 20
                q2_accept_eta = np.abs(self.p4_cache['genV2prod2'].eta) < 5.0
                cuts['q2_accept'] = np.logical_and(q2_accept_pt, q2_accept_eta)

                cuts['lep_accept'] = df['genV1prod1_pt'] > 5.0
            case _:
                raise RuntimeError(f'Illegal channel {self.channel}, only SL or DL are allowed')
        return cuts

    def ComputeRecoLevelCuts(self, df, cuts):
        # compute signal specific cuts on gen level objects in order they will be applied
        has_reco_jets = None
        match self.channel:
            case 'DL':
                ak4_selection = np.logical_and(self.p4_cache['centralJet'].pt > 20.0, np.abs(self.p4_cache['centralJet'].eta) < 2.5)
                has_ak4_jets = ak.count_nonzero(ak4_selection, axis=1) >= 2
                ak8_selection = np.logical_and(self.p4_cache['SelectedFatJet'].pt > 200.0, np.abs(self.p4_cache['SelectedFatJet'].eta) < 2.5)
                has_ak8_jets = ak.count_nonzero(ak8_selection, axis=1) >= 1
                has_reco_jets = np.logical_or(has_ak8_jets, has_ak4_jets)
                cuts['has_reco_leptons'] = np.logical_and(self.p4_cache['lep1'].pt > 0.0, self.p4_cache['lep2'].pt > 0.0)
            case 'SL':
                ak4_jet_pt_spec = self.p4_cache['centralJet'].pt > 20.0
                ak4_bjet_eta_spec = np.abs(self.p4_cache['centralJet'].eta) < 2.5
                ak4_light_jet_eta_spec = np.abs(self.p4_cache['centralJet'].eta) < 5.0
        
                ak8_jet_pt_spec = self.p4_cache['SelectedFatJet'].pt > 200.0
                ak8_bjet_eta_spec = np.abs(self.p4_cache['SelectedFatJet'].eta) < 2.5
                ak8_light_jet_eta_spec = np.abs(self.p4_cache['SelectedFatJet'].eta) < 5.0

                valid_ak8q = np.logical_and(ak8_jet_pt_spec, ak8_light_jet_eta_spec)
                valid_ak8b = np.logical_and(ak8_jet_pt_spec, ak8_bjet_eta_spec)
                valid_ak4q = np.logical_and(ak4_jet_pt_spec, ak4_light_jet_eta_spec)
                valid_ak4b = np.logical_and(ak4_jet_pt_spec, ak4_b_jet_eta_spec)
                num_valid_ak4q = ak.count_nonzero(valid_ak4q, axis=1)
                num_valid_ak4b = ak.count_nonzero(valid_ak4b, axis=1)
                num_valid_ak8q = ak.count_nonzero(valid_ak8q, axis=1)
                num_valid_ak8b = ak.count_nonzero(valid_ak8b, axis=1)

                # >= 2 ak4 b jet cand with pt > 20 and eta < 2.5 and >= 2 ak4 light cand with pt > 20 and eta < 5  
                # valid_ak4b and valid_ak4q have the SAME dimensions => CAN or them
                has_2ak4b_2ak4q = ak.count_nonzero(np.logical_or(valid_ak4b, valid_ak4q), axis=1) >= 4

                # >= 2 ak4 b jet cand with pt > 20 and eta < 2.5 and >= 1 ak8 light cand with pt > 200 and eta < 5
                # valid_ak4b and valid_ak8q have DIFFERENT dimensions => CAN NOT or them
                has_2ak4b_1ak8q = np.logical_and(num_valid_ak4b >= 2, num_valid_ak8q >= 1)

                # >= 1 ak8 b jet cand with pt > 200 and eta < 2.5 and >= 2 ak4 light cand with pt > 20 and eta < 5 
                has_1ak8b_2ak4q = np.logical_and(num_valid_ak4q >=2, num_valid_ak8b >= 1)

                # >= 1 ak8 b jet cand with pt > 200 and eta < 2.5 and >= 1 ak8 light cand with pt > 200 and eta < 5
                has_1ak8b_1ak8q = ak.count_nonzero(np.logical_or(valid_ak8b, valid_ak8q), axis=1) >= 2

                # combine 4 cases: either (2ak4b, 2ak4q) or (2ak4b, 1ak8q) or (1ak8b, 2ak4q) or (1ak8b, 1ak8q)
                has_reco_jets = np.any(np.stack([has_2ak4b_2ak4q, has_2ak4b_1ak8q, has_1ak8b_2ak4q, has_1ak8b_1ak8q], axis=1), axis=1)

                cuts['has_reco_leptons'] = self.p4_cache['lep1'].pt > 0.0
            case _:
                raise RuntimeError(f'Illegal channel {self.channel}, only SL or DL are allowed')
        
        assert has_reco_jets is not None
        cuts['has_reco_jets'] = has_reco_jets
        return cuts

    def _MatchBQuarks(self, df):
        bquark_cuts = {}

        bb_dr = self.p4_cache['genb1'].deltaR(self.p4_cache['genb2'])
        boosted_bb = bb_dr < 0.8

        mindR_b1_Jet = ak.min(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        mindR_b2_Jet = ak.min(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        self._b1_match_idx = ak.argmin(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1) 
        self._b2_match_idx = ak.argmin(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        has_two_matches = np.logical_and(mindR_b1_Jet < 0.4, mindR_b2_Jet < 0.4)
        reco_2b = np.logical_and(has_two_matches, self._b1_match_idx != self._b2_match_idx)
        mindR_Hbb_FatJet = ak.min(self.p4_cache['genHbb'].deltaR(self.p4_cache['SelectedFatJet']), axis=1)
        self._fatbb_idx = ak.argmin(self.p4_cache['genHbb'].deltaR(self.p4_cache['SelectedFatJet']), axis=1) 

        match self.bb_topology:
            case 'resolved':
                bquark_cuts['resolved_bb'] = ~boosted_bb
                bquark_cuts['reco_resolved_bb'] = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
            case 'boosted':
                bquark_cuts['boosted_bb'] = boosted_bb
                bquark_cuts['reco_boosted_bb'] = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
            case 'mixed':
                as_fat = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
                as_slim = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
                bquark_cuts['reco_bb'] = np.logical_or(as_fat, as_slim)
            case _:
                raise RuntimeError(f'Illegal use_topology value: {self.bb_topology}')
        return bquark_cuts


    def _MatchLightQuarks(self, df):
        light_quark_cuts = {}

        qq_dr = self.p4_cache['genV2prod1'].deltaR(self.p4_cache['genV2prod2'])
        boosted_qq = qq_dr < 0.8

        mindR_q1_Jet = ak.min(self.p4_cache['genV2prod1'].deltaR(self.p4_cache['centralJet']), axis=1)
        mindR_q2_Jet = ak.min(self.p4_cache['genV2prod2'].deltaR(self.p4_cache['centralJet']), axis=1)
        q1_match_idx = ak.argmin(self.p4_cache['genV2prod1'].deltaR(self.p4_cache['centralJet']), axis=1) 
        q2_match_idx = ak.argmin(self.p4_cache['genV2prod2'].deltaR(self.p4_cache['centralJet']), axis=1)
        has_two_matches = np.logical_and(mindR_q1_Jet < 0.4, mindR_q2_Jet < 0.4)
        reco_2q = np.logical_and(has_two_matches, q1_match_idx != q2_match_idx)
        mindR_Hvv_FatJet = ak.min(self.p4_cache['genV2'].deltaR(self.p4_cache['SelectedFatJet']), axis=1)
        fatW_idx = ak.argmin(self.p4_cache['genV2'].deltaR(self.p4_cache['SelectedFatJet']), axis=1) 

        match self.qq_topology:
            case 'resolved':
                bquark_cuts['resolved_qq'] = ~boosted_qq
                bquark_cuts['reco_resolved_qq'] = np.logical_and(reco_2q, df['ncentralJet'] >= 2)
            case 'boosted':
                bquark_cuts['boosted_bb'] = boosted_bb
                bquark_cuts['reco_boosted_bb'] = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
            case 'mixed':
                as_fat = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
                as_slim = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
                bquark_cuts['reco_bb'] = np.logical_or(as_fat, as_slim)
            case _:
                raise RuntimeError(f'Illegal use_topology value: {topology}')

        return light_quark_cuts
        

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

        match self.channel:
            case 'DL':   
                match self.bb_topology:
                    case 'resolved':
                        cuts['is_resolved'] = ~is_boosted
                        cuts['reco_as_resolved'] = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
                    case 'boosted':
                        cuts['is_boosted'] = is_boosted
                        cuts['reco_as_boosted'] = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
                    case 'mixed':
                        as_fat = np.logical_and(mindR_Hbb_FatJet < 0.8, df['nSelectedFatJet'] >= 1)
                        as_slim = np.logical_and(reco_2b, df['ncentralJet'] >= 2)
                        cuts['successful_b_jet_reco'] = np.logical_or(as_fat, as_slim)
                    case _:
                        raise RuntimeError(f'Illegal use_topology value: {topology}')

                reco_lep1 = df['lep1_legType'] == df['lep1_gen_kind']
                reco_lep2 = df['lep2_legType'] == df['lep2_gen_kind']
                cuts['successful_lep_reco'] = np.logical_and(reco_lep1, reco_lep2)
            case 'SL':
                match (self.bb_topology, self.qq_topology):
                    pass
            case _:
                raise RuntimeError(f'Illegal channel {self.channel}, only SL or DL are allowed')

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
                cuts['successful_b_jet_reco'] = np.logical_or(as_fat, as_slim)
            case _:
                raise RuntimeError(f'Illegal use_topology value: {topology}')

        # lepton matching
        match self.channel:
            case 'DL':
                reco_lep1 = df['lep1_legType'] == df['lep1_gen_kind']
                reco_lep2 = df['lep2_legType'] == df['lep2_gen_kind']
                cuts['successful_lep_reco'] = np.logical_and(reco_lep1, reco_lep2)
            case 'SL':
                # qq_dr = self.p4_cache['genV2prod1'].deltaR(self.p4_cache['genV2prod2'])
                # is_boosted = qq_dr < 0.8

                # mindR_q1_Jet = ak.min(self.p4_cache['genV2prod1'].deltaR(self.p4_cache['centralJet']), axis=1)
                # mindR_q2_Jet = ak.min(self.p4_cache['genV2prod2'].deltaR(self.p4_cache['centralJet']), axis=1)
                # q1_match_idx = ak.argmin(self.p4_cache['genV2prod1'].deltaR(self.p4_cache['centralJet']), axis=1)
                # q2_match_idx = ak.argmin(self.p4_cache['genV2prod2'].deltaR(self.p4_cache['centralJet']), axis=1)
                # has_two_matches = np.logical_and(mindR_q1_Jet < 0.4, mindR_q2_Jet < 0.4)
                # reco_2q = np.logical_and(has_two_matches, q1_match_idx != q2_match_idx)
                # mindR_HVV_FatJet = ak.min(self.p4_cache['genV2'].deltaR(self.p4_cache['SelectedFatJet']), axis=1)

                # as_fat = np.logical_and(mindR_HVV_FatJet < 0.8, df['nSelectedFatJet'] >= 2)
                # as_slim = np.logical_and(reco_2q, df['ncentralJet'] >= 4)
                # cuts['successful_light_jet_reco'] = np.logical_or(as_fat, as_slim)

                reco_lep1 = df['lep1_legType'] == df['lep1_gen_kind']
                cuts['successful_lep_reco'] = reco_lep1
            case _:
                raise RuntimeError(f'Illegal channel {self.channel}, only SL or DL are allowed')
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
            title = f"Cutflow {self.channel} {sample_type} M={mass}" if mp and mp > 0 else f"Cutflow  {self.channel} {sample_type}"
            ax.set(ylabel='Number of events passed', title=title)
            if self.loader_cfg['use_percentages']:
                ax.bar_label(bar_container, fmt=lambda x: f"{x:.1f}%", fontsize=8)
            else:
                ax.bar_label(bar_container, fmt=lambda x: int(x), fontsize=8)
            plt.xticks(rotation=45, ha='right')

            cutlfow_name = f"cutflow_{self.channel}_{sample_type}"
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

    def GetObjData(self, objects=None, as_df=True, unravel=False, selector=None, *args, **kwargs):
        """
        Args:
            objects: list of particles whose data to return
            as_df: bool controlling return format
            unravel: if true, if target shape > 1, arrays will be separated 
            filter: callable describing how to filter rows from df by event id
        if as_df is True returns dict obj_name -> df with this object's data
        else returns dict obj_name -> (obj_feature_names, obj_feature_arrays)
        """
        obj_dict = self.object_feature_names if objects is None else {obj: self.object_feature_names[obj] for obj in objects}
        df = self.df[selector(self.df, *args)] if selector else self.df

        if as_df and unravel:
            raise RuntimeError(f'Incompatible arguments: `as_df`={as_df} and `unravel`={unravel} cannot be true simultaneously')

        out = {}
        if as_df:
            out = {obj_name: df[obj_features] for obj_name, obj_features in obj_dict.items()}
        else:
            if unravel:
                for obj_name, obj_features in obj_dict.items():
                    shape = self.obj_shapes[obj_name]
                    if shape == 1:
                        out[obj_name] = (obj_features, df[obj_features].values)
                    else:
                        for i in range(shape):
                            # format of features is objname_id_feature
                            # for jets, there is centralJet_1_* and centralJet_10_*
                            # condition f'{obj_name}_{i + 1}' in f (without trailing underscore) is not working correctly
                            features = [f for f in obj_features if f'{obj_name}_{i + 1}_' in f]
                            out[f'{obj_name}_{i + 1}'] = (features, df[features].values) 
            else:
                out = {obj_name: (obj_features, df[obj_features].values) for obj_name, obj_features in obj_dict.items()}
        return out