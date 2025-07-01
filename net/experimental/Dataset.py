import numpy as np
import awkward as ak
import pandas as pd
import uproot
import yaml
import vector
import matplotlib.pyplot as plt
import os


sample_type_map = {1: "GluGluToRadion",
                   2: "GluGluToBulkGraviton", 
                   3 : "VBFToRadion", 
                   4 : "VBFToBulkGraviton"}


class Dataset:
    def __init__(self, branch_loading_cfg):
        self.dataset_cfg = None
        self.df = pd.DataFrame()
        self.p4_cache = {}
        with open (branch_loading_cfg, 'r') as branch_loading_cfg_file:
            self.dataset_cfg = yaml.safe_load(branch_loading_cfg_file)
        assert self.dataset_cfg, "Must contain mapping of branches to load for each object"

    def MakeDataFrame(self, tree):
        branch_collections = self.dataset_cfg['branches_to_load']
        branch_dict = {}

        if 'X_mass' in tree.keys():
            branch_dict['masspoint'] = tree['X_mass'].array().to_numpy()

        for obj_coll in branch_collections.keys():
            branches_cfg = branch_collections[obj_coll]
            branch_list = branches_cfg['names']
            arrays = tree.arrays(branch_list)
            shape = branches_cfg['shape']

            momentum_branches = {}
            for branch in branch_list:
                min_depth, max_depth = arrays[branch].layout.minmax_depth
                branch_var = branch.split('_')[-1]
                if branch_var in ['pt', 'eta', 'phi', 'mass']:
                    particle = branch.split('_')[0]
                    if particle not in momentum_branches:
                        momentum_branches[particle] = {}
                    if max_depth > 1:
                        momentum_branches[particle][branch_var] = ak.fill_none(ak.pad_none(arrays[branch], shape), 0.0)
                    else:
                        momentum_branches[particle][branch_var] = arrays[branch].to_numpy()
                else:
                    if max_depth > 1:
                        arr = ak.fill_none(ak.pad_none(arrays[branch], shape), 0)
                        arr = arr[:, :shape]
                        for idx in range(shape):
                            branch_dict[f'{branch}_{idx}'] = arr[:, idx].to_numpy()
                    else:
                        branch_dict[branch] = arrays[branch].to_numpy()

            if momentum_branches:
                for particle in momentum_branches.keys():
                    p4 = vector.zip(momentum_branches[particle])
                    self.p4_cache[particle] = p4
                    if shape > 1:
                        p4 = p4[:, :shape]
                    # met p4 only has 2 keys - pt and phi
                    if len(momentum_branches[particle].keys()) > 2:
                        for new_var in ['px', 'py', 'pz', 'E']:
                            if shape > 1:
                                for idx in range(shape):
                                    branch_dict[f'{particle}_{new_var}_{idx}'] = (lambda x: getattr(x, new_var))(p4[:, idx]).to_numpy()
                            else:
                                branch_dict[f'{particle}_{new_var}'] = (lambda x: getattr(x, new_var))(p4).to_numpy()
                    else:
                        for new_var in ['px', 'py']:
                            branch_dict[f'{particle}_{new_var}'] = (lambda x: getattr(x, new_var))(p4).to_numpy()
    
        # compute high-level kinematics needed for selections using p4_cache
        branch_dict['b1_accept_pt'] = self.p4_cache['genb1'].pt > 20
        branch_dict['b2_accept_pt'] = self.p4_cache['genb2'].pt > 20
        branch_dict['b1_accept_eta'] = np.abs(self.p4_cache['genb1'].eta) < 2.5
        branch_dict['b2_accept_eta'] = np.abs(self.p4_cache['genb2'].eta) < 2.5
        branch_dict['genb1_accept'] = np.logical_and(branch_dict['b1_accept_pt'], branch_dict['b1_accept_eta'])
        branch_dict['genb2_accept'] = np.logical_and(branch_dict['b2_accept_pt'], branch_dict['b2_accept_eta'])
        branch_dict['genb_accept'] = np.logical_and(branch_dict['genb1_accept'], branch_dict['genb2_accept'])
        branch_dict['bb_dr'] = self.p4_cache['genb1'].deltaR(self.p4_cache['genb2'])
        branch_dict['mindR_b1_Jet'] = ak.min(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        branch_dict['mindR_b2_Jet'] = ak.min(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        branch_dict['b1_match_idx'] = ak.argmin(self.p4_cache['genb1'].deltaR(self.p4_cache['centralJet']), axis=1)
        branch_dict['b2_match_idx'] = ak.argmin(self.p4_cache['genb2'].deltaR(self.p4_cache['centralJet']), axis=1)
        branch_dict['has_two_matches'] = np.logical_and(branch_dict['mindR_b1_Jet'] < 0.4, branch_dict['mindR_b2_Jet'] < 0.4)
        # branch_dict['one_to_one_match'] = np.logical_and(branch_dict['has_two_matches'], branch_dict['b1_match_idx'] != branch_dict['b2_match_idx'])
        branch_dict['mindR_Hbb_FatJet'] = ak.min(self.p4_cache['genHbb'].deltaR(self.p4_cache['SelectedFatJet']), axis=1)
        branch_dict['is_boosted'] = branch_dict['bb_dr'] < 0.8
        branch_dict['is_resolved'] = ~branch_dict['is_boosted']
        branch_dict['reco_as_boosted'] = np.logical_and(branch_dict['mindR_Hbb_FatJet'] < 0.8, branch_dict['nSelectedFatJet'] >= 1)
        branch_dict['reco_2b'] = np.logical_and(branch_dict['has_two_matches'], branch_dict['b1_match_idx'] != branch_dict['b2_match_idx'])
        branch_dict['reco_as_resolved'] = np.logical_and(branch_dict['reco_2b'], branch_dict['ncentralJet'] >= 2)
        branch_dict['successful_jet_reco'] = np.logical_or(branch_dict['reco_as_resolved'], branch_dict['reco_as_boosted'])
        branch_dict['reco_lep1'] = np.logical_and(self.p4_cache['lep1'].pt > 0.0, branch_dict['lep1_legType'] == branch_dict['lep1_gen_kind'])
        branch_dict['reco_lep2'] = np.logical_and(self.p4_cache['lep2'].pt > 0.0, branch_dict['lep2_legType'] == branch_dict['lep2_gen_kind'])
        branch_dict['successful_lep_reco'] = np.logical_and(branch_dict['reco_lep1'], branch_dict['reco_lep2'])
        branch_dict['successful_reco'] = np.logical_and(branch_dict['successful_lep_reco'], branch_dict['successful_jet_reco'])
        
        return pd.DataFrame.from_dict(branch_dict)

    def Load(self):
        for file_name in self.dataset_cfg['file_names']:
            tree_name = self.dataset_cfg['tree_name']
            tree = uproot.open(f'{file_name}:{tree_name}')

            df = self.MakeDataFrame(tree)
            if self.dataset_cfg['apply_selections']:
                df = self.ApplySelections(df, self.dataset_cfg['plot_cutflow'])
                
            self.df = pd.concat([self.df, df], axis=0)

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

            if self.dataset_cfg['use_topology'] == 'boosted':
                df = df[df['is_boosted']]
                cutflow['is_boosted'] = len(df)
                df = df[df['reco_lep1']]
                cutflow['reco_lep1'] = len(df)
                df = df[df['reco_lep2']]
                cutflow['reco_lep2'] = len(df)
                df = df[df['reco_as_boosted']]
                cutflow['reco_as_boosted'] = len(df)
            elif self.dataset_cfg['use_topology'] == 'resolved':
                df = df[df['is_resolved']]
                cutflow['is_resolved'] = len(df)
                df = df[df['reco_lep1']]
                cutflow['reco_lep1'] = len(df)
                df = df[df['reco_lep2']]
                cutflow['reco_lep2'] = len(df)
                df = df[df['reco_as_resolved']]
                cutflow['reco_as_resolved'] = len(df)
            elif self.dataset_cfg['use_topology'] == 'mixed':
                df = df[df['successful_lep_reco']]
                cutflow['successful_lep_reco'] = len(df)
                df = df[df['successful_jet_reco']]
                cutflow['successful_jet_reco'] = len(df)
            else:
                raise RuntimeError(f"Got unknown topology {self.dataset_cfg['use_topology']}")

            cut_names = cutflow.keys()
            passed_evt_info = list(cutflow.values())
            if self.dataset_cfg['use_percentages']:
                passed_evt_info = [val/passed_evt_info[0]*100 for val in passed_evt_info]
            fig, ax = plt.subplots()
            bar_container = ax.bar(cut_names, passed_evt_info)
            sample_type = sample_type_map[df['sample_type'].iloc[0]]
            mp = df['masspoint'].iloc[0]
            mass = int(mp) if not np.isnan(mp) else None
            title = f"Cutflow {sample_type} M={mass}" if mp else f"Cutflow {sample_type}"
            ax.set(ylabel='Number of events passed', title=title)
            if self.dataset_cfg['use_percentages']:
                ax.bar_label(bar_container, fmt=lambda x: f"{x:.1f}%", fontsize=8)
            else:
                ax.bar_label(bar_container, fmt=lambda x: int(x), fontsize=8)
            plt.xticks(rotation=45, ha='right')

            plt.savefig("cutflow.pdf", format='pdf', bbox_inches='tight')
            plt.close()

            return df

                
