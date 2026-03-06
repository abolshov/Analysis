import tensorflow as tf

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print(f"Enabled memory growth for {len(gpus)} GPU(s).")
    except RuntimeError as e:
        # Memory growth must be set before GPUs have been initialized
        print(e)

import psutil
import uproot
import pathlib
import awkward as ak
import numpy as np
import os
import re

from typing import List

def to_numpy(*,
             ak_array: ak.Array, 
             dtype=np.float64) -> np.ndarray:
    """
    Pre-allocate numpy array for better memory efficiency
    """
    fields = ak.fields(ak_array)
    M = len(ak_array)
    N = len(fields)
    
    result = np.empty((M, N), dtype=dtype)
    
    for i, field in enumerate(fields):
        result[:, i] = ak.to_numpy(ak_array[field])
    
    return result

def load_file(*,
              tree_name: str,
              file_path: str | os.PathLike | pathlib.Path,
              list_of_branches: List[str],
              convert_to_numpy: bool):

    file = uproot.open(file_path)
    tree = file[tree_name]
    branches = tree.arrays(list_of_branches)

    if convert_to_numpy:
        return to_numpy(ak_array=branches)
    
    return branches

def main():
    print("Training DNN for event classification")

    process = psutil.Process(os.getpid())
    memory_bytes = process.memory_info().rss
    memory_mb = memory_bytes / (1024 * 1024)
    print(f"Memory usage {memory_mb:.2f} MB")

    weight_file_pattern = re.compile(r"nParity(\d)_Merged_weight.root")
    event_file_pattern = re.compile(r"nParity(\d)_Merged.root")

    train_file_dir = "/home/artem/Desktop/CMS/data/DNN/SL/resolved/Run3_2022/Dataset"
    event_files = {}
    weight_files = {}
    for f in os.listdir(train_file_dir):
        weight_match = weight_file_pattern.search(f)
        event_match = event_file_pattern.search(f)
        abs_path = os.path.abspath(os.path.join(train_file_dir, f))
        if weight_match:
            parity = int(weight_match.group(1))
            weight_files[parity] = abs_path
        elif event_match:
            parity = int(event_match.group(1))
            event_files[parity] = abs_path
        else:
            continue

    parity_file_map = { p: (event_files[p], weight_files[p]) for p in event_files.keys()}
    for p, (event_file, weight_file) in parity_file_map.items():
        print(f"Parity {p}: event_file={event_file}, weight_file={weight_file}")

    train_parity = 0
    train_event_file_path, train_weight_file_path = parity_file_map[train_parity]

    use_extra_variables = False

    base_branches = [
        'other_jet1_mass', 
        'deta_bb', 
        'dPhi_jet1_jet2', 
        'other_jet1_btagPNetB', 
        'nSelBtag_fatjets',  
        'dPhi_MET_dibjet',
        'fatbjet_eta', 
        'light_jet2_mass', 
        'other_jet2_mass',
        'other_jet2_pt', 
        'other_jet1_pt', 
        'lep1_mass', 
        'bjet2_btagPNetB', 
        'other_jet2_eta', 
        'bjet1_phi', 
        'MT2_bb',  
        'bjet2_pt', 
        'fatbjet_phi',
        'other_jet2_btagPNetB', 
        'bjet1_mass', 
        'dphi_bb', 
        'bjet2_mass', 
        'PuppiMET_pt', 
        'light_jet2_phi',
        'other_jet2_phi', 
        'bjet1_eta', 
        'fatbjet_mass', 
        'bb_mass_PNetRegPtRawCorr_PNetRegPtRawCorrNeutrino', 
        'light_jet1_phi', 
        'other_jet1_eta', 
        'fatbjet_particleNetWithMass_HbbvsQCD', 
        'SingleLep_DeepHME_mass', 
        'nSelBtag_jets', 
        'Lep1Jet1Jet2_mass', 
        'bjet2_eta', 
        'PuppiMET_phi', 
        'fatbjet_particleNet_XbbVsQCD', 
        'fatbjet_mass_PNetCorr', 
        'SingleLep_DeepHME_mass_error', 
        'CosTheta_bb', 
        'bjet1_pt', 
        'MT', 
        'fatbjet_pt', 
        'bjet2_phi', 
        'fatbjet_msoftdrop', 
        'lep1_phi', 
        'dR_dibjet', 
        'HT', 
        'lep1_pt', 
        'bjet1_btagPNetB', 
    ]

    extra_branches = [
        'dphi_hadT_hadW', 
        'fatbjet_muEF', 
        'min_dphi_b_hadW', 
        'dphi_Hbb_Hww', 
        'mT_Hww', 
        'min_dphi_bl', 
        'lep2_pt', 
        'm_Hbb', 
        'lep2_mass', 
        'fatbjet_neEmEF', 
        'pT_Hbb', 
        'min_dphi_b_lepW', 
        'mass_hadT', 
        'light_jet2_pt', 
        'light_jet2_btagPNetB', 
        'light_jet1_mass',
        'light_jet1_eta', 
        'pT_lepT', 
        'fatbjet_neMultiplicity',
        'min_dR_b_hadW', 
        'ratio_pT_hadT_const', 
        'fatbjet_tau3', 
        'light_jet1_btagPNetB', 
        'pTtoE_ratio_Hww', 
        'event', 
        'm_Hww',
        'dphi_hadW_lepW', 
        'm_hadW',
        'dR_hadW_lepW', 
        'dR_Hbb_Hww', 
        'light_jet1_pt', 
        'fatbjet_tau1', 
        'mT_lepT', 
        'fatbjet_tau2', 
        'mT_hadT', 
        'deta_leadBjet_hadW', 
        'pT_lepW', 
        'lep2_eta', 
        'dR_hadW_l', 
        'm_b_lepT',
        'min_dR_b_lepW', 
        'pT_hadW', 
        'deta_hadT_hadW', 
        'dphi_hadW_l', 
        'min_dR_bl',  
        'mT_lepW', 
        'lep2_phi', 
        'fatbjet_neHEF', 
        'pT_Hww', 
        'pT_hadT',
        'fatbjet_tau4', 
        'dR_leadBjet_hadW', 
        'pTtoE_ratio_HH', 
        'fatbjet_nConstituents', 
        'light_jet2_eta', 
        'dphi_leadBjet_hadW', 
        'pTtoE_ratio_Hbb',
        'dR_hadT_hadW', 
        'mT_hadW'
    ]

    branches_to_load = base_branches
    if use_extra_variables:
        branches_to_load.extend(extra_branches)

    X_train = load_file(tree_name="Events", 
                        file_path=train_event_file_path, 
                        list_of_branches=branches_to_load,
                        convert_to_numpy=True)

    weights_train = load_file(tree_name="weight_tree", 
                        file_path=train_weight_file_path, 
                        list_of_branches=["class_weight", "class_target"],
                        convert_to_numpy=False)

    train_weights = weights_train["class_weight"].to_numpy()
    y_train = weights_train["class_target"].to_numpy()

    memory_bytes = process.memory_info().rss
    memory_mb = memory_bytes / (1024 * 1024)
    print(f"Loaded {X_train.shape[1]} variables for {X_train.shape[0]} events for training set")
    print(f"Memory usage {memory_mb:.2f} MB")

    val_parity = 0
    val_event_file_path, val_weight_file_path = parity_file_map[val_parity]
    X_val = load_file(tree_name="Events", 
                      file_path=val_event_file_path, 
                      list_of_branches=branches_to_load,
                      convert_to_numpy=True)
    
    y_val = load_file(tree_name="weight_tree", 
                      file_path=val_weight_file_path, 
                      list_of_branches=["class_target"],
                      convert_to_numpy=True)
    
    memory_bytes = process.memory_info().rss
    memory_mb = memory_bytes / (1024 * 1024)
    print(f"Loaded {X_val.shape[1]} variables for {X_val.shape[0]} events for validation set")
    print(f"Memory usage {memory_mb:.2f} MB")

if __name__ == "__main__":
    main()