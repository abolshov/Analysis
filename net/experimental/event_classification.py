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
import awkward as ak
import numpy as np
import os
import re

def main():
    print("Training DNN for event classification")

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
    train_event_file = uproot.open(train_event_file_path)
    train_weight_file = uproot.open(train_weight_file_path)

    event_tree_name = "Events"
    weight_tree_name = "weight_tree"
    event_tree = train_event_file[event_tree_name]
    weight_tree = train_weight_file[weight_tree_name]

    branches_to_load = [
        'other_jet1_mass', 
        'dphi_hadT_hadW', 
        'fatbjet_muEF', 
        'min_dphi_b_hadW', 
        'deta_bb', 
        'dPhi_jet1_jet2', 
        'other_jet1_btagPNetB', 
        'nSelBtag_fatjets', 
        'dphi_Hbb_Hww', 
        'mT_Hww', 
        'dPhi_MET_dibjet',
        'fatbjet_eta', 
        'min_dphi_bl', 
        'lep2_pt', 
        'fatbjet_particleNet_XbbVsQCD_bin', 
        'light_jet2_mass', 
        'other_jet2_mass',
        'm_Hbb', 
        'lep2_mass', 
        'fatbjet_neEmEF', 
        'pT_Hbb', 
        'min_dphi_b_lepW', 
        'other_jet2_pt', 
        'other_jet1_pt', 
        'lep1_mass', 
        'bjet2_btagPNetB', 
        'mass_hadT', 
        'other_jet2_eta', 
        'light_jet2_pt', 
        'light_jet2_btagPNetB', 
        'light_jet1_mass',
        'bjet1_phi', 
        'MT2_bb',  
        'light_jet1_eta', 
        'pT_lepT', 
        'bjet2_pt', 
        'fatbjet_neMultiplicity',
        'min_dR_b_hadW', 
        'ratio_pT_hadT_const', 
        'fatbjet_tau3', 
        'light_jet1_btagPNetB', 
        'pTtoE_ratio_Hww', 
        'fatbjet_phi',
        'other_jet2_btagPNetB', 
        'event', 
        'bjet1_mass', 
        'm_Hww',
        'dphi_hadW_lepW', 
        'm_hadW',
        'dR_hadW_lepW', 
        'dR_Hbb_Hww', 
        'dphi_bb', 
        'light_jet1_pt', 
        'bjet2_mass', 
        'fatbjet_tau1', 
        'PuppiMET_pt', 
        'light_jet2_phi',
        'other_jet2_phi', 
        'bjet1_eta', 
        'mT_lepT', 
        'fatbjet_mass', 
        'fatbjet_tau2', 
        'mT_hadT', 
        'bb_mass_PNetRegPtRawCorr_PNetRegPtRawCorrNeutrino', 
        'light_jet1_phi', 
        'other_jet1_eta', 
        'deta_leadBjet_hadW', 
        'fatbjet_particleNetWithMass_HbbvsQCD', 
        'pT_lepW', 
        'SingleLep_DeepHME_mass', 
        'lep2_eta', 
        'dR_hadW_l', 
        'm_b_lepT',
        'nSelBtag_jets', 
        'min_dR_b_lepW', 
        'Lep1Jet1Jet2_mass', 
        'pT_hadW', 
        'deta_hadT_hadW', 
        'bjet2_eta', 
        'dphi_hadW_l', 
        'PuppiMET_phi', 
        'min_dR_bl', 
        'fatbjet_particleNet_XbbVsQCD', 
        'fatbjet_mass_PNetCorr', 
        'SingleLep_DeepHME_mass_error', 
        'CosTheta_bb', 
        'mT_lepW', 
        'lep2_phi', 
        'fatbjet_neHEF', 
        'pT_Hww', 
        'bjet1_pt', 
        'MT', 
        'fatbjet_pt', 
        'bjet2_phi', 
        'fatbjet_msoftdrop', 
        'lep1_phi', 
        'pT_hadT',
        'fatbjet_tau4', 
        'dR_dibjet', 
        'dR_leadBjet_hadW', 
        'pTtoE_ratio_HH', 
        'fatbjet_nConstituents', 
        'HT', 
        'light_jet2_eta', 
        'lep1_pt', 
        'dphi_leadBjet_hadW', 
        'bjet1_btagPNetB', 
        'pTtoE_ratio_Hbb',
        'dR_hadT_hadW', 
        'mT_hadW'
    ]

    weight_branches = weight_tree.arrays()
    event_branches = event_tree.arrays(branches_to_load)

    process = psutil.Process(os.getpid())
    memory_bytes = process.memory_info().rss
    memory_mb = memory_bytes / (1024 * 1024)
    print(f"Loaded {len(event_branches)} events, current memory usage {memory_mb:.2f} MB")

if __name__ == "__main__":
    main()