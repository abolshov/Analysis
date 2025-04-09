import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
import argparse
import os
from itertools import groupby
import re


def MapIterTrees(iter_file_name):
    """
    maps event_id to iteration data
    parameters: name of file with iteration data
    returns: dictionary (event_id, list of trees with iter data for this event)
    """
    with uproot.open(iter_file_name) as iter_file:
        iter_tree_names = iter_file.keys()
        
        tree_name_dict = {}
        for event_id, group in groupby(iter_tree_names, lambda tree_name: int(tree_name.split('_')[1])):
            tree_name_dict[event_id] = [iter_tree_name for iter_tree_name in group]
        
        return tree_name_dict


def ExtractJetIndices(tree_name, sort=True):
    """
    extracts indices of jets from corresponding tree name
    parameters: name of tree (string), flag indicating whether sort indices or not
    returns: sorted list of indices
    """
    # problem: this function cannot be used in SL because in SL jet combination contains 4 indices
    jet_comb = re.findall(r'b\db\d', tree_name)
    if len(jet_comb) > 1:
        raise RuntimeError(f"Obtained impossible jet combinations {jet_comb} when parsing tree name {tree_name}")
    return sorted([int(c) for c in re.findall(r'\d+', jet_comb[0])]) if sort else [int(c) for c in re.findall(r'\d+', jet_comb[0])]


def ComputePeakValue(data, weights=None):
    """
    peak value from computing histogram of the input array and extracting bin with max count
    parameters: 1d array with data and optional 1d array of weights
    returns: x-position of the center of the bin with maximal count
    """
    counts, bin_edges = np.histogram(data, bins=1800, range=(200, 2000), weights=weights)
    peak_index = np.argmax(counts)
    peak_bin_edge_left = bin_edges[peak_index]
    peak_bin_edge_right = bin_edges[peak_index + 1]
    return (peak_bin_edge_right + peak_bin_edge_left)/2


def main():
    # make directory for plotting
    plotting_dir = "plotting_dl_M800"
    os.makedirs(plotting_dir, exist_ok=True)

    hme_out_true_name = "hme_out_true_dl_M800.root" # hme output for reco data matched to partons
    hme_out_name = "hme_out_dl_M800.root" # hme ran on everything

    hme_true_out_tree = uproot.open(f"{hme_out_true_name}:hme_tree")
    hme_out_tree = uproot.open(f"{hme_out_name}:hme_tree")

    hme_true_mass = hme_true_out_tree["hme_mass"].array()
    hme_mass = hme_out_tree["hme_mass"].array()

    hme_true_event_id = hme_true_out_tree["event_id"].array()
    hme_event_id = hme_out_tree["event_id"].array()

    # find indices of events such that true_evt_id is in list of all event ids hme_event_id
    common_evt_id, hme_true_idx, hme_idx = np.intersect1d(hme_true_event_id, hme_event_id, return_indices=True)

    plt.hist(hme_true_mass, bins=100, range=(200, 2000), histtype='step', label="matched inputs")
    plt.hist(hme_mass[hme_idx], bins=100, range=(200, 2000), histtype='step', label="any inputs")
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.title(f"HME mass")
    plt.xlabel("mX, [GeV]")
    plt.savefig(os.path.join(plotting_dir, "hme.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    true_iter_file_name = "hme_iter_true_dl_M800.root"
    true_tree_name_dict = MapIterTrees(true_iter_file_name)

    iter_file_name = "hme_iter_dl_M800.root"
    tree_name_dict = MapIterTrees(iter_file_name)
    iter_file = uproot.open(iter_file_name)
    
    # nanoAOD file used to get HME
    nano_file_name = "nano_dl_M800_v2.root"
    nano_tree = uproot.open(f"{nano_file_name}:Events")
    nano_branches = nano_tree.arrays()

    for event_id, tree_name_list in tree_name_dict.items():
        evt_plt_dir = f"event_{event_id}"
        os.makedirs(os.path.join(plotting_dir, evt_plt_dir), exist_ok=True)

        # in the map event_id -> true_tree_name find name of the true tree and extract indices
        true_tree_name = true_tree_name_dict[event_id][0] if true_tree_name_dict[event_id] else None
        if not true_tree_name:
            print(f"In event {event_id} reconstructed jets do not match MC quarks")
            continue

        match_indices = ExtractJetIndices(true_tree_name)
        match_b1_idx, match_b2_idx = ExtractJetIndices(true_tree_name, sort=False)

        match_jet1_pt = nano_branches[nano_branches["event"] == event_id]["centralJet_pt"][0][match_b1_idx]
        match_jet2_pt = nano_branches[nano_branches["event"] == event_id]["centralJet_pt"][0][match_b2_idx]

        quark1_pt = nano_branches[nano_branches["event"] == event_id]["genb1_pt"][0]
        quark2_pt = nano_branches[nano_branches["event"] == event_id]["genb2_pt"][0]

        lead_quark_pt = max(quark1_pt, quark2_pt)
        sublead_quark_pt = min(quark1_pt, quark2_pt)

        lead_reco_pt = max(match_jet1_pt, match_jet2_pt)
        sublead_reco_pt = min(match_jet1_pt, match_jet2_pt)

        nano_hbb_pt = nano_branches[nano_branches["event"] == event_id]["genHbb_pt"][0]
        nano_hbb_eta = nano_branches[nano_branches["event"] == event_id]["genHbb_eta"][0]
        nano_hbb_phi = nano_branches[nano_branches["event"] == event_id]["genHbb_phi"][0]

        nano_hww_pt = nano_branches[nano_branches["event"] == event_id]["genHVV_pt"][0]
        nano_hww_eta = nano_branches[nano_branches["event"] == event_id]["genHVV_eta"][0]
        nano_hww_phi = nano_branches[nano_branches["event"] == event_id]["genHVV_phi"][0]

        nano_nu1_pt = nano_branches[nano_branches["event"] == event_id]["genV1prod2_pt"][0]
        nano_nu1_eta = nano_branches[nano_branches["event"] == event_id]["genV1prod2_eta"][0]
        nano_nu1_phi = nano_branches[nano_branches["event"] == event_id]["genV1prod2_phi"][0]

        nano_nu2_pt = nano_branches[nano_branches["event"] == event_id]["genV2prod2_pt"][0]
        nano_nu2_eta = nano_branches[nano_branches["event"] == event_id]["genV2prod2_eta"][0]
        nano_nu2_phi = nano_branches[nano_branches["event"] == event_id]["genV2prod2_phi"][0]

        cumulative_mass = []
        cumulative_weights = []

        comb_mass = []
        comb_label = []
        comb_weight = []

        for tree_idx, tree_name in enumerate(tree_name_list):
            jet_indices = ExtractJetIndices(tree_name)
            comb_tag = tree_name.split('_')[-1][:-2] 
            label = comb_tag
            if jet_indices == match_indices:
                label += ": true"

            tree = iter_file[tree_name]
            iter_branches = tree.arrays()

            # plotting comparisons of jet passed to HME vs b quark vs reco jet
            bjet1_pt = iter_branches["bjet_1_pt"].to_numpy()
            bjet2_pt = iter_branches["bjet_2_pt"].to_numpy()

            counts, edges = np.histogram(bjet1_pt, bins=50, range=(0, 500))
            height = np.max(counts)
            plt.hist(bjet1_pt, bins=50, range=(0, 500), histtype='step', label="leading HME jet")
            plt.hist(lead_reco_pt*np.ones(height), bins=50, range=(0, 500), label="leading reco jet")
            plt.hist(lead_quark_pt*np.ones(height), bins=50, range=(0, 500), label="gen quark")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} leading jet pt")
            plt.xlabel("pt, [GeV]")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_lead_jet_pt.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            counts, edges = np.histogram(bjet2_pt, bins=50, range=(0, 500))
            height = np.max(counts)
            plt.hist(bjet2_pt, bins=50, range=(0, 500), histtype='step', label="subleading HME jet")
            plt.hist(sublead_reco_pt*np.ones(height), bins=50, range=(0, 500), label="subleading reco jet")
            plt.hist(sublead_quark_pt*np.ones(height), bins=50, range=(0, 500), label="gen quark")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} subleading jet pt")
            plt.xlabel("pt, [GeV]")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_sublead_jet_pt.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            # plotting HME mass for different combinations
            # branch mass is array of lenght 4
            mass = iter_branches["mass"]
            # branch weight is single float
            # need to transform it to shape (n_iters, 4) = mass.shape
            weights = iter_branches["weight"]
            weights = ak.unflatten(weights, 1).to_numpy()
            weights = np.repeat(weights, 4, axis=1).flatten()
            mass = ak.flatten(mass).to_numpy()

            succesful_iters = mass > 0.0
            mass = mass[succesful_iters]
            weights = weights[succesful_iters]

            cumulative_mass.extend(mass)
            cumulative_weights.extend(weights)

            estimated_mass = ComputePeakValue(mass, weights=weights)
            if jet_indices == match_indices:
                label += f", m={estimated_mass:.2f}"
            else :
                label += f": m={estimated_mass:.2f}"
            
            comb_mass.append(mass)
            comb_weight.append(weights)
            comb_label.append(label)

            # H->bb comparisons
            hbb_pt = iter_branches["Hbb_pt"].to_numpy()
            counts, edges = np.histogram(hbb_pt, bins=50, range=(0, 500))
            height = np.max(counts)
            plt.hist(hbb_pt, bins=50, range=(0, 500), histtype='step', label="HME H->bb")
            plt.hist(nano_hbb_pt*np.ones(height), bins=50, range=(0, 500), label="gen H->bb")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} H->bb pt")
            plt.xlabel("pt, [GeV]")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_hbb_pt.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            hbb_eta = iter_branches["Hbb_eta"].to_numpy()
            counts, edges = np.histogram(hbb_eta, bins=50, range=(-6, 6))
            height = np.max(counts)
            plt.hist(hbb_eta, bins=50, range=(-6, 6), histtype='step', label="HME H->bb")
            plt.hist(nano_hbb_eta*np.ones(height), bins=50, range=(-6, 6), label="gen H->bb")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} H->bb eta")
            plt.xlabel("eta")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_hbb_eta.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            hbb_phi = iter_branches["Hbb_phi"].to_numpy()
            counts, edges = np.histogram(hbb_phi, bins=50, range=(-4, 4))
            height = np.max(counts)
            plt.hist(hbb_phi, bins=50, range=(-4, 4), histtype='step', label="HME H->bb")
            plt.hist(nano_hbb_phi*np.ones(height), bins=50, range=(-4, 4), label="gen H->bb")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} H->bb phi")
            plt.xlabel("phi")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_hbb_phi.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            # H->WW comparisons
            # dump all iterations in single histogram
            hww_pt = ak.flatten(iter_branches["Hww_pt"]).to_numpy()
            hww_pt = hww_pt[succesful_iters]
            counts, edges = np.histogram(hww_pt, weights=weights, bins=50, range=(0, 500))
            height = np.max(counts)
            plt.hist(hww_pt, weights=weights, bins=50, range=(0, 500), histtype='step', label="HME H->WW")
            plt.hist(nano_hww_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen H->WW")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} H->WW pt")
            plt.xlabel("pt, [GeV]")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_hww_pt.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            hww_eta = ak.flatten(iter_branches["Hww_eta"]).to_numpy()
            hww_eta = hww_eta[succesful_iters]
            counts, edges = np.histogram(hww_eta, weights=weights, bins=50, range=(-6, 6))
            height = np.max(counts)
            plt.hist(hww_eta, weights=weights, bins=50, range=(-6, 6), histtype='step', label="HME H->WW")
            plt.hist(nano_hww_eta*np.ones(int(height)), bins=50, range=(-6, 6), label="gen H->WW")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} H->WW eta")
            plt.xlabel("eta")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_hww_eta.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            hww_phi = ak.flatten(iter_branches["Hww_phi"]).to_numpy()
            hww_phi = hww_phi[succesful_iters]
            counts, edges = np.histogram(hww_phi, weights=weights, bins=50, range=(-4, 4))
            height = np.max(counts)
            plt.hist(hww_phi, weights=weights, bins=50, range=(-4, 4), histtype='step', label="HME H->WW")
            plt.hist(nano_hww_phi*np.ones(int(height)), bins=50, range=(-4, 4), label="gen H->WW")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} H->WW phi")
            plt.xlabel("phi")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_hww_phi.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            # # split H->WW into categories by control index now 
            # hww_pt = ak.flatten(iter_branches["Hww_pt"]).to_numpy()
            # weights = iter_branches["weight"]
            # weights = ak.unflatten(weights, 1).to_numpy()
            # weights = np.repeat(weights, 4, axis=1).flatten()
            # control_mask = np.arange(len(hww_pt)) % 4
            # from_onshellW = np.logical_and(succesful_iters, control_mask // 2 == 0)
            # from_offshellW = np.logical_and(succesful_iters, control_mask // 2 == 1)

            # counts, edges = np.histogram(hww_pt[from_onshellW], weights=weights[from_onshellW], bins=50, range=(0, 500))
            # height = np.max(counts)
            # plt.hist(hww_pt[from_onshellW], weights=weights[from_onshellW], bins=50, range=(0, 500), histtype='step', label="HME H->WW")
            # plt.hist(nano_hww_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen H->WW")
            # plt.legend(loc='upper right')
            # plt.grid(True)
            # plt.title(f"Event {event_id} comb {comb_tag} H->WW pt")
            # plt.xlabel("pt, [GeV]")
            # plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_from_onshell_hww_pt.pdf"), format='pdf', bbox_inches='tight')
            # plt.close()

            # counts, edges = np.histogram(hww_pt[from_offshellW], weights=weights[from_offshellW], bins=50, range=(0, 500))
            # height = np.max(counts)
            # plt.hist(hww_pt[from_offshellW], weights=weights[from_offshellW], bins=50, range=(0, 500), histtype='step', label="HME H->WW")
            # plt.hist(nano_hww_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen H->WW")
            # plt.legend(loc='upper right')
            # plt.grid(True)
            # plt.title(f"Event {event_id} comb {comb_tag} H->WW pt")
            # plt.xlabel("pt, [GeV]")
            # plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_from_offshell_hww_pt.pdf"), format='pdf', bbox_inches='tight')
            # plt.close()

            # neutrino comparisons
            # compare both gen neutrinos with nu_offshell
            nu_offshell_pt = ak.flatten(iter_branches["nu_offshell_pt"]).to_numpy()
            nu_offshell_pt = nu_offshell_pt[succesful_iters]
            counts, edges = np.histogram(nu_offshell_pt, weights=weights, bins=50, range=(0, 500))
            height = np.max(counts)
            plt.hist(nu_offshell_pt, weights=weights, bins=50, range=(0, 500), histtype='step', label="HME nu offshell")
            plt.hist(nano_nu1_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen nu1")
            plt.hist(nano_nu2_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen nu2")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} nu from offshell W pt")
            plt.xlabel("pt, [GeV]")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_nu_offshell_pt.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            nu_offshell_eta = ak.flatten(iter_branches["nu_offshell_eta"]).to_numpy()
            nu_offshell_eta = nu_offshell_eta[succesful_iters]
            counts, edges = np.histogram(nu_offshell_eta, weights=weights, bins=50, range=(-6, 6))
            height = np.max(counts)
            plt.hist(nu_offshell_eta, weights=weights, bins=50, range=(-6, 6), histtype='step', label="HME nu offshell")
            plt.hist(nano_nu1_eta*np.ones(int(height)), bins=50, range=(-6, 6), label="gen nu1")
            plt.hist(nano_nu2_eta*np.ones(int(height)), bins=50, range=(-6, 6), label="gen nu2")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} nu from offshell W eta")
            plt.xlabel("eta")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_nu_offshell_eta.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            nu_offshell_phi = ak.flatten(iter_branches["nu_offshell_phi"]).to_numpy()
            nu_offshell_phi = nu_offshell_phi[succesful_iters]
            counts, edges = np.histogram(nu_offshell_phi, weights=weights, bins=50, range=(-4, 4))
            height = np.max(counts)
            plt.hist(nu_offshell_phi, weights=weights, bins=50, range=(-4, 4), histtype='step', label="HME nu offshell")
            plt.hist(nano_nu1_phi*np.ones(int(height)), bins=50, range=(-4, 4), label="gen nu1")
            plt.hist(nano_nu2_phi*np.ones(int(height)), bins=50, range=(-4, 4), label="gen nu2")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} nu from offshell W phi")
            plt.xlabel("phi")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_nu_offshell_phi.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            # compare both gen neutrinos with nu_onshell
            nu_onshell_pt = ak.flatten(iter_branches["nu_onshell_pt"]).to_numpy()
            nu_onshell_pt = nu_onshell_pt[succesful_iters]
            counts, edges = np.histogram(nu_onshell_pt, weights=weights, bins=50, range=(0, 500))
            height = np.max(counts)
            plt.hist(nu_onshell_pt, weights=weights, bins=50, range=(0, 500), histtype='step', label="HME nu onshell")
            plt.hist(nano_nu1_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen nu1")
            plt.hist(nano_nu2_pt*np.ones(int(height)), bins=50, range=(0, 500), label="gen nu2")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} nu from onshell W pt")
            plt.xlabel("pt, [GeV]")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_nu_onshell_pt.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            nu_onshell_eta = ak.flatten(iter_branches["nu_onshell_eta"]).to_numpy()
            nu_onshell_eta = nu_onshell_eta[succesful_iters]
            counts, edges = np.histogram(nu_onshell_eta, weights=weights, bins=50, range=(-6, 6))
            height = np.max(counts)
            plt.hist(nu_onshell_eta, weights=weights, bins=50, range=(-6, 6), histtype='step', label="HME nu onshell")
            plt.hist(nano_nu1_eta*np.ones(int(height)), bins=50, range=(-6, 6), label="gen nu1")
            plt.hist(nano_nu2_eta*np.ones(int(height)), bins=50, range=(-6, 6), label="gen nu2")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} nu from onshell W eta")
            plt.xlabel("eta")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_nu_onshell_eta.pdf"), format='pdf', bbox_inches='tight')
            plt.close()

            nu_onshell_phi = ak.flatten(iter_branches["nu_onshell_phi"]).to_numpy()
            nu_onshell_phi = nu_onshell_phi[succesful_iters]
            counts, edges = np.histogram(nu_offshell_phi, weights=weights, bins=50, range=(-4, 4))
            height = np.max(counts)
            plt.hist(nu_onshell_phi, weights=weights, bins=50, range=(-4, 4), histtype='step', label="HME nu onshell")
            plt.hist(nano_nu1_phi*np.ones(int(height)), bins=50, range=(-4, 4), label="gen nu1")
            plt.hist(nano_nu2_phi*np.ones(int(height)), bins=50, range=(-4, 4), label="gen nu2")
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.title(f"Event {event_id} comb {comb_tag} nu from onshell W phi")
            plt.xlabel("phi")
            plt.savefig(os.path.join(plotting_dir, evt_plt_dir, f"{comb_tag}_nu_onshell_phi.pdf"), format='pdf', bbox_inches='tight')
            plt.close()
            
        
        for i, mass in enumerate(comb_mass):
            plt.hist(mass, weights=comb_weight[i], bins=20, range=(200, 2000), histtype='step', label=comb_label[i])
        estimated_mass = ComputePeakValue(cumulative_mass, weights=cumulative_weights)
        plt.hist(cumulative_mass, weights=cumulative_weights, bins=20, range=(200, 2000), histtype='step', label=f"cumulative, m={estimated_mass:.2f}", linestyle='--')
        plt.grid(True)
        plt.legend(loc='upper right')
        plt.title(f"Event {event_id} likelihoods")
        plt.xlabel("mX, [GeV]")
        plt.savefig(os.path.join(plotting_dir, evt_plt_dir, "mass.pdf"), format='pdf', bbox_inches='tight')
        plt.close()
        break


if __name__ == '__main__':
    main()