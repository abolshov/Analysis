import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
from CutTracker import CutTracker
import os
import itertools


def MaskIndices(indices, array_to_slice):
    """
    https://stackoverflow.com/questions/40894594/extract-elements-from-numpy-array-that-are-not-in-list-of-indexes
    """
    whole_set, in_set = ak.unzip(ak.cartesian([ak.local_index(array_to_slice), indices], nested=True))
    return ak.any(whole_set == in_set, axis=-1)


def MinDeltaRs(branches, gen_branch_name):
    parton1_p4 = vector.zip({'pt': branches[f'{gen_branch_name}1_pt'],
                             'eta': branches[f'{gen_branch_name}1_eta'],
                             'phi': branches[f'{gen_branch_name}1_phi'],
                             'mass': branches[f'{gen_branch_name}1_mass']})

    parton2_p4 = vector.zip({'pt': branches[f'{gen_branch_name}2_pt'],
                             'eta': branches[f'{gen_branch_name}2_eta'],
                             'phi': branches[f'{gen_branch_name}2_phi'],
                             'mass': branches[f'{gen_branch_name}2_mass']})

    pt = branches['centralJet_pt']
    eta = branches['centralJet_eta']
    phi = branches['centralJet_phi']
    mass = branches['centralJet_mass']
    jets = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})

    dr1 = ak.min(parton1_p4.deltaR(jets), axis=1)
    dr2 = ak.min(parton2_p4.deltaR(jets), axis=1)
    return dr1, dr2

def HasMatching(branches, gen_branch_name):
    dr1, dr2 = MinDeltaRs(branches, gen_branch_name)
    return np.logical_and(dr1 < 0.4, dr2 < 0.4)

def HasOneToOneMatching(branches, gen_branch_name):
    match1_idx, match_idx2 = MatchedJetIdx(branches, gen_branch_name)
    return match1_idx != match_idx2

def MaxMatchIdx(branches, gen_branch_name):
    return np.maximum(*MatchedJetIdx(branches, gen_branch_name))

def DeltaRBetweenQuarks(branches, gen_branch_name):
    parton1_p4 = vector.zip({'pt': branches[f'{gen_branch_name}1_pt'],
                             'eta': branches[f'{gen_branch_name}1_eta'],
                             'phi': branches[f'{gen_branch_name}1_phi'],
                             'mass': branches[f'{gen_branch_name}1_mass']})

    parton2_p4 = vector.zip({'pt': branches[f'{gen_branch_name}2_pt'],
                             'eta': branches[f'{gen_branch_name}2_eta'],
                             'phi': branches[f'{gen_branch_name}2_phi'],
                             'mass': branches[f'{gen_branch_name}2_mass']})

    return parton1_p4.deltaR(parton2_p4)

def MatchedJetIdx(branches, gen_branch_name):
    parton1_p4 = vector.zip({'pt': branches[f'{gen_branch_name}1_pt'],
                             'eta': branches[f'{gen_branch_name}1_eta'],
                             'phi': branches[f'{gen_branch_name}1_phi'],
                             'mass': branches[f'{gen_branch_name}1_mass']})

    parton2_p4 = vector.zip({'pt': branches[f'{gen_branch_name}2_pt'],
                             'eta': branches[f'{gen_branch_name}2_eta'],
                             'phi': branches[f'{gen_branch_name}2_phi'],
                             'mass': branches[f'{gen_branch_name}2_mass']})

    pt = branches['centralJet_pt']
    eta = branches['centralJet_eta']
    phi = branches['centralJet_phi']
    mass = branches['centralJet_mass']
    jets = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})

    match1_idx = ak.argmin(parton1_p4.deltaR(jets), axis=1)
    match2_idx = ak.argmin(parton2_p4.deltaR(jets), axis=1)
    return match1_idx, match2_idx

def main():
    # file_name = "/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root"
    file_name = "/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-800/nano_0.root"
    channel = "SL"
    plotting_dir = f"jet_btag_study/{channel}"
    os.makedirs(plotting_dir, exist_ok=True)

    tracker = CutTracker(file_name)
    tracker.Apply(lambda branches: branches['genb1_pt'] > 20, "b1 quark pt")
    tracker.Apply(lambda branches: np.abs(branches['genb1_eta']) < 2.5, "b1 quark eta")
    tracker.Apply(lambda branches: branches['genb2_pt'] > 20, "b2 quark pt")
    tracker.Apply(lambda branches: np.abs(branches['genb2_eta']) < 2.5, "b2 quark eta")
    if channel == "DL":
        tracker.Apply(lambda branches: branches['ncentralJet'] >= 3, "jet multiplicity")
    elif channel == "SL":
        tracker.Apply(lambda branches: branches['ncentralJet'] >= 4, "jet multiplicity")
    tracker.Apply(lambda branches: branches['lep1_pt'] > 0, "has reco lep1")
    if channel == "DL":
        tracker.Apply(lambda branches: branches['lep2_pt'] > 0, "has reco lep2")
    tracker.Apply(lambda branches: DeltaRBetweenQuarks(branches, "genb") > 0.4, "b quarks dR")
    tracker.Apply(lambda branches: HasMatching(branches, "genb"), "b matching")
    tracker.Apply(lambda branches: HasOneToOneMatching(branches, "genb"), "1 to 1 b matching")
    # tracker.Apply(lambda branches: MaxMatchIdx(branches, "genb") == 2, "jet == 3 is match")
    # tracker.Apply(lambda branches: MaxMatchIdx(branches, "genb") >= 2, "jet >= 3 is match")
    # tracker.Apply(lambda branches: MaxMatchIdx(branches, "genb") < 2, "jets < #3 are matches")
    
    tagger_names = ["PNet", "DeepFlav"]
    tagger_variables = {"PNet": ["B", "QvG", "CvL", "CvB"],
                        "DeepFlav": ["B", "QG", "CvL", "CvB"]}

    if channel == "DL":
        print("this is DL channel")

        tracker.PlotCutflow()
        selected_branches = tracker.GetSelectedEvents()

        m1, m2 = MatchedJetIdx(selected_branches, "genb")

        for tagger_name in tagger_names:
            tagger_variable_names = tagger_variables[tagger_name]
            for tagger_variable_name in tagger_variable_names:
                score_array = selected_branches[f"centralJet_btag{tagger_name}{tagger_variable_name}"]
                for num, match_jet_idx in enumerate([m1, m2]):    
                    plt.hist(ak.flatten(score_array[ak.unflatten(match_jet_idx, 1)]), bins=20, range=(0, 1), histtype='step')
                    plt.title(f"B jet #{num + 1} {tagger_name}{tagger_variable_name} score")
                    plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                    plt.ylabel("Count")
                    plt.grid(True)
                    plt.savefig(os.path.join(plotting_dir, f"btag{tagger_name}{tagger_variable_name}_b{num + 1}.pdf"), format='pdf', bbox_inches='tight')
                    plt.close()

                    plt.hist(score_array[:, 2], bins=20, range=(0, 1), histtype='step', label="jet 3")
                    plt.hist(score_array[:, 1], bins=20, range=(0, 1), histtype='step', label="jet 2")
                    plt.hist(score_array[:, 0], bins=20, range=(0, 1), histtype='step', label="jet 1")
                    plt.title(f"{tagger_name}{tagger_variable_name} scores for first 3 jets")
                    plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                    plt.ylabel("Count")
                    plt.legend(loc='best')
                    plt.grid(True)
                    plt.savefig(os.path.join(plotting_dir, f"btag{tagger_name}{tagger_variable_name}_first_3.pdf"), format='pdf', bbox_inches='tight')
                    plt.close()

                for jet1_idx, jet2_idx in itertools.combinations([0, 1, 2], 2):
                    plt.hist2d(score_array[:, jet1_idx].to_numpy(), score_array[:, jet2_idx].to_numpy(), bins=(20, 20), cmap=plt.cm.viridis, range=[[0, 1], [0, 1]])
                    plt.title(f"btag{tagger_name}{tagger_variable_name} scores")
                    plt.xlabel(f"jet {jet1_idx + 1} btag{tagger_name}{tagger_variable_name} score")
                    plt.ylabel(f"jet {jet2_idx + 1} btag{tagger_name}{tagger_variable_name} score")
                    plt.colorbar()
                    plt.savefig(os.path.join(plotting_dir, f"btag{tagger_name}{tagger_variable_name}_{jet1_idx + 1}vs{jet2_idx + 1}.pdf"), format='pdf', bbox_inches='tight')
                    plt.close()
        
        # num_evts_medium_wp = ak.count_nonzero(btag_score[:, 2] > 0.5, axis=0)
        # frac_med_wp = num_evts_medium_wp/len(btag_score[:, 2])*100
        # print(f"Scenario: {tracker.GetLastCutDescription()}")
        # print(f"Percentage of events with jet #3 passing medium btag WP: {frac_med_wp:.2f}%")
    elif channel == "SL":
        print("this is SL channel")

        exclude_bjets = True

        tracker.Apply(lambda branches: branches['genV2prod1_pt'] > 20, "q1 pt")
        tracker.Apply(lambda branches: np.abs(branches['genV2prod1_eta']) < 5.0, "q1 eta")
        tracker.Apply(lambda branches: branches['genV2prod2_pt'] > 20, "q2 pt")
        tracker.Apply(lambda branches: np.abs(branches['genV2prod2_eta']) < 5.0, "q2 eta")
        tracker.Apply(lambda branches: DeltaRBetweenQuarks(branches, "genV2prod") > 0.4, "light quarks dR")
        tracker.Apply(lambda branches: HasMatching(branches, "genV2prod"), "light matching")
        tracker.Apply(lambda branches: HasOneToOneMatching(branches, "genV2prod"), "1 to 1 light matching")
        tracker.PlotCutflow(percentages=True)
        
        selected_branches = tracker.GetSelectedEvents()
        l1, l2 = MatchedJetIdx(selected_branches, "genV2prod")
        matched_jet_indices = ak.concatenate([ak.unflatten(l1, 1), ak.unflatten(l2, 1)], axis=1)
        if exclude_bjets:
            b1, b2 = MatchedJetIdx(selected_branches, "genb")
            matched_jet_indices = ak.concatenate([matched_jet_indices, 
                                                  ak.unflatten(b1, 1), 
                                                  ak.unflatten(b2, 1)], axis=1)
        matched_jet_indices_mask = MaskIndices(matched_jet_indices, selected_branches["centralJet_pt"])

        for tagger_name in tagger_names:
            tagger_variable_names = tagger_variables[tagger_name]
            for tagger_variable_name in tagger_variable_names:
                score_array = selected_branches[f"centralJet_btag{tagger_name}{tagger_variable_name}"]
                for num, match_jet_idx in enumerate([l1, l2]):    
                    plt.hist(ak.flatten(score_array[ak.unflatten(match_jet_idx, 1)]), bins=20, range=(0, 1), histtype='step')
                    plt.title(f"light jet #{num + 1} {tagger_name}{tagger_variable_name} score")
                    plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                    plt.ylabel("Count")
                    plt.grid(True)
                    plt.savefig(os.path.join(plotting_dir, f"btag{tagger_name}{tagger_variable_name}_q{num + 1}.pdf"), format='pdf', bbox_inches='tight')
                    plt.close()

                max_other_jets_score = ak.max(score_array[~matched_jet_indices_mask], axis=1, mask_identity=False)
                min_other_jets_score = ak.min(score_array[~matched_jet_indices_mask], axis=1, mask_identity=False)

                first_jet_score = ak.flatten(score_array[ak.unflatten(l1, 1)])
                second_jet_score = ak.flatten(score_array[ak.unflatten(l2, 1)])
                max_score = np.maximum(first_jet_score, second_jet_score)
                min_score = np.minimum(first_jet_score, second_jet_score)

                plt.hist(max_score, bins=20, range=(0, 1), histtype='step')
                plt.title(f"Max light jet {tagger_name} {tagger_variable_name} score")
                plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                plt.ylabel("Count")
                plt.grid(True)
                plt.savefig(os.path.join(plotting_dir, f"max_light_btag{tagger_name}{tagger_variable_name}.pdf"), format='pdf', bbox_inches='tight')
                plt.close()

                plt.hist(min_score, bins=20, range=(0, 1), histtype='step')
                plt.title(f"Min light jet {tagger_name} {tagger_variable_name} score")
                plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                plt.ylabel("Count")
                plt.grid(True)
                plt.savefig(os.path.join(plotting_dir, f"min_light_btag{tagger_name}{tagger_variable_name}.pdf"), format='pdf', bbox_inches='tight')
                plt.close()

                plt.hist(max_score, bins=20, range=(0, 1), histtype='step', label="max light")
                plt.hist(max_other_jets_score, bins=20, range=(0, 1), histtype='step', label="max other")
                plt.title(f"Max {tagger_name}{tagger_variable_name} score comparison")
                plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                plt.ylabel("Count")
                plt.legend(loc="best")
                plt.grid(True)
                plt.savefig(os.path.join(plotting_dir, f"max_btag{tagger_name}{tagger_variable_name}_light_vs_other.pdf"), format='pdf', bbox_inches='tight')
                plt.close()

                plt.hist(min_score, bins=20, range=(0, 1), histtype='step', label="min light")
                plt.hist(min_other_jets_score, bins=20, range=(0, 1), histtype='step', label="min other")
                plt.title(f"Min {tagger_name}{tagger_variable_name} score comparison")
                plt.xlabel(f"btag{tagger_name}{tagger_variable_name} score")
                plt.ylabel("Count")
                plt.legend(loc="best")
                plt.grid(True)
                plt.savefig(os.path.join(plotting_dir, f"min_btag{tagger_name}{tagger_variable_name}_light_vs_other.pdf"), format='pdf', bbox_inches='tight')
                plt.close()

                plt.hist2d(max_score.to_numpy(), max_other_jets_score.to_numpy(), bins=(20, 20), cmap=plt.cm.viridis, range=[[0, 1], [0, 1]])
                plt.title(f"Max {tagger_name}{tagger_variable_name} score")
                plt.xlabel(f"Max light jet {tagger_name}{tagger_variable_name} score")
                plt.ylabel(f"Max other jet {tagger_name}{tagger_variable_name} score")
                plt.colorbar()
                plt.savefig(os.path.join(plotting_dir, f"max_btag{tagger_name}{tagger_variable_name}_light_vs_other_2D.pdf"), format='pdf', bbox_inches='tight')
                plt.close()

                plt.hist2d(min_score.to_numpy(), min_other_jets_score.to_numpy(), bins=(20, 20), cmap=plt.cm.viridis, range=[[0, 1], [0, 1]])
                plt.title(f"Min {tagger_name}{tagger_variable_name} score")
                plt.xlabel(f"Min light jet {tagger_name}{tagger_variable_name} score")
                plt.ylabel(f"Min other jet {tagger_name}{tagger_variable_name} score")
                plt.colorbar()
                plt.savefig(os.path.join(plotting_dir, f"min_btag{tagger_name}{tagger_variable_name}_light_vs_other_2D.pdf"), format='pdf', bbox_inches='tight')
                plt.close()


if __name__ == '__main__':
    main()