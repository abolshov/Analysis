import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
from CutTracker import CutTracker
import os


def MinDeltaRs(branches):
    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})

    pt = branches['centralJet_pt']
    eta = branches['centralJet_eta']
    phi = branches['centralJet_phi']
    mass = branches['centralJet_mass']
    jets = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})

    dr_b1 = ak.min(b1_p4.deltaR(jets), axis=1)
    dr_b2 = ak.min(b2_p4.deltaR(jets), axis=1)
    return dr_b1, dr_b2

def HasMatching(branches):
    dr_b1, dr_b2 = MinDeltaRs(branches)
    return np.logical_and(dr_b1 < 0.4, dr_b2 < 0.4)

def HasOneToOneMatching(branches):
    b1_match_idx, b2_match_idx = MatchedJetIdx(branches)
    return b1_match_idx != b2_match_idx

def MaxMatchIdx(branches):
    return np.maximum(*MatchedJetIdx(branches))

def DeltaRBetweenBQuarks(branches):
    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})
    return b1_p4.deltaR(b2_p4)

def MatchedJetIdx(branches):
    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})

    pt = branches['centralJet_pt']
    eta = branches['centralJet_eta']
    phi = branches['centralJet_phi']
    mass = branches['centralJet_mass']
    jets = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})

    b1_match_idx = ak.argmin(b1_p4.deltaR(jets), axis=1)
    b2_match_idx = ak.argmin(b2_p4.deltaR(jets), axis=1)
    return b1_match_idx, b2_match_idx

def main():
    file_name = "/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root"
    plotting_dir = "jet_btag_study"
    os.makedirs(plotting_dir, exist_ok=True)

    tracker = CutTracker(file_name)
    tracker.Apply(lambda branches: branches['genb1_pt'] > 20, "b1 quark pt")
    tracker.Apply(lambda branches: np.abs(branches['genb1_eta']) < 2.5, "b1 quark eta")
    tracker.Apply(lambda branches: branches['genb2_pt'] > 20, "b2 quark pt")
    tracker.Apply(lambda branches: np.abs(branches['genb2_eta']) < 2.5, "b2 quark eta")
    tracker.Apply(lambda branches: branches['ncentralJet'] >= 3, "jet multiplicity")
    tracker.Apply(lambda branches: branches['lep1_pt'] > 0, "has reco lep1")
    tracker.Apply(lambda branches: branches['lep2_pt'] > 0, "has reco lep2")
    tracker.Apply(lambda branches: DeltaRBetweenBQuarks(branches) > 0.4, "quarks dR")
    tracker.Apply(lambda branches: HasMatching(branches), "matching")
    tracker.Apply(lambda branches: HasOneToOneMatching(branches), "1 to 1 matching")
    tracker.Apply(lambda branches: MaxMatchIdx(branches) == 2, "jet #3 is match")
    # tracker.Apply(lambda branches: MaxMatchIdx(branches) < 2, "jets < #3 are matches")

    selected_branches = tracker.GetSelectedEvents()
    btag_score = selected_branches["centralJet_btagPNetB"]

    tracker.Print()

    # plot btag score of 3d jet when it is a match
    plt.hist(btag_score[:, 2], bins=20, range=(0, 1), histtype='step')
    plt.title("Jet #3 btag score")
    plt.xlabel("btagPNetB score")
    plt.ylabel("Count")
    plt.grid(True)
    plt.savefig(os.path.join(plotting_dir, "btag_3.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    plt.hist(btag_score[:, 2], bins=20, range=(0, 1), histtype='step', label="jet 3")
    plt.hist(btag_score[:, 1], bins=20, range=(0, 1), histtype='step', label="jet 2")
    plt.hist(btag_score[:, 0], bins=20, range=(0, 1), histtype='step', label="jet 1")
    plt.title("Btag scores")
    plt.xlabel("btagPNetB score")
    plt.ylabel("Count")
    plt.legend(loc='best')
    plt.grid(True)
    plt.savefig(os.path.join(plotting_dir, "btag_cmp.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    plt.hist2d(btag_score[:, 2].to_numpy(), btag_score[:, 1].to_numpy(), bins=(10, 10), cmap=plt.cm.viridis, range=[[0, 1], [0, 1]])
    plt.title("Btag scores")
    plt.xlabel("jet 3 btag")
    plt.ylabel("jet 2 btag")
    plt.colorbar()
    plt.savefig(os.path.join(plotting_dir, "btag3btag2.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    plt.hist2d(btag_score[:, 2].to_numpy(), btag_score[:, 0].to_numpy(), bins=(10, 10), cmap=plt.cm.viridis, range=[[0, 1], [0, 1]])
    plt.title("Btag scores")
    plt.xlabel("jet 3 btag")
    plt.ylabel("jet 1 btag")
    plt.colorbar()
    plt.savefig(os.path.join(plotting_dir, "btag3btag1.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    num_evts_medium_wp = ak.count_nonzero(btag_score[:, 2] > 0.5, axis=0)
    frac_med_wp = num_evts_medium_wp/len(btag_score[:, 2])*100
    print(f"Scenario: {tracker.GetLastCutDescription()}")
    print(f"Percentage of events with jet #3 passing medium btag WP: {frac_med_wp:.2f}%")

if __name__ == '__main__':
    main()