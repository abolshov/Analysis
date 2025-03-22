import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
from CutTracker import CutTracker
import argparse
import os


def MissingTrueJetFraction(branches):
    true_jet_tag = branches['centralJet_TrueBjetTag']

    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})

    resolved = b1_p4.deltaR(b2_p4) >= 0.4
    pt_cut = np.logical_and(b1_p4.pt > 20, b2_p4.pt > 20)
    eta_cut = np.logical_and(np.abs(b1_p4.eta) < 2.5, np.abs(b2_p4.eta) < 2.5)
    acceptance = np.logical_and(pt_cut, eta_cut)
    recoverable = np.logical_and(acceptance, resolved)
    true_jet_tag = true_jet_tag[recoverable]
    
    max_jet = ak.max(ak.count(true_jet_tag, axis=1))
    true_jet_tag = ak.fill_none(ak.pad_none(true_jet_tag, max_jet), False)
    n_events = len(true_jet_tag)
    
    return ak.count_nonzero(ak.all(~true_jet_tag, axis=1))/n_events


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


def FirstNContainTrue(branches, n):
    true_jet_tag = branches['centralJet_TrueBjetTag']

    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})

    resolved = b1_p4.deltaR(b2_p4) >= 0.4
    pt_cut = np.logical_and(b1_p4.pt > 20, b2_p4.pt > 20)
    eta_cut = np.logical_and(np.abs(b1_p4.eta) < 2.5, np.abs(b2_p4.eta) < 2.5)
    acceptance = np.logical_and(pt_cut, eta_cut)
    recoverable = np.logical_and(acceptance, resolved)
    true_jet_tag = true_jet_tag[recoverable]
    
    max_jet = ak.max(ak.count(true_jet_tag, axis=1))
    true_jet_tag = ak.fill_none(ak.pad_none(true_jet_tag, max_jet), False)
    n_events = len(true_jet_tag)
    
    first_n = true_jet_tag[:, 0]
    for i in range(n):
        first_n = np.logical_or(first_n, true_jet_tag[:, i])
        
    return ak.count_nonzero(first_n)/n_events


def MakePlot(x, y, plot_params):
    plt.plot(x, y, plot_params['marker'], color=plot_params['color']) 
    plt.title(plot_params['title'])
    plt.xlabel(plot_params['xlabel'])
    plt.ylabel(plot_params['ylabel'])
    plt.xscale('log')
    plt.grid(True)
    plt.savefig(plot_params['name'], format='png', bbox_inches='tight')
    plt.close()


colors = {0: "blue",
          1: "red",
          5: "magenta",
          2: "orange",
          4: "cyan",
          3: "green"}


def main():
    file_list = []
    with open("dl_files.txt", 'r') as input_files:
        file_list = [name for name in input_files if "Radion" in name]    
    all_efficiencies = []
    masspoints = []
    num_jets = [2, 3, 4]

    plotting_dir = "cutflow_radion/"
    os.makedirs(plotting_dir, exist_ok=True)

    for file_name in file_list:
        print(f"processing {file_name}")
        tracker = CutTracker(file_name)
        tracker.Apply(lambda branches: branches['genb1_pt'] > 20, "b1 quark pt")
        tracker.Apply(lambda branches: np.abs(branches['genb1_eta']) < 2.5, "b1 quark eta")
        tracker.Apply(lambda branches: branches['genb2_pt'] > 20, "b2 quark pt")
        tracker.Apply(lambda branches: np.abs(branches['genb2_eta']) < 2.5, "b2 quark eta")
        tracker.Apply(lambda branches: branches['nJet'] >= 2, "jet multiplicity")
        tracker.Apply(lambda branches: branches['lep1_pt'] > 0, "has reco lep1")
        tracker.Apply(lambda branches: branches['lep2_pt'] > 0, "has reco lep2")
        tracker.Apply(lambda branches: DeltaRBetweenBQuarks(branches) > 0.4, "quarks dR")
        tracker.Apply(lambda branches: HasMatching(branches), "matching")
        tracker.Apply(lambda branches: HasOneToOneMatching(branches), "1 to 1 matching")
        tracker.PlotCutflow(percentages=True, save_path=plotting_dir)

        denom = tracker.NumEventsSelected()
        max_match_idx = MaxMatchIdx(tracker.branches)
        efficiencies_for_sample = []
        for n in num_jets:    
            numer = ak.count_nonzero(max_match_idx < n)
            eff = numer/denom
            efficiencies_for_sample.append(eff)
        all_efficiencies.append(efficiencies_for_sample)
        masspoints.append(tracker.masspoint)

    all_efficiencies = np.array(all_efficiencies)
    for n in num_jets:
            plt.plot(masspoints, all_efficiencies[:, n - 2], 's', color=colors[n], label=f"{n}")

    plt.title("B jet selection efficiency for different number of jets")
    plt.xlabel("mX, [GeV]")
    plt.ylabel("Efficiency")
    # plt.legend(loc='upper right')
    plt.legend(loc='best')
    plt.xscale('log')
    plt.grid(True)
    plt.savefig("eff_log.pdf", format='pdf', bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()