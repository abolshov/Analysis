import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
from CutTracker import CutTracker
from MatchingUtils import *
import argparse
import os


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

    all_jet0_true = []
    all_jet1_true = []
    all_jet2_true = []
    all_jet_greater2_true = []

    channel = "DL"
    plotting_dir = f"bjet_select_eff_study/{channel}"
    os.makedirs(plotting_dir, exist_ok=True)
    os.makedirs(os.path.join(plotting_dir, "cutflow_radion"), exist_ok=True)

    for file_name in file_list:
        print(f"processing {file_name}")

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
        tracker.PlotCutflow(percentages=True, save_path=os.path.join(plotting_dir, "cutflow_radion"))

        selected_branches = tracker.GetSelectedEvents()
        m1, m2 = MatchedJetIdx(selected_branches, "genb")
        max_match_idx = np.maximum(m1, m2)

        denom = tracker.NumEventsSelected()
        efficiencies_for_sample = []
        for n in num_jets:    
            numer = ak.count_nonzero(max_match_idx < n)
            eff = numer/denom
            efficiencies_for_sample.append(eff)
        all_efficiencies.append(efficiencies_for_sample)
        masspoints.append(tracker.masspoint)

        have_1to1_matching = tracker.cutflow["1 to 1 b matching"]
        tracker.Apply(lambda branches: JetNIsMatch(branches, "genb", 0), "jet 0 is true")
        jet0_true = tracker.cutflow["jet 0 is true"]

        all_jet0_true.append(jet0_true/have_1to1_matching)

        selected_branches = tracker.GetSelectedEvents()
        m1, m2 = MatchedJetIdx(selected_branches, "genb")
        jet1_true = ak.count_nonzero(np.logical_or(m1 == 1, m2 == 1))
        jet2_true = ak.count_nonzero(np.logical_or(m1 == 2, m2 == 2))
        jet_greater2_true = ak.count_nonzero(np.maximum(m1, m2) > 2)

        all_jet1_true.append(jet1_true/jet0_true)
        all_jet2_true.append(jet2_true/jet0_true)
        all_jet_greater2_true.append(jet_greater2_true/jet0_true)
        

    all_efficiencies = np.array(all_efficiencies)
    for n in num_jets:
        plt.plot(masspoints, all_efficiencies[:, n - 2], 's', color=colors[n], label=f"{n}")

    plt.title("B jet selection efficiency for different number of jets")
    plt.xlabel("mX, [GeV]")
    plt.ylabel("Efficiency")
    plt.legend(loc='best')
    plt.xscale('log')
    plt.grid(True)
    plt.savefig(os.path.join(plotting_dir, "eff_log.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    plt.plot(masspoints, np.array(all_jet0_true), 's', color=colors[1])
    plt.title("Fraction of events where jet with highest btag is true")
    plt.xlabel("mX, [GeV]")
    plt.ylabel("Fraction")
    # plt.xscale('log')
    plt.grid(True)
    plt.savefig(os.path.join(plotting_dir, "highest_btag_true.pdf"), format='pdf', bbox_inches='tight')
    plt.close()

    true_jet_label = {0: "1", 1: "2", 2: ">2"}
    sums = [all_jet1_true[i] + all_jet2_true[i] + all_jet_greater2_true[i] for i in range(len(masspoints))]
    for i, values in enumerate([np.array(all_jet1_true), np.array(all_jet2_true), np.array(all_jet_greater2_true)]):
        plt.plot(masspoints, values, 's', color=colors[i], label=true_jet_label[i])

    plt.plot(masspoints, sums, 's', color=colors[4], label="sum")
    plt.title("Match index distribuion")
    plt.xlabel("mX, [GeV]")
    plt.ylabel("Fraction")
    plt.legend(loc='best')
    plt.xscale('log')
    plt.grid(True)
    plt.savefig(os.path.join(plotting_dir, "match_idx_dsit.pdf"), format='pdf', bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()