import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
import argparse


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
    true_jet_tag = true_jet_tag[resolved]
    
    max_jet = ak.max(ak.count(true_jet_tag, axis=1))
    true_jet_tag = ak.fill_none(ak.pad_none(true_jet_tag, max_jet), False)
    n_events = len(true_jet_tag)
    
    return ak.count_nonzero(ak.all(~true_jet_tag, axis=1))/n_events


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
    true_jet_tag = true_jet_tag[resolved]
    
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
          2: "green",
          3: "orange",
          4: "cyan",
          5: "magenta"}


def main():
    parser = argparse.ArgumentParser(prog='bjet_select_eff_study')
    parser.add_argument('files', type=str, help="List of files")
    parser.add_argument('n_jets', type=int, help="Number of jets to check")
    parser.add_argument('x_type', type=str, help="Type of X (Graviton/Radion)")
    parser.add_argument('channel', type=str, help="Channel (DL/SL)")

    args = parser.parse_args()

    files = args.files
    n_jets = args.n_jets
    x_type = args.x_type
    ch = args.channel

    if n_jets > 5:
        raise RuntimeError("Maximum number of jets is too large")

    with open(files, 'r') as file_list:
        masspoints = []
        # selection_eff = []
        missing_fraction = []

        selection_efficiencies = [[] for _ in range(n_jets)]

        for name in file_list:
            if x_type not in name:
                continue

            if ch == 'DL' and '2B2JLNu' in name:
                continue

            if ch == 'SL' and '2B2L2Nu' in name:
                continue

            file = uproot.open(name)
            tree = file['Events']
            branches = tree.arrays()

            mp = branches['X_mass'][0]
            masspoints.append(mp)

            # eff = FirstNContainTrue(branches, n_jets)
            # selection_eff.append(eff)

            for n in range(n_jets):
                eff = FirstNContainTrue(branches, n + 1)
                selection_efficiencies[n].append(eff)

            mis_frac = MissingTrueJetFraction(branches)
            missing_fraction.append(mis_frac)

        # MakePlot(masspoints, selection_eff, {'marker': 's', 
        #                                      'title': f'B jet selection efficiency: first {n_jets} jets', 
        #                                      'xlabel': 'mX, [GeV]', 
        #                                      'ylabel': 'Efficiency',
        #                                      'name': f'selection_{x_type}_{ch}.png',
        #                                      'color': 'orange'})
        MakePlot(masspoints, missing_fraction, {'marker': 's', 
                                                'title': 'Fraction of events with no b jets', 
                                                'xlabel': 'mX, [GeV]', 
                                                'ylabel': 'Missing fraction',
                                                'name': f'missing_fraction_{x_type}_{ch}.png',
                                                'color': 'blue'})

        for n in range(n_jets):
            plt.plot(masspoints, selection_efficiencies[n], 's', color=colors[n], label=f"{n + 1}")

        plt.title("B jet selection efficiency for different number of jets")
        plt.xlabel("mX, [GeV]")
        plt.ylabel("Efficiency")
        plt.legend(loc='upper right')
        plt.xscale('log')
        plt.grid(True)
        plt.savefig("selection_cmp.png", format='png', bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main()