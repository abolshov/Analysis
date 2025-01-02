import uproot 
import awkward as ak 
import vector
import numpy as np


HIGGS_MASS = 125.0


def ComputeMass(bj1_p4, bj2_p4, lj1_p4, lj2_p4, lep_p4, met_p4):
    vis_p4 = bj1_p4 + bj2_p4 + lj1_p4 + lj2_p4

    a = HIGGS_MASS*HIGGS_MASS - vis.mass()*vis.mass() + 2.0*vis.px()*met_p4.px() + 2.0*vis.px()*met_p4.py()
    A = 4.0*(vis.E()*vis.E() - vis.pz()*vis.pz())
    B = -4.0*a*vis.pz()
    C = -4.0*vis.E()*vis.E()*(met_p4.px()*met_p4.px() + met_p4.px()*met_p4.px()) - a*a
    delta = B*B - 4.0*A*C

    double pz1 = (-B + np.sqrt(delta))/(2.0*A)
    double pz2 = (-B - np.sqrt(delta))/(2.0*A)

    E1 = np.sqrt(met_p4.px**2 + met_p4.py**2 + pz1**2)
    E2 = np.sqrt(met_p4.px**2 + met_p4.py**2 + pz2**2)

    nu1 = vector.zip({"px": met_p4.px,
                      "py": met_p4.py,
                      "pz": pz1,
                      "E": E1})

    nu2 = vector.zip({"px": met_p4.px,
                      "py": met_p4.py,
                      "pz": pz2,
                      "E": E2})

    sol1 = vis + mu1
    sol2 = vis + nu2

    m1 = ak.nan_to_num(sol1.mass)
    m2 = ak.nan_to_num(sol2.mass)

    return m1, m2


def main():
    print("analytical solution")

    file = uproot.open("nano_0.root")
    tree = file["Events"]
    branches = tree.arrays()

    jet_pt = branches['centralJet_pt']
    jet_eta = branches['centralJet_eta']
    jet_phi =branches['centralJet_phi']
    jet_mass = branches['centralJet_mass']
    jet_p4 = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})

    lep_p4 = vector.zip({'pt':branches['lep1_pt'],
                         'eta':branches['lep1_eta'],
                         'phi':branches['lep1_phi'],
                         'mass':branches['lep1_mass']})

    met_p4 = vector.zip({'pt':branches['PuppiMET_pt'],
                         'eta':0.0,
                         'phi':branches['PuppiMET_phi'],
                         'mass':0.0})

    


if __name__ == "__main__":
    main()