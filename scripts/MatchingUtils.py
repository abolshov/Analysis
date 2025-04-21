import vector
import awkward as ak
import numpy as np

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

def JetNIsMatch(branches, gen_branch_name, n):
    m1, m2 = MatchedJetIdx(branches, gen_branch_name)
    return np.logical_or(m1 == n, m2 == n)