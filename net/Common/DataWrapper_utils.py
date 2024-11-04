import awkward as ak
import vector
import numpy as np


def Px(obj_p4):
    return obj_p4.px


def Py(obj_p4):
    return obj_p4.py


def Pz(obj_p4):
    return obj_p4.pz


def E(obj_p4):
    return obj_p4.E


def GetNumPyArray(awk_arr, tot_length, i):
    return ak.to_numpy(ak.fill_none(ak.pad_none(awk_arr[:, :tot_length], tot_length), 0.0))[:, i]


def EventTopology(branches, bb_res, bb_boost, qq_res, qq_boost):
    bb_selection = ak.Array(np.ones(len(branches)) == 1)
    b1_p4 = vector.zip({'pt': branches['genb1_pt'], 
                        'eta': branches['genb1_eta'], 
                        'phi': branches['genb1_phi'], 
                        'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'], 
                        'eta': branches['genb2_eta'], 
                        'phi': branches['genb2_phi'], 
                        'mass': branches['genb2_mass']})

    bb_dr = b1_p4.deltaR(b2_p4)

    if bb_res and not bb_boost:
        bb_selection = bb_dr >= 0.4
    elif bb_boost and not bb_res:
        bb_selection = bb_dr < 0.4

    qq_selection = ak.Array(np.ones(len(branches)) == 1)
    q1_p4 = vector.zip({'pt': branches['genV2prod1_pt'], 
                        'eta': branches['genV2prod1_eta'], 
                        'phi': branches['genV2prod1_phi'], 
                        'mass': branches['genV2prod1_mass']})

    q2_p4 = vector.zip({'pt': branches['genV2prod2_pt'], 
                        'eta': branches['genV2prod2_eta'], 
                        'phi': branches['genV2prod2_phi'], 
                        'mass': branches['genV2prod2_mass']})
    
    qq_dr = q1_p4.deltaR(q2_p4)

    if qq_res and not qq_boost:
        bb_selection = qq_dr >= 0.4
    elif qq_boost and not qq_res:
        bb_selection = qq_dr < 0.4

    event_topology = bb_selection & qq_selection
    return event_topology


def IsCorrectLepton(branches, lep):
    if lep not in [1, 2]:
        raise RuntimeError("Can have at most 2 leptons")

    reco_lep_type = branches[f'lep{lep}_type']
    gen_lep_type = branches[f'lep{lep}_genLep_kind']

    reco_lep_mu = reco_lep_type == 2
    reco_lep_ele = reco_lep_type == 1

    gen_lep_mu = (gen_lep_type == 2) | (gen_lep_type == 4)
    gen_lep_ele = (gen_lep_type == 1) | (gen_lep_type == 3)

    lep_correct = (reco_lep_mu & gen_lep_mu) | (reco_lep_ele & gen_lep_ele)
    return lep_correct


def Acceptance(branches):
    jet_multiplicity = (branches['ncentralJet'] >= 4)

    lep1_correct = IsCorrectLepton(branches, 1)

    light_quarks_accept = ((branches['genV2prod1_pt'] > 20.0) & (np.abs(branches['genV2prod1_eta']) < 5.0)) & ((branches['genV2prod2_pt'] > 20.0) & (np.abs(branches['genV2prod2_eta']) < 5.0))
    b_quarks_accept = ((branches['genb1_pt'] > 20.0) & (np.abs(branches['genb1_eta']) < 2.5)) & ((branches['genb2_pt'] > 20.0) & (np.abs(branches['genb2_eta']) < 2.5))
    quarks_accept = light_quarks_accept & b_quarks_accept

    matching = (branches['genV2prod1_vis_pt'] > 0.0) & (branches['genV2prod2_vis_pt'] > 0.0) & (branches['genb1_vis_pt'] > 0.0) & (branches['genb2_vis_pt'] > 0.0)

    return quarks_accept & matching & jet_multiplicity & lep1_correct