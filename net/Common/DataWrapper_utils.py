import awkward as ak
import vector
import numpy as np
import pandas as pd

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def IsCorrectLepton(df, lep):
    if lep not in [1, 2]:
        raise RuntimeError("Can have at most 2 leptons")

    reco_lep_type = df[f'lep{lep}_type']
    gen_lep_type = df[f'lep{lep}_genLep_kind']

    reco_lep_mu = reco_lep_type == 2
    reco_lep_ele = reco_lep_type == 1

    gen_lep_mu = np.logical_or(gen_lep_type == 2, gen_lep_type == 4)
    gen_lep_ele = np.logical_or(gen_lep_type == 1, gen_lep_type == 3)
    
    lep_correct = np.logical_or(np.logical_and(reco_lep_mu, gen_lep_mu), np.logical_and(reco_lep_ele, gen_lep_ele))
    
    return lep_correct


def TransformPNetFactorsToResolutions(df, n_jets):
    for i in range(n_jets):
        name = f"jet{i + 1}"
        pt = np.sqrt(df[f"{name}_px"]*df[f"{name}_px"] + df[f"{name}_py"]*df[f"{name}_py"])
        pt *= df[f"{name}_PNetRegPtRawCorr"]
        df[f"{name}_PNetRegPtRawRes"] = df[f"{name}_PNetRegPtRawRes"]*pt

    return df


# add px, py, pz, E of n objects with base name to df from array of their 4-vectors
# pandas is issuing performance warnings about adding columns in a loop
# may switch to creating a tmp df and using pd.concat
def AddKinematics(df, p4, base_name, n):
    objects = [f"{base_name}{i + 1}" for i in range(n)]
    variables = ['px', 'py', 'pz', 'E']
    func_map = {'px': lambda x: x.px,
                'py': lambda x: x.py,
                'pz': lambda x: x.pz,
                'E': lambda x: x.E}

    tmp_dict = {}
    for obj in objects:
        idx = int(obj[-1]) - 1
        obj_p4 = p4[:, idx]
        for var in variables:
            name = f'{obj}_{var}'
            tmp_dict[name] = func_map[var](obj_p4).to_numpy()

    df = pd.concat([df, pd.DataFrame.from_dict(tmp_dict)], axis=1)
    return df


# create df with px, py, pz, E of jets, leptons and met
def AddKinematicFeatures(df, branches, n_lep, n_jet):
    # make an array of jet p4 with resulting shape (n_events, n_jet)
    pt = ak.fill_none(ak.pad_none(branches['centralJet_pt'], n_jet), 0.0)
    eta = ak.fill_none(ak.pad_none(branches['centralJet_eta'], n_jet), 0.0)
    phi = ak.fill_none(ak.pad_none(branches['centralJet_phi'], n_jet), 0.0)
    mass = ak.fill_none(ak.pad_none(branches['centralJet_mass'], n_jet), 0.0)
    jet = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})
    jet = jet[:, :n_jet]
    df = AddKinematics(df, jet, "jet", n_jet)
    
    lep1_p4 = vector.zip({'pt': branches['lep1_pt'],
                          'eta': branches['lep1_eta'],
                          'phi': branches['lep1_phi'],
                          'mass': branches['lep1_mass']})
    sz = len(lep1_p4)
    zero_arr = ak.Array(np.zeros(sz))
    lep2_p4 = vector.zip({'pt': zero_arr, 'eta': zero_arr, 'phi': zero_arr, 'mass': zero_arr})
    if n_lep == 2:
        lep2_p4 = vector.zip({'pt': branches['lep2_pt'],
                              'eta': branches['lep2_eta'],
                              'phi': branches['lep2_phi'],
                              'mass': branches['lep2_mass']})
    lep = ak.concatenate([ak.unflatten(lep1_p4, 1), ak.unflatten(lep2_p4, 1)], axis=1)
    df = AddKinematics(df, lep, "lep", n_lep)

    met_p4 = lep1_p4 = vector.zip({'pt': branches['PuppiMET_pt'], 'eta': 0.0, 'phi': branches['PuppiMET_phi'], 'mass': 0.0})
    df['met_px'] = met_p4.px.to_numpy()
    df['met_py'] = met_p4.py.to_numpy()

    return df


def IsKinematic(feature_name):
    for var in ['px', 'py', 'pz', 'E']:
        if var in feature_name:
            return True
    return False


def AddJetFeatures(df, branches, feature_list, n_jet):
    for feature in feature_list:
        if IsKinematic(feature):
            continue
        
        branch_name = f"centralJet_{feature}"
        feature_branch = ak.fill_none(ak.pad_none(branches[branch_name], n_jet), 0.0)

        tmp_dict = {f"jet{i + 1}_{feature}": feature_branch[:, i].to_numpy() for i in range(n_jet)}
        df = pd.concat([df, pd.DataFrame.from_dict(tmp_dict)], axis=1)
        
    return df


def AddMindR(df, branches, n_lep, n_jet):
    pt = ak.fill_none(ak.pad_none(branches['centralJet_pt'], n_jet), 0.0)
    eta = ak.fill_none(ak.pad_none(branches['centralJet_eta'], n_jet), 0.0)
    phi = ak.fill_none(ak.pad_none(branches['centralJet_phi'], n_jet), 0.0)
    mass = ak.fill_none(ak.pad_none(branches['centralJet_mass'], n_jet), 0.0)
    jets_p4 = vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})
    jets_p4 = jets_p4[:, :n_jet]

    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})

    mindR_b1 = ak.min(b1_p4.deltaR(jets_p4), axis=1)
    mindR_b2 = ak.min(b2_p4.deltaR(jets_p4), axis=1)
    df['mindR_b1'] = mindR_b1.to_numpy()
    df['mindR_b2'] = mindR_b2.to_numpy()

    if n_lep == 1:
        q1_p4 = vector.zip({'pt': branches['genV2prod1_pt'],
                            'eta': branches['genV2prod1_eta'],
                            'phi': branches['genV2prod1_phi'],
                            'mass': branches['genV2prod1_mass']})

        q2_p4 = vector.zip({'pt': branches['genV2prod2_pt'],
                            'eta': branches['genV2prod2_eta'],
                            'phi': branches['genV2prod2_phi'],
                            'mass': branches['genV2prod2_mass']})

        mindR_q1 = ak.min(q1_p4.deltaR(jets_p4), axis=1)
        mindR_q2 = ak.min(q2_p4.deltaR(jets_p4), axis=1)
        df['mindR_q1'] = mindR_q1.to_numpy()
        df['mindR_q2'] = mindR_q2.to_numpy()

    return df    


def ApplyFiducialSelection(df, branches, n_lep, n_jet, apply_gen_reco_match):
    b1_p4 = vector.zip({'pt': branches['genb1_pt'],
                       'eta': branches['genb1_eta'],
                       'phi': branches['genb1_phi'],
                       'mass': branches['genb1_mass']})

    b2_p4 = vector.zip({'pt': branches['genb2_pt'],
                        'eta': branches['genb2_eta'],
                        'phi': branches['genb2_phi'],
                        'mass': branches['genb2_mass']})

    col_names = ['b1_pt', 'b2_pt', 'b1_eta', 'b2_eta', 'b1_vis_pt', 'b2_vis_pt',
                 'n_jet', 'lep1_type', 'lep2_type', 'lep1_genLep_kind', 'lep2_genLep_kind']
    branch_name_map = {'b1_pt': 'genb1_pt',
                       'b1_eta': 'genb1_eta',
                       'b2_pt': 'genb2_pt',
                       'b2_eta': 'genb2_eta',
                       'b1_vis_pt': 'genb1_vis_pt',
                       'b2_vis_pt': 'genb2_vis_pt',
                       'n_jet': 'ncentralJet',
                       'lep1_type': 'lep1_type', 
                       'lep2_type': 'lep2_type', 
                       'lep1_genLep_kind': 'lep1_genLep_kind', 
                       'lep2_genLep_kind': 'lep2_genLep_kind'}
    tmp_dict = {name: branches[branch_name_map[name]].to_numpy() for name in col_names}
    tmp_dict['bb_dr'] = b1_p4.deltaR(b2_p4)
    df = pd.concat([df, pd.DataFrame.from_dict(tmp_dict)], axis=1)

    if n_lep == 1:
        q1_p4 = vector.zip({'pt': branches['genV2prod1_pt'],
                            'eta': branches['genV2prod1_eta'],
                            'phi': branches['genV2prod1_phi'],
                            'mass': branches['genV2prod1_mass']})

        q2_p4 = vector.zip({'pt': branches['genV2prod2_pt'],
                            'eta': branches['genV2prod2_eta'],
                            'phi': branches['genV2prod2_phi'],
                            'mass': branches['genV2prod2_mass']})

        col_names = ['q1_pt', 'q2_pt', 'q1_eta', 'q2_eta', 'q1_vis_pt', 'q2_vis_pt']
        branch_name_map = {'q1_pt': 'genV2prod1_pt',
                           'q1_eta': 'genV2prod1_eta',
                           'q2_pt': 'genV2prod2_pt',
                           'q2_eta': 'genV2prod2_eta',
                           'q1_vis_pt': 'genV2prod1_vis_pt',
                           'q2_vis_pt': 'genV2prod2_vis_pt'}

        tmp_dict = {name: branches[branch_name_map[name]].to_numpy() for name in col_names}
        tmp_dict['qq_dr'] = q1_p4.deltaR(q2_p4)
        df = pd.concat([df, pd.DataFrame.from_dict(tmp_dict)], axis=1)

    if apply_gen_reco_match:
        df = AddMindR(df, branches, n_lep, n_jet)

    # apply cuts
    # cuts to b quarks are applied anyway
    df = df[df['b1_vis_pt'] > 0.0]
    df = df[df['b2_vis_pt'] > 0.0]
    df = df[df['bb_dr'] > 0.4]
    df = df[df['b1_pt'] > 20.0]
    df = df[df['b2_pt'] > 20.0]
    df = df[np.abs(df['b1_eta']) < 2.5]
    df = df[np.abs(df['b2_eta']) < 2.5]
    if apply_gen_reco_match:
        # I want to apply this cut in the very end to make sure I already have correct objects
        df = df[df['mindR_b1'] < 0.4]
        df = df[df['mindR_b2'] < 0.4]
        # cut = np.logical_or(df['mindR_b1'] >= 0.4, df['mindR_b2'] >= 0.4)
        # df = df[cut]

    if n_lep == 2:
        # cuts specific to DL channel
        df = df[df['n_jet'] >= 2]
        df = df[IsCorrectLepton(df, 1)]
        df = df[IsCorrectLepton(df, 2)]
    elif n_lep == 1:
        # cuts specific to SL channel
        df = df[df['n_jet'] >= 4]
        df = df[df['q1_vis_pt'] > 0.0]
        df = df[df['q2_vis_pt'] > 0.0]
        df = df[df['qq_dr'] > 0.4]
        df = df[df['q1_pt'] > 20.0]
        df = df[df['q2_pt'] > 20.0]
        df = df[np.abs(df['q1_eta']) < 5.0]
        df = df[np.abs(df['q2_eta']) < 5.0]
        df = df[IsCorrectLepton(df, 1)]
        if apply_gen_reco_match:
            df = df[df['mindR_q1'] < 0.4]
            df = df[df['mindR_q2'] < 0.4]
    else:
        raise RuntimeError("Wrong number of leptons")

    return df