import tensorflow as tf
import pandas as pd
import numpy as np
import vector
import awkward as ak
import os
import yaml
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

from Dataloader import Dataloader
from ErrorProp import ErrorPropagator

from MiscUtils import *
from PlotUtils import PlotHist, PlotCompare2D, PlotCovarMtrx, PlotHistStack
from ModelUtils import LoadModel

def PreparePredictions(model, training_params, X_test):
    ys_pred = None
    if training_params['standardize']:
        input_means = training_params['input_train_means']
        input_scales = training_params['input_train_scales']
        X_test -= input_means
        X_test /= input_scales

        print(X_test.shape)

        ys_pred = model.predict(X_test)

        if training_params['add_mass_loss']:
            # drop dimension containing concatenation of all elements before it
            ys_pred = ys_pred[:-1]

        # now ys_pred has shape (n_heads, n_events, n_quantiles) 
        ys_pred = np.array(ys_pred)
        target_scales = np.array(training_params['target_train_scales'])
        target_means = np.array(training_params['target_train_means'])

        # broadcast arrays with means from shape (n_heads,) to shape of ys_pred by inserting new axes
        ys_pred *= target_scales[:, None, None]
        ys_pred += target_means[:, None, None]
        # swap heads and events dimension so events is first
        # (n_heads, n_events, n_quantiles) → (n_events, n_heads, n_quantiles)
        ys_pred = ys_pred.transpose(1, 0, 2)
    else:
        ys_pred = model.predict(X)
        if training_params['add_mass_loss']:
            ys_pred = ys_pred[:-1]
        ys_pred = np.array(ys_pred)
        ys_pred = ys_pred.transpose(1, 0, 2)

    return ys_pred

def ComputeMass(central):
    hvv_en = central[:, 3]
    hbb_en = central[:, -1]

    hvv_p3 = central[:, :3]
    hbb_p3 = central[:, 4:-1]

    x_en = hvv_en + hbb_en
    x_p3 = hvv_p3 + hbb_p3
    x_mass_sqr = np.square(x_en) - np.sum(np.square(x_p3), axis=1)
    neg_mass = x_mass_sqr <= 0.0
    x_mass = np.sqrt(np.abs(x_mass_sqr))
    x_mass = np.where(neg_mass, -1.0, x_mass)
    return x_mass

def main():
    # load signal data
    file = '../train_data/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M_800/nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)

    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)

    # def SignalRegionSelection(df, mod, parity, mass, sample_type):
    #     selection = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
    #     selection = np.logical_and(selection, df['sample_type'] == sample_type)
    #     two_btagged_jets = np.logical_and(df['centralJet_1_btagPNetB'] > 0.245, df['centralJet_2_btagPNetB'] > 0.245)
    #     lep1_p4 = vector.zip({'px': df['lep1_px'], 'py': df['lep1_py'], 'pz': df['lep1_pz'], 'E': df['lep1_E']})
    #     lep2_p4 = vector.zip({'px': df['lep2_px'], 'py': df['lep2_py'], 'pz': df['lep2_pz'], 'E': df['lep2_E']})
    #     same_flavor = df['lep1_legType'] == df['lep2_legType']
    #     mll = (lep1_p4 + lep2_p4).mass.to_numpy()
    #     z_peak_veto = np.logical_and(same_flavor, np.abs(mll - 90.0) < 10.0)
    #     selection = np.logical_and(selection, two_btagged_jets)
    #     # selection = np.logical_and(selection, mll < 12.0)
    #     selection = np.logical_and(selection, z_peak_veto)
    #     return selection

    # X_sig, input_names, y_sig, target_names = dataloader.Get(SignalRegionSelection, 2, 1, 800, 1)  
    X_sig, input_names, y_sig, target_names = dataloader.Get(TestSelection, 2, 1, 800, 1)

    # load model
    model, training_params = LoadModel('model_cfg.yaml')
    print(model.summary())

    ys_pred_sig = PreparePredictions(model, training_params, X_sig)

    # extract prediction for central quantile for all targets
    central_sig = ys_pred_sig[:, :, 1]
    assert y_sig.shape == central_sig.shape

    # hvv_PtEtaPhiE_pred = ToPtEtaPhiE(central_sig[:, :4])
    # hvv_PtEtaPhiE_true = ToPtEtaPhiE(y[:, :4])

    # hbb_PtEtaPhiE_pred = ToPtEtaPhiE(central_sig[:, 4:])
    # hbb_PtEtaPhiE_true = ToPtEtaPhiE(y[:, 4:])

    # pred_PtEtaPhiE = np.concatenate([hvv_PtEtaPhiE_pred, hbb_PtEtaPhiE_pred], axis=-1)
    # true_PtEtaPhiE = np.concatenate([hvv_PtEtaPhiE_true, hbb_PtEtaPhiE_true], axis=-1)

    # var_names = ['genHVV_pt', 'genHVV_eta', 'genHVV_phi', 'genHVV_E', 
    #              'genHbb_pt', 'genHbb_eta', 'genHbb_phi', 'genHbb_E']
    # for i, name in enumerate(var_names):
    #     var = name.split('_')[-1]
    #     obj = name.split('_')[0]

    #     bin_left = np.min([np.quantile(true_PtEtaPhiE[:, i], 0.01), np.quantile(pred_PtEtaPhiE[:, i], 0.01)]) - 1.0
    #     bin_right = np.max([np.quantile(true_PtEtaPhiE[:, i], 0.99), np.quantile(pred_PtEtaPhiE[:, i], 0.99)]) + 1.0
    #     bins = np.linspace(bin_left, bin_right, 100)

    #     PlotCompare2D(true_PtEtaPhiE[:, i], 
    #                   pred_PtEtaPhiE[:, i], 
    #                   var,
    #                   obj,
    #                   bins=bins, 
    #                   title=f'{objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} comparison',
    #                   xlabel=f'True {ground_truth_map[name]}',
    #                   ylabel=f'Predicted {ground_truth_map[name]}',
    #                   plotting_dir=os.path.join(training_params['model_dir'], 'plots'))

    # compute errors as up quantile (-1) minus down quantile (0)
    pred_errors_sig = ys_pred_sig[:, :, -1] - ys_pred_sig[:, :, 0]
    true_errors = y_sig - central_sig

    if training_params['use_energy_layer']:
        pred_errors_sig = np.delete(pred_errors_sig, [3, 7], axis=1)

    global_corr_mtrx = np.array(training_params['global_corr_mtrx'])
    ep = ErrorPropagator(global_corr_mtrx, central_sig[:, :3], central_sig[:, 4:7])
    prop_errors_sig = ep.Propagate(pred_errors_sig)

    pretty_labels = [ground_truth_map[name] for name in target_names if 'E' not in name]
    PlotCovarMtrx(global_corr_mtrx, 'empirical', pretty_labels, os.path.join(training_params['model_dir'], 'plots'))

    bin_left = np.min(prop_errors_sig) - 1.0
    bin_right = np.max(prop_errors_sig) + 1.0
    bins = np.linspace(bin_left, bin_right, 50)
    PlotHist(data=prop_errors_sig, 
             bins=bins,
             title="Predicted X->HH mass errors",
             ylabel='Count',
             xlabel='Mass error, [GeV]',
             plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
             file_name='sig_mx_error')


    # load bkg data
    file = 'TTbar.root'
    # file = 'DY.root'
    bkg_latex = '$T\\bar{T}$'
    bkg_name = 'ttbar'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)

    # def SignalRegionSelection(df):
    #     selection = np.logical_and(df['centralJet_1_btagPNetB'] > 0.245, df['centralJet_2_btagPNetB'] > 0.245)
    #     lep1_p4 = vector.zip({'px': df['lep1_px'], 'py': df['lep1_py'], 'pz': df['lep1_pz'], 'E': df['lep1_E']})
    #     lep2_p4 = vector.zip({'px': df['lep2_px'], 'py': df['lep2_py'], 'pz': df['lep2_pz'], 'E': df['lep2_E']})
    #     same_flavor = df['lep1_legType'] == df['lep2_legType']
    #     mll = (lep1_p4 + lep2_p4).mass.to_numpy()
    #     z_peak_veto = np.logical_and(same_flavor, np.abs(mll - 90.0) < 10.0)
    #     # selection = np.logical_and(selection, mll < 12.0)
    #     selection = np.logical_and(selection, z_peak_veto)
    #     return selection

    X_bkg, input_names, y_bkg, target_names = dataloader.Get()
    ys_pred_bkg = PreparePredictions(model, training_params, X_bkg)
    central_bkg = ys_pred_bkg[:, :, 1]
    assert y_bkg.shape == central_bkg.shape

    pred_errors_bkg = ys_pred_bkg[:, :, -1] - ys_pred_bkg[:, :, 0]
    if training_params['use_energy_layer']:
        pred_errors_bkg = np.delete(pred_errors_bkg, [3, 7], axis=1)

    ep = ErrorPropagator(global_corr_mtrx, central_bkg[:, :3], central_bkg[:, 4:7])
    prop_errors_bkg = ep.Propagate(pred_errors_bkg)

    bin_left = np.min(prop_errors_bkg) - 1.0
    bin_right = np.max(prop_errors_bkg) + 1.0
    bins = np.linspace(bin_left, bin_right, 50)
    PlotHist(data=prop_errors_bkg, 
             bins=bins,
             title="Predicted X->HH mass errors",
             ylabel='Count',
             xlabel='Mass error, [GeV]',
             plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
             file_name='bkg_mx_error')

    PlotHistStack([prop_errors_sig, prop_errors_bkg], 
                  {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
                   f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
                  f'sig_vs_{bkg_name}_mx_error',
                  bins=np.linspace(0, 250, 50),
                  plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
                  density=True,
                  xlabel='Mass error, [GeV]',
                  ylabel='Density',
                  title='Signal vs Background mass errors')


    mass_sig = ComputeMass(central_sig)
    mass_bkg = ComputeMass(central_bkg)

    rel_sig_err = prop_errors_sig/mass_sig
    rel_bkg_err = prop_errors_bkg/mass_bkg

    PlotHistStack([rel_sig_err, rel_bkg_err], 
                  {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
                   f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
                  f'sig_vs_{bkg_name}_rel_mx_error',
                  val_range=(0, 0.5),
                  plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
                  density=True,
                  xlabel='Mass relative error',
                  ylabel='Density',
                  title='Signal vs Background relative mass errors')

    mass_low = np.quantile(mass_sig, 0.16)
    mass_high = np.quantile(mass_sig, 0.84)
    sig_core_mask = np.logical_and(mass_sig > mass_low, mass_sig < mass_high)
    bkg_core_mask = np.logical_and(mass_bkg > mass_low, mass_bkg < mass_high)

    PlotHistStack([prop_errors_sig[sig_core_mask], prop_errors_bkg[bkg_core_mask]], 
                  {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
                   f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
                  f'core_sig_vs_{bkg_name}_mx_error',
                  bins=np.linspace(0, 250, 50),
                  plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
                  density=True,
                  xlabel='Mass error, [GeV]',
                  ylabel='Density',
                  title='Signal vs Background mass errors')

    PlotHistStack([rel_sig_err[sig_core_mask], rel_bkg_err[bkg_core_mask]], 
                  {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
                   f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
                  f'core_sig_vs_{bkg_name}_rel_mx_error',
                  val_range=(0, 0.5),
                  plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
                  density=True,
                  xlabel='Mass relative error',
                  ylabel='Density',
                  title='Signal vs Background mass errors')

if __name__ == '__main__':
    main()