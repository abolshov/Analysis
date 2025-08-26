import tensorflow as tf
import pandas as pd
import numpy as np
import vector
import awkward as ak
import os
import uproot

from Dataloader import Dataloader
from ErrorProp import ErrorPropagator

from MiscUtils import *
from PlotUtils import PlotHist, PlotCompare2D, PlotCovarMtrx, PlotHistStack
from ModelUtils import LoadModel
from PrepUtils import *

def main():
    channel = 'DL'
    masspoint = 350
    # load signal data
    suffix = '2B2Vto2B2JLNu' if channel == 'SL' else '2B2Vto2B2L2Nu'
    file = f'../train_data/Run3_2022/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)

    def TestSelection(df, mod, parity, sample_type):
        return np.logical_and(df['event'] % mod, df['sample_type'] == sample_type)

    X_sig, input_names, y_sig, target_names = dataloader.Get(TestSelection, 2, 1, 1)

    # load model
    model, training_params = LoadModel('model_cfg.yaml')
    print(model.summary())

    ys_pred_sig = PreparePredictions(model, training_params, X_sig)

    # extract prediction for central quantile for all targets
    central_sig = ys_pred_sig[:, :, 1]
    assert y_sig.shape == central_sig.shape

    pred_mass = ComputeMass(central_sig)
    with uproot.recreate(os.path.join(training_params['model_dir'], f'nn_out_{channel}_M{masspoint}.root')) as output_file:
        df = dataloader.Get(TestSelection, 2, 1, 1, df=True)
        vis_mass = ComputeVisMass(df, channel)

        data = {'nn_mass': pred_mass,
                'vis_mass': vis_mass,   
                'event': df['event']}

        output_file['Events'] = data

    # compute errors as up quantile (-1) minus down quantile (0)
    pred_errors_sig = ys_pred_sig[:, :, -1] - ys_pred_sig[:, :, 0]
    true_errors = y_sig - central_sig

    if training_params['use_energy_layer']:
        pred_errors_sig = np.delete(pred_errors_sig, [3, 7], axis=1)

    global_corr_mtrx = np.array(training_params['global_corr_mtrx'])
    ep = ErrorPropagator(global_corr_mtrx, central_sig[:, :3], central_sig[:, 4:7])
    prop_errors_sig = ep.Propagate(pred_errors_sig)

    pretty_labels = [ground_truth_map[name] for name in target_names if 'E' not in name]
    # PlotCovarMtrx(global_corr_mtrx, 'empirical', pretty_labels, os.path.join(training_params['model_dir'], 'plots'))

    # bin_left = np.min(prop_errors_sig) - 1.0
    # bin_right = np.max(prop_errors_sig) + 1.0
    # bins = np.linspace(bin_left, bin_right, 50)
    # PlotHist(data=prop_errors_sig, 
    #          bins=bins,
    #          title="Predicted X->HH mass errors",
    #          ylabel='Count',
    #          xlabel='Mass error, [GeV]',
    #          plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
    #          file_name='sig_mx_error')


    # # load bkg data
    # file = 'TTbar.root'
    # # file = 'DY.root'
    # bkg_latex = '$T\\bar{T}$'
    # bkg_name = 'ttbar'
    # dataloader = Dataloader('dataloader_config.yaml')
    # dataloader.Load(file)

    # X_bkg, input_names, y_bkg, target_names = dataloader.Get()
    # ys_pred_bkg = PreparePredictions(model, training_params, X_bkg)
    # central_bkg = ys_pred_bkg[:, :, 1]
    # assert y_bkg.shape == central_bkg.shape

    # pred_errors_bkg = ys_pred_bkg[:, :, -1] - ys_pred_bkg[:, :, 0]
    # if training_params['use_energy_layer']:
    #     pred_errors_bkg = np.delete(pred_errors_bkg, [3, 7], axis=1)

    # ep = ErrorPropagator(global_corr_mtrx, central_bkg[:, :3], central_bkg[:, 4:7])
    # prop_errors_bkg = ep.Propagate(pred_errors_bkg)

    # bin_left = np.min(prop_errors_bkg) - 1.0
    # bin_right = np.max(prop_errors_bkg) + 1.0
    # bins = np.linspace(bin_left, bin_right, 50)
    # PlotHist(data=prop_errors_bkg, 
    #          bins=bins,
    #          title="Predicted X->HH mass errors",
    #          ylabel='Count',
    #          xlabel='Mass error, [GeV]',
    #          plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
    #          file_name='bkg_mx_error')

    # PlotHistStack([prop_errors_sig, prop_errors_bkg], 
    #               {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
    #                f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
    #               f'sig_vs_{bkg_name}_mx_error',
    #               bins=np.linspace(0, 250, 50),
    #               plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
    #               density=True,
    #               xlabel='Mass error, [GeV]',
    #               ylabel='Density',
    #               title='Signal vs Background mass errors')


    # mass_sig = ComputeMass(central_sig)
    # mass_bkg = ComputeMass(central_bkg)

    # rel_sig_err = prop_errors_sig/mass_sig
    # rel_bkg_err = prop_errors_bkg/mass_bkg

    # PlotHistStack([rel_sig_err, rel_bkg_err], 
    #               {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
    #                f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
    #               f'sig_vs_{bkg_name}_rel_mx_error',
    #               val_range=(0, 0.5),
    #               plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
    #               density=True,
    #               xlabel='Mass relative error',
    #               ylabel='Density',
    #               title='Signal vs Background relative mass errors')

    # mass_low = np.quantile(mass_sig, 0.16)
    # mass_high = np.quantile(mass_sig, 0.84)
    # sig_core_mask = np.logical_and(mass_sig > mass_low, mass_sig < mass_high)
    # bkg_core_mask = np.logical_and(mass_bkg > mass_low, mass_bkg < mass_high)

    # PlotHistStack([prop_errors_sig[sig_core_mask], prop_errors_bkg[bkg_core_mask]], 
    #               {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
    #                f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
    #               f'core_sig_vs_{bkg_name}_mx_error',
    #               bins=np.linspace(0, 250, 50),
    #               plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
    #               density=True,
    #               xlabel='Mass error, [GeV]',
    #               ylabel='Density',
    #               title='Signal vs Background mass errors')

    # PlotHistStack([rel_sig_err[sig_core_mask], rel_bkg_err[bkg_core_mask]], 
    #               {'Radion $M_X=800$': {'linewidth': 2, 'color': 'blue'},
    #                f'{bkg_latex}': {'linewidth': 2, 'color': 'red'}},
    #               f'core_sig_vs_{bkg_name}_rel_mx_error',
    #               val_range=(0, 0.5),
    #               plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
    #               density=True,
    #               xlabel='Mass relative error',
    #               ylabel='Density',
    #               title='Signal vs Background mass errors')

if __name__ == '__main__':
    main()