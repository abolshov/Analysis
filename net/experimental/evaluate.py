import tensorflow as tf
import pandas as pd
import numpy as np
import os
import yaml
from sklearn.preprocessing import StandardScaler

from Dataloader import Dataloader
from ErrorProp import ErrorPropagator

from MiscUtils import *
from PlotUtils import PlotHist, PlotCompare2D, PlotCovarMtrx
from ModelUtils import LoadModel

def PreparePredictions(model, training_params, X_test):
    ys_pred = None
    if training_params['standardize']:
        input_means = training_params['input_train_means']
        input_scales = training_params['input_train_scales']
        X_test -= input_means
        X_test /= input_scales
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

def main():
    # load signal data
    file = 'nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)

    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)
        
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

    # pretty_labels = [ground_truth_map[name] for name in target_names if 'E' not in name]
    # PlotCovarMtrx(global_corr_mtrx, 'empirical', pretty_labels, os.path.join(training_params['model_dir'], 'plots'))

    bin_left = np.min(prop_errors_sig) - 1.0
    bin_right = np.max(prop_errors_sig) + 1.0
    bins = np.linspace(bin_left, bin_right, 50)
    PlotHist(data=pred_errors_sig, 
             bins=bins,
             title="Predicted X->HH mass errors",
             ylabel='Count',
             xlabel='Mass error, [GeV]',
             plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
             file_name='sig_mx_error')


    # load bkg data
    file = 'DY.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)

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

    plt.hist(pred_errors_sig, bins=np.linspace(0, 250, 50), density=True, color='blue', label='sig', histtype='step', linewidth=3)
    plt.hist(prop_errors_bkg, bins=np.linspace(0, 250, 50), density=True, color='red', label='bkg', histtype='step', linewidth=3)
    plt.xlabel('Mass error, [GeV]')
    plt.ylabel('Density')
    plt.title('Signal vs Background Mass errors')
    plt.legend()
    plt.grid()
    plt.savefig('sig_vs_bkg_mx_error.pdf', bbox_inches='tight')

if __name__ == '__main__':
    main()