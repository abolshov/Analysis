import tensorflow as tf
import pandas as pd
import numpy as np
import os
import yaml
from sklearn.preprocessing import StandardScaler
from sklearn.covariance import EmpiricalCovariance, MinCovDet

from Dataloader import Dataloader

from MiscUtils import *
from PlotUtils import PlotHist, PlotCompare2D, PlotCovarMtrx
from ModelUtils import LoadModel


def main():
    # load dataset for model evaluation
    # file = 'DY.root'
    file = '../train_data/Run3_2022/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M_800/nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)

    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)
        
    X, input_names, y, target_names = dataloader.Get(TestSelection, 2, 1, 800, 1)

    # load model
    model, training_params = LoadModel('model_cfg.yaml')
    print(model.summary())

    ys_pred = None
    if training_params['standardize']:
        input_means = training_params['input_train_means']
        input_scales = training_params['input_train_scales']
        X -= input_means
        X /= input_scales
        ys_pred = model.predict(X)

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

    # extract prediction for central quantile for all targets
    central = ys_pred[:, :, 1]
    assert y.shape == central.shape

    hvv_PtEtaPhiE_pred = ToPtEtaPhiE(central[:, :4])
    hvv_PtEtaPhiE_true = ToPtEtaPhiE(y[:, :4])

    hbb_PtEtaPhiE_pred = ToPtEtaPhiE(central[:, 4:])
    hbb_PtEtaPhiE_true = ToPtEtaPhiE(y[:, 4:])

    pred_PtEtaPhiE = np.concatenate([hvv_PtEtaPhiE_pred, hbb_PtEtaPhiE_pred], axis=-1)
    true_PtEtaPhiE = np.concatenate([hvv_PtEtaPhiE_true, hbb_PtEtaPhiE_true], axis=-1)

    var_names = ['genHVV_pt', 'genHVV_eta', 'genHVV_phi', 'genHVV_E', 
                 'genHbb_pt', 'genHbb_eta', 'genHbb_phi', 'genHbb_E']
    for i, name in enumerate(var_names):
        var = name.split('_')[-1]
        obj = name.split('_')[0]

        bin_left = np.min([np.quantile(true_PtEtaPhiE[:, i], 0.01), np.quantile(pred_PtEtaPhiE[:, i], 0.01)]) - 1.0
        bin_right = np.max([np.quantile(true_PtEtaPhiE[:, i], 0.99), np.quantile(pred_PtEtaPhiE[:, i], 0.99)]) + 1.0
        bins = np.linspace(bin_left, bin_right, 100)

        PlotCompare2D(true_PtEtaPhiE[:, i], 
                      pred_PtEtaPhiE[:, i], 
                      var,
                      obj,
                      bins=bins, 
                      title=f'{objects[obj] if obj in objects else obj} {pretty_vars[var] if var in pretty_vars else var} comparison',
                      xlabel=f'True {ground_truth_map[name]}',
                      ylabel=f'Predicted {ground_truth_map[name]}',
                      plotting_dir=os.path.join(training_params['model_dir'], 'plots'))

    # compute errors as up quantile (-1) minus down quantile (0)
    pred_errors = ys_pred[:, :, -1] - ys_pred[:, :, 0]
    true_errors = y - central

    n_events, n_targets = pred_errors.shape
    # sigma_mtrx will be used to transform all-dataset correlation matrix to n_events per-event covariance matrices
    # that are needed for error propagation
    I = np.eye(n_targets)[None, :, :] # shape (1, n_targets, n_targets)
    sigma_mtrx = pred_errors[:, :, None]*I # broadcast 

    # PROBLEM: for now some elements of sigma_mtrx are negative (either due to quantile crossing or bc it includes energy)
    # need to erase energy in case when higgs mass constraint was used for neural net training
    # erase when fixed
    assert np.all(sigma_mtrx > 0.0), "Diagonal matrix of variances for event contains negative variances"

    normalize_cov_mtrx = True
    # technically normalized covariance matrix is correlation matrix
    # to get it, compute covariance matrix of standardized dataset
    normalzied_true_errors = None
    if normalize_cov_mtrx:
        normalzied_true_errors = StandardScaler().fit_transform(true_errors)
    
    # actually this matrix must be computed across entire training dataset
    # and loaded here from training config
    method = 'empirical'
    global_covar = None
    if method == 'empirical':
        global_covar = EmpiricalCovariance().fit(normalzied_true_errors if normalize_cov_mtrx else true_errors)
    elif mehtod == 'mcd':
        global_covar = MinCovDet().fit(normalzied_true_errors if normalize_cov_mtrx else true_errors)
    else:
        raise RuntimeError(f'Illegal covariance matrix estimation method {method}')

    pretty_labels = [ground_truth_map[name] for name in target_names]
    global_covar_mtrx = global_covar.covariance_
    PlotCovarMtrx(global_covar_mtrx, method, pretty_labels, os.path.join(training_params['model_dir'], 'plots'))

    # covariance matrix for each event
    # will be used for error propagation
    event_covar_mtrx = np.sqrt(sigma_mtrx)*global_covar_mtrx*np.sqrt(sigma_mtrx)

    assert event_covar_mtrx.shape[1:] == global_covar_mtrx.shape

    # pred_dict = {}
    # quantiles = [0.16, 0.5, 0.84]
    # for i, name in enumerate(target_names):
    #     pred_array = ys_pred[i]
    #     for q_idx, quantile in enumerate(quantiles):
    #         q_pred_array = pred_array[:, q_idx]
    #         q = int(100*quantile)
    #         pred_dict[f'{name}_q{q}'] = q_pred_array

    # bins = np.linspace(0, 2000, 100)
    # for quantile in quantiles:
    #     q = int(100*quantile)

    #     X_E_pred = pred_df[f"genHbb_E_q{q}"] + pred_df[f"genHVV_E_q{q}"]
    #     X_px_pred = pred_df[f"genHbb_px_q{q}"] + pred_df[f"genHVV_px_q{q}"]
    #     X_py_pred = pred_df[f"genHbb_py_q{q}"] + pred_df[f"genHVV_py_q{q}"]
    #     X_pz_pred = pred_df[f"genHbb_pz_q{q}"] + pred_df[f"genHVV_pz_q{q}"]

    #     X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
    #     X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
    #     print(f"quantile {q}: positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

    #     X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
    #     PlotHist(data=X_mass_pred, 
    #              bins=bins,
    #              title="Predicted X->HH mass",
    #              ylabel='Count',
    #              xlabel='Mass, [GeV]',
    #              plotting_dir=model_dir,
    #              peak=True,
    #              width=True,
    #              count=True,
    #              file_name=f'eval_mass_q{q}')

    # # plot errors
    # for i, name in enumerate(target_names):
    #     up = pred_df[f'{name}_q84']
    #     down = pred_df[f'{name}_q16']
    #     err = up - down
    #     pred = pred_df[f'{name}_q50']
    #     print(f'Target {name}:')
    #     print(f'\tquantile crossing fraction: {len(err[err < 0.0])/len(err):.2f}')

    #     bins = np.linspace(np.min(err) - 10, np.max(err) + 10, 50)
    #     PlotHist(data=err, 
    #              bins=bins,
    #              title=f'Predicted {name} error',
    #              ylabel='Count',
    #              xlabel='Error, [GeV]',
    #              plotting_dir=model_dir,
    #              file_name=f'eval_error_{name}')


if __name__ == '__main__':
    main()