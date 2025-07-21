import tensorflow as tf
import pandas as pd
import numpy as np
import os


from Dataloader import Dataloader
from NetUtils import *
from LossUtils import *
from LayerUtils import *


from MiscUtils import *
from PlotUtils import PlotHist


def main():
    # load dataset for model evaluation
    file = 'nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)
    X, input_names, y, target_names = dataloader.Get()

    # load model
    model_dir = 'model_test'
    model_name = 'model_test'
    path_to_model = os.path.join(model_dir, f'{model_name}.keras')
    custom_objects = {'TrainableSiLU': TrainableSiLU,
                      'QuantileLoss': QuantileLoss,
                      'QuantileOrderingLayer': QuantileOrderingLayer}
    model = tf.keras.models.load_model(path_to_model, custom_objects=custom_objects)

    # make predictions and plot them
    ys_pred = model.predict(X_test)
    
    pred_dict = {}
    for i, name in enumerate(target_names):
        pred_array = ys_pred[i]
        for q_idx, quantile in enumerate(quantiles):
            q_pred_array = pred_array[:, q_idx]
            q = int(100*quantile)
            pred_dict[f'{name}_q{q}'] = q_pred_array

    pred_df = pd.DataFrame.from_dict(pred_dict)

    bins = np.linspace(0, 2000, 100)
    for quantile in quantiles:
        q = int(100*quantile)

        X_E_pred = pred_df[f"genHbb_E_q{q}"] + pred_df[f"genHVV_E_q{q}"]
        X_px_pred = pred_df[f"genHbb_px_q{q}"] + pred_df[f"genHVV_px_q{q}"]
        X_py_pred = pred_df[f"genHbb_py_q{q}"] + pred_df[f"genHVV_py_q{q}"]
        X_pz_pred = pred_df[f"genHbb_pz_q{q}"] + pred_df[f"genHVV_pz_q{q}"]

        X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
        X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
        print(f"quantile {q}: positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

        X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
        PlotHist(data=X_mass_pred, 
                 bins=bins,
                 title="Predicted X->HH mass",
                 ylabel='Count',
                 xlabel='Mass, [GeV]',
                 plotting_dir=model_dir,
                 peak=True,
                 width=True,
                 count=True,
                 file_name=f'eval_mass_q{q}')

    # plot errors
    for i, name in enumerate(target_names):
        up = pred_df[f'{name}_q84']
        down = pred_df[f'{name}_q16']
        ground_truth = y_test[:, i]
        err = up - down
        pred = pred_df[f'{name}_q50']
        correct_predictions = np.abs(pred - ground_truth) < err
        proba = len(ground_truth[correct_predictions])/len(ground_truth)
        print(f'Target {name}:')
        print(f'\tprediction coverage probability: {proba:.2f}')
        print(f'\tquantile crossing fraction: {len(err[err < 0.0])/len(err):.2f}')

        bins = np.linspace(np.min(err) - 10, np.max(err) + 10, 50)
        PlotHist(data=err, 
                 bins=bins,
                 title=f'Predicted {name} error',
                 ylabel='Count',
                 xlabel='Error, [GeV]',
                 plotting_dir=model_dir,
                 file_name=f'eval_error_{name}')


if __name__ == '__main__':
    main()