import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt
import vector


from Dataloader import Dataloader
from NetUtils import TrainableSiLU
from LayerUtils import PtLayer, EtaLayer, PhiLayer, EnergyErrorLayer
from LossUtils import DeltaPhiLoss, CombinedEnergyLoss, LogPtLoss


from MiscUtils import *
from PlotUtils import PlotMetric, PlotHist


def main():
    files = []
    with open('files_Run3_2022.txt', 'r') as file_cfg:
        files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(files)
    X_train, input_names, y_train, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
    os.environ['TF_DETERMINISTIC_OPS'] = '1'
    tf.random.set_seed(42)

    model_name = "model_hh_dl_curv"
    model_dir = model_name
    os.makedirs(model_dir, exist_ok=True)

    learning_rate_scheduler = tf.keras.callbacks.LearningRateScheduler(Scheduler)

    input_shape = X_train.shape[1:]
    num_units = 384
    num_layers = 12
    l2_strength = 0.01
    epochs = 30
    batch_size = 256

    # prepare base part of the model
    inputs = tf.keras.layers.Input(shape=input_shape)
    x = tf.keras.layers.Dense(num_units, 
                              activation='silu', 
                              kernel_initializer='random_normal', 
                              bias_initializer='random_normal',
                              kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(inputs)
    for _ in range(num_layers - 1):
        x = tf.keras.layers.Dense(num_units, 
                                  activation='silu', 
                                  kernel_initializer='random_normal', 
                                  bias_initializer='random_normal',
                                  kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(x)

    # separate heads to H->bb and H->VV   
    hbb_head = tf.keras.layers.Dense(4)(x)
    hvv_head = tf.keras.layers.Dense(4)(x)

    # compute output variables and assign losses
    losses = {}

    # prepare H->bb outputs
    hbb_pt_head = PtLayer(name='genHbb_pt')(hbb_head)
    losses['genHbb_pt'] = LogPtLoss()
    hbb_eta_head = EtaLayer(name='genHbb_eta')(hbb_head)
    losses['genHbb_eta'] = tf.keras.losses.LogCosh()
    hbb_phi_head = PhiLayer(name='genHbb_phi')(hbb_head)
    losses['genHbb_phi'] = DeltaPhiLoss()
    hbb_en_head = EnergyErrorLayer(name='genHbb_E')(hbb_head)
    losses['genHbb_E'] = CombinedEnergyLoss(log=True)

    # prepare H->VV outputs
    hvv_pt_head = PtLayer(name='genHVV_pt')(hvv_head)
    losses['genHVV_pt'] = LogPtLoss()
    hvv_eta_head = EtaLayer(name='genHVV_eta')(hvv_head)
    losses['genHVV_eta'] = tf.keras.losses.LogCosh()
    hvv_phi_head = PhiLayer(name='genHVV_phi')(hvv_head)
    losses['genHVV_phi'] = DeltaPhiLoss()
    hvv_en_head = EnergyErrorLayer(name='genHVV_E')(hvv_head)
    losses['genHVV_E'] = CombinedEnergyLoss(log=True)
    
    outputs = [hvv_pt_head, hvv_eta_head, hvv_phi_head, hvv_en_head,
               hbb_pt_head, hbb_eta_head, hbb_phi_head, hbb_en_head]

    # assign weights to losses: not working
    # for some reason at the moment applying weights drags loss to nan
    loss_weights = {'genHbb_pt': 1.0,
                    'genHbb_eta': 1.0,
                    'genHbb_phi': 1.0,
                    'genHbb_E': 0.1,
                    'genHVV_pt': 1.0,
                    'genHVV_eta': 1.0,
                    'genHVV_phi': 1.0,
                    'genHVV_E': 0.1}

    ys_train = []
    for idx, target_name in enumerate(target_names):
        target_array = y_train[:, idx]
        target_array = np.reshape(target_array, (-1, 1))
        ys_train.append(target_array)

    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    model.compile(loss=losses,
                #   loss_weights=loss_weights,
                  optimizer=tf.keras.optimizers.Adam(3e-4))

    tf.keras.utils.plot_model(model, os.path.join(model_dir, f"summary_{model.name}.pdf"), show_shapes=True)

    history = model.fit(X_train,
                        ys_train,
                        validation_split=0.2,
                        verbose=1,
                        batch_size=batch_size,
                        epochs=epochs,
                        shuffle=True,
                        callbacks=[learning_rate_scheduler])

    tf.keras.utils.plot_model(model, os.path.join(model_dir, f"summary_{model.name}.pdf"), show_shapes=True)
    model.save(os.path.join(model_dir, f"{model_name}.keras"))

    history_keys = list(history.history.keys())
    drawable_metrics = [key for key in history_keys if 'loss' in key and 'val' not in key]
    for metric in drawable_metrics:
        PlotMetric(history, model_name, metric, plotting_dir=model_dir)    

    # form testing dataset
    def TestSelection(df, mod, parity, mass, sample_type):
        tmp = np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)
        return np.logical_and(tmp, df['sample_type'] == sample_type)

    X_test, _, y_test, _ = dataloader.Get(TestSelection, 2, 1, 800, 1)
    ys_pred = model.predict(X_test)

    assert len(ys_pred) == len(target_names), f"mismatch between number of predicted values and ground truth values: ({len(ys_pred)} vs {len(target_names)})"

    #plot predicted errors
    var_dict = {}
    for i, name in enumerate(target_names):
        ground_truth = y_test[:, i]
        
        if 'E' in name: 
            PlotCompare2D(ground_truth, ys_pred[i][:, 0], name, plotting_dir=model_dir)

            err = ys_pred[i][:, 1]
            bins = np.linspace(np.min(err) - 10, np.max(err) + 10, 50)
            PlotHist(data=err, 
                    bins=bins,
                    title=f'Predicted {name} error',
                    ylabel='Count',
                    xlabel='Error, [GeV]',
                    plotting_dir=model_dir,
                    file_name=f'pred_error_{name}')
            
            var_dict[name] = ys_pred[i][:, 0]
        else:
            var_dict[name] = ys_pred[i]
            PlotCompare2D(ground_truth, ys_pred[i], name, plotting_dir=model_dir)

    hbb_p4 = vector.zip({name.split('_')[-1]: var_dict[name] for name in var_dict.keys() if 'Hbb' in name})
    hvv_p4 = vector.zip({name.split('_')[-1]: var_dict[name] for name in var_dict.keys() if 'HVV' in name})
    X_p4 = hbb_p4 + hvv_p4
    pred_mass = np.array(X_p4.mass)
    
    PlotHist(data=pred_mass, 
             bins=np.linspace(0, 2500, 100),
            #  peak=True,
            #  width=True,
            #  count=True,
             title='Predicted X->HH error',
             ylabel='Count',
             xlabel='Mass, [GeV]',
             plotting_dir=model_dir,
             file_name='pred_mass')


if __name__ == '__main__':
    main()