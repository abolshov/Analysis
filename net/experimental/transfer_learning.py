import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os

from DataWrapper import DataWrapper
from NetUtils import TrainableSiLU, MassLoss, CombinedLoss
from PlotUtils import PlotCompare2D, PlotHist, PlotMetric

def main():
    path_to_pretrained = "model_hh_dl_resolved_aug_train_silu_l2reg/model_hh_dl_resolved_aug_train_silu_l2reg.keras"
    custom_objects = {"TrainableSiLU": TrainableSiLU}
    base_model = tf.keras.models.load_model(path_to_pretrained, custom_objects=custom_objects)

    base_model.trainable = False

    input_file_names = []
    with open("dl_train_files.txt", 'r') as file:
        input_file_names = [line.rstrip() for line in file.readlines()]

    dw = DataWrapper()
    dw.ReadFiles(input_file_names)
    dw.AugmentDataset()
    dw.FormTrainSet()

    input_shape = dw.train_features.shape
    num_units = 8

    inputs = tf.keras.Input(shape=input_shape[1:])
    x = base_model(inputs, training=False)
    x = tf.keras.layers.Dense(units=num_units, 
                              activation='relu',
                              kernel_initializer='random_normal', 
                              bias_initializer='random_normal')(x)
    x = tf.keras.layers.Dense(units=num_units, 
                              activation='relu',
                              kernel_initializer='random_normal', 
                              bias_initializer='random_normal')(x)
    x = tf.keras.layers.Dense(units=num_units, 
                              activation='relu',
                              kernel_initializer='random_normal', 
                              bias_initializer='random_normal')(x)
    outputs = tf.keras.layers.Dense(units=num_units)(x)
    composite_model = tf.keras.Model(inputs, outputs)

    composite_model.compile(loss='mse', 
                            optimizer=tf.keras.optimizers.Adam(3e-4))

    history = composite_model.fit(dw.train_features,
                                  dw.train_labels,
                                  validation_split=0.2,
                                  verbose=1,
                                  batch_size=1024,
                                  epochs=30)

    transfer_learning_dir = "transfer_learning"
    model_name = "model_hh_dl_transf_relu"
    model_dir = os.path.join(transfer_learning_dir, model_name)
    os.makedirs(model_dir, exist_ok=True)
    composite_model.save(os.path.join(model_dir, f"{model_name}.keras"))
    PlotMetric(history, model_name, "loss", plotting_dir=model_dir)
    tf.keras.utils.plot_model(composite_model, 
                            to_file=os.path.join(model_dir, f'{model_name}_arch.pdf'), 
                            show_shapes=True)

    dw.Clear()
    dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
    dw.FormTestSet()

    # print(composite_model.summary())

    # first_input = dw.test_features.iloc[0].to_numpy()
    # first_input = first_input[None, :]
    # print(base_model.predict(first_input)[0])
    # print(composite_model.predict(first_input)[0])
    # print(dw.test_labels.iloc[0].to_numpy())

    res = composite_model.predict(dw.test_features)
    pred_dict = {label: res[:, i] for i, label in enumerate(dw.test_labels.columns)}
    pred_df = pd.DataFrame.from_dict(pred_dict)

    for col in pred_df.columns:
        PlotCompare2D(dw.test_labels[col], pred_df[col], col, plotting_dir=model_dir)

    bins = np.linspace(0, 2000, 100)

    X_E_pred = pred_df["genHbb_E"] + pred_df["genHVV_E"]
    X_px_pred = pred_df["genHbb_px"] + pred_df["genHVV_px"]
    X_py_pred = pred_df["genHbb_py"] + pred_df["genHVV_py"]
    X_pz_pred = pred_df["genHbb_pz"] + pred_df["genHVV_pz"]

    X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
    X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
    pos_frac = len(X_mass_sqr_pos)/len(X_mass_pred_sqr)
    print(f"positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={pos_frac:.3f}")

    if pos_frac > 0.0:
        PlotHist(data=np.sqrt(X_mass_sqr_pos), 
                 title="Predicted X->HH mass",
                 bins=np.linspace(0, 2000, 100),
                 ylabel='Count',
                 xlabel='mass, [GeV]',
                 pos_frac=pos_frac,
                 plotting_dir=model_dir,
                 file_name='pred_mass')

    Hvv_mass_sqr = pred_df["genHVV_E"]**2 - pred_df["genHVV_px"]**2 - pred_df["genHVV_py"]**2 - pred_df["genHVV_pz"]**2
    pos_frac = len(Hvv_mass_sqr[Hvv_mass_sqr > 0.0])/(len(Hvv_mass_sqr))
    Hvv_mass = np.sqrt(Hvv_mass_sqr[Hvv_mass_sqr > 0.0])
    if pos_frac > 0.0:
        PlotHist(data=Hvv_mass, 
                 title="Predicted H->VV mass",
                 ylabel='Count',
                 xlabel='mass, [GeV]',
                 pos_frac=pos_frac,
                 plotting_dir=model_dir,
                 file_name='Hvv_mass')

    Hbb_mass_sqr = pred_df["genHbb_E"]**2 - pred_df["genHbb_px"]**2 - pred_df["genHbb_py"]**2 - pred_df["genHbb_pz"]**2
    pos_frac = len(Hbb_mass_sqr[Hbb_mass_sqr > 0.0])/(len(Hbb_mass_sqr))
    Hbb_mass = np.sqrt(Hbb_mass_sqr[Hbb_mass_sqr > 0.0])
    if pos_frac > 0.0:
        PlotHist(data=Hbb_mass, 
                 title="Predicted H->bb mass",
                 ylabel='Count',
                 xlabel='mass, [GeV]',
                 pos_frac=pos_frac,
                 plotting_dir=model_dir,
                 file_name='Hbb_mass')


if __name__ == '__main__':
    main()