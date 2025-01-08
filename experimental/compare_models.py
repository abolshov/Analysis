import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from DataWrapper import DataWrapper


m1_name = "model_hh_augmented_new"
m2_name = "model_hh_augmented_old"
m1 = tf.keras.models.load_model(f"{m1_name}.keras")
m2 = tf.keras.models.load_model(f"{m2_name}.keras")

dw = DataWrapper()
dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
dw.FormTestSet()

res1 = m1.predict(dw.test_features)
pred_dict_1 = {label: res1[:, i] for i, label in enumerate(dw.test_labels.columns)}
pred_df_1 = pd.DataFrame.from_dict(pred_dict_1)

res2 = m2.predict(dw.test_features)
pred_dict_2 = {label: res2[:, i] for i, label in enumerate(dw.test_labels.columns)}
pred_df_2 = pd.DataFrame.from_dict(pred_dict_2)

E1 = pred_df_1["genHbb_E"] + pred_df_1["genHVV_E"]
px1 = pred_df_1["genHbb_px"] + pred_df_1["genHVV_px"]
py1 = pred_df_1["genHbb_py"] + pred_df_1["genHVV_py"]
pz1 = pred_df_1["genHbb_pz"] + pred_df_1["genHVV_pz"]
mass1 = np.sqrt(E1**2 - px1**2 - py1**2 - pz1**2)
print(f"mass1_median={np.median(mass1):.2f}")

E2 = pred_df_2["genHbb_E"] + pred_df_2["genHVV_E"]
px2 = pred_df_2["genHbb_px"] + pred_df_2["genHVV_px"]
py2 = pred_df_2["genHbb_py"] + pred_df_2["genHVV_py"]
pz2 = pred_df_2["genHbb_pz"] + pred_df_2["genHVV_pz"]
mass2 = np.sqrt(E2**2 - px2**2 - py2**2 - pz2**2)
print(f"mass2_median={np.median(mass2):.2f}")

bins = np.linspace(200, 2000, 100)
plt.hist(mass1, bins=bins, label=f'{m1_name}', alpha=0.75)
plt.hist(mass2, bins=bins, label=f'{m2_name}', alpha=0.75)
plt.title('Model comparison')
plt.ylabel('Count')
plt.xlabel(f'mass, [GeV]')
plt.legend(loc='upper right')
plt.grid(True)
plt.savefig(f"{m1_name}_vs_{m2_name}.pdf", bbox_inches='tight')
plt.clf()