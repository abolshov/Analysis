import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from DataWrapper import DataWrapper


def CustomLoss(y_true, y_pred):
    log_cosh = tf.keras.losses.log_cosh(y_true, y_pred)

    hbb_p3 = y_pred[:, 0:3]
    hbb_E = y_pred[:, 3]
    mh_bb_sqr = hbb_E*hbb_E - tf.reduce_sum(hbb_p3**2, axis=1, keepdims=True)

    hww_p3 = y_pred[:, 4:7]
    hww_E = y_pred[:, 7]
    mh_ww_sqr = hww_E*hww_E - tf.reduce_sum(hww_p3**2, axis=1, keepdims=True)

    mh_bb = tf.sign(mh_bb_sqr)*tf.sqrt(tf.abs(mh_bb_sqr))
    mh_ww = tf.sign(mh_ww_sqr)*tf.sqrt(tf.abs(mh_ww_sqr))

    custom = 0.5*tf.keras.losses.log_cosh(mh_bb, mh) + 0.5*tf.keras.losses.log_cosh(mh_ww, mh)
    return log_cosh + 0.1*custom


def PredWidth(pred_mass):
    q_84 = np.quantile(pred_mass, 0.84)
    q_16 = np.quantile(pred_mass, 0.16)
    width = q_84 - q_16
    return width 


def PredPeak(pred_mass):
    counts = np.bincount(pred_mass)
    peak = np.argmax(counts)
    return peak 


def PlotCompare(target, output, quantity):
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist(target, bins=bins, label=f'target_{quantity}')
    plt.hist(output, bins=bins, label=f'output_{quantity}')
    plt.title(f'{quantity}')
    plt.ylabel('Count')
    plt.xlabel(f'{quantity}')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig(f"cmp_{quantity}.pdf", bbox_inches='tight')
    plt.clf()


def PlotCompare2D(target, output, quantity):
    plt.grid(False)
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist2d(target, output, bins=bins)
    plt.title(f'{quantity}')
    plt.ylabel(f'output_{quantity}')
    plt.xlabel(f'target_{quantity}')
    plt.savefig(f"cmp2d_{quantity}.pdf", bbox_inches='tight')
    plt.clf()

# model = tf.keras.models.load_model("model_hh_custloss.keras", custom_objects={"CustomLoss": CustomLoss})
model = tf.keras.models.load_model("model_hh_custloss.keras")
print(model.summary())

dw = DataWrapper()
dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
dw.FormTestSet()
res = model.predict(dw.test_features)
pred_dict = {label: res[:, i] for i, label in enumerate(dw.test_labels.columns)}
pred_df = pd.DataFrame.from_dict(pred_dict)

for col in pred_df.columns:
    # PlotCompare(dw.test_labels[col], pred_df[col], col)
    PlotCompare2D(dw.test_labels[col], pred_df[col], col)

X_E_pred = pred_df["genHbb_E"] + pred_df["genHVV_E"]
X_px_pred = pred_df["genHbb_px"] + pred_df["genHVV_px"]
X_py_pred = pred_df["genHbb_py"] + pred_df["genHVV_py"]
X_pz_pred = pred_df["genHbb_pz"] + pred_df["genHVV_pz"]

Hbb_en = np.sqrt(125.0**2 + pred_df["genHbb_px"]**2 + pred_df["genHbb_py"]**2 + pred_df["genHbb_pz"]**2)
Hww_en = np.sqrt(125.0**2 + pred_df["genHVV_px"]**2 + pred_df["genHVV_py"]**2 + pred_df["genHVV_pz"]**2)
X_en = Hbb_en + Hww_en
X_mass = np.sqrt(X_en**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2)
bins = np.linspace(200, 2000, 100)
plt.hist(X_mass, bins=bins)
plt.title('Predicted X->HH mass with mH = 125.0 GeV constraint')
plt.ylabel('Count')
plt.xlabel(f'mass, [GeV]')
plt.figtext(0.75, 0.8, f"peak: {PredPeak(X_mass):.2f}")
plt.figtext(0.75, 0.75, f"width: {PredWidth(X_mass):.2f}")
plt.grid(True)
plt.savefig(f"constr_pred_mass.pdf", bbox_inches='tight')
plt.clf()

Hbb_mass_sqr = pred_df["genHbb_E"]**2 - pred_df["genHbb_px"]**2 - pred_df["genHbb_py"]**2 - pred_df["genHbb_pz"]**2

pos_frac = len(Hbb_mass_sqr[Hbb_mass_sqr > 0.0])/(len(Hbb_mass_sqr))
Hbb_mass = np.sqrt(Hbb_mass_sqr[Hbb_mass_sqr > 0.0])

bins = np.linspace(0, 250, 100)
plt.hist(Hbb_mass, bins=bins)
plt.title('Predicted H->bb mass')
plt.ylabel('Count')
plt.xlabel(f'mass, [GeV]')
plt.figtext(0.75, 0.8, f"peak: {PredPeak(Hbb_mass):.2f}")
plt.figtext(0.75, 0.75, f"width: {PredWidth(Hbb_mass):.2f}")
plt.figtext(0.75, 0.7, f"pos: {pos_frac:.2f}")
plt.grid(True)
plt.savefig(f"Hbb_mass.pdf", bbox_inches='tight')
plt.clf()

Hvv_mass_sqr = pred_df["genHVV_E"]**2 - pred_df["genHVV_px"]**2 - pred_df["genHVV_py"]**2 - pred_df["genHVV_pz"]**2
pos_frac = len(Hvv_mass_sqr[Hvv_mass_sqr > 0.0])/(len(Hvv_mass_sqr))
Hvv_mass = np.sqrt(Hvv_mass_sqr[Hvv_mass_sqr > 0.0])

plt.hist(Hvv_mass, bins=bins)
plt.title('Predicted H->VV mass')
plt.ylabel('Count')
plt.xlabel(f'mass, [GeV]')
plt.figtext(0.75, 0.8, f"peak: {PredPeak(Hvv_mass):.2f}")
plt.figtext(0.75, 0.75, f"width: {PredWidth(Hvv_mass):.2f}")
plt.figtext(0.75, 0.7, f"pos: {pos_frac:.2f}")
plt.grid(True)
plt.savefig(f"Hvv_mass.pdf", bbox_inches='tight')
plt.clf()

# eventIds = np.array(dw.test_events.astype(int))
# pred_df["eventId"] = eventIds
# print(pred_df.head())

# true_df = dw.test_labels
# true_df["eventIds"] = np.array(dw.test_events.astype(int))

# np.savetxt("pred_data_logcosh.txt", pred_df.values, fmt='%10.5f')
# np.savetxt("true_data.txt", true_df.values, fmt='%10.5f')

X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
print(f"positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
bins = np.linspace(200, 2000, 100)
plt.hist(X_mass_pred, bins=bins)
plt.title('Predicted X->HH mass')
plt.ylabel('Count')
plt.xlabel(f'mass, [GeV]')
plt.figtext(0.75, 0.8, f"peak: {PredPeak(X_mass_pred):.2f}")
plt.figtext(0.75, 0.75, f"width: {PredWidth(X_mass_pred):.2f}")
plt.grid(True)
plt.savefig(f"pred_mass.pdf", bbox_inches='tight')
plt.clf()

# pred_offshellW_E = np.sqrt(pred_df["offshellV_mass"]**2
#                           + pred_df["offshellV_px"]**2 
#                           + pred_df["offshellV_py"]**2 
#                           + pred_df["offshellV_pz"]**2)

# true_offshellW_E = np.sqrt(dw.test_labels["offshellV_mass"]**2
#                            + dw.test_labels["offshellV_px"]**2 
#                            + dw.test_labels["offshellV_py"]**2 
#                            + dw.test_labels["offshellV_pz"]**2)

# PlotCompare(true_offshellW_E, pred_offshellW_E, "offshellW_E")
# PlotCompare2D(true_offshellW_E, pred_offshellW_E, "offshellW_E")

# bins = np.linspace(0, 1000, 100)
# plt.hist(pred_offshellW_E, bins=bins)
# plt.title('Predicted energy of offshell W')
# plt.ylabel('Count')
# plt.xlabel(f'energy, [GeV]')
# plt.grid(True)
# plt.savefig(f"pred_energy.pdf", bbox_inches='tight')
# plt.clf()

# bins = np.linspace(0, 100, 50)
# plt.hist(pred_df["offshellV_mass"], bins=bins, label="pred")
# plt.hist(dw.test_labels["offshellV_mass"], bins=bins, label="true")
# plt.title('Predicted mass of offshell W')
# plt.ylabel('Count')
# plt.xlabel(f'mass, [GeV]')
# plt.grid(True)
# plt.legend(loc='upper right')
# plt.savefig(f"pred_ofshell_mass.pdf", bbox_inches='tight')
# plt.clf()

# bins = np.linspace(0, 1000, 100)
# plt.hist(pred_offshellW_E, bins=bins, label="pred")
# plt.hist(true_offshellW_E, bins=bins, label="true")
# plt.title('True E vs predicted E')
# plt.ylabel('Count')
# plt.xlabel(f'energy, [GeV]')
# plt.grid(True)
# plt.legend(loc='upper right')
# plt.savefig(f"energy_cmp.pdf", bbox_inches='tight')
# plt.clf()

# offshellW_mass_sqr = np.array(pred_df["offshellV_E"]**2
#                               - pred_df["offshellV_px"]**2 
#                               - pred_df["offshellV_py"]**2 
#                               - pred_df["offshellV_pz"]**2)

# bins = np.linspace(-1000, 1000, 100)
# plt.hist(offshellW_mass_sqr, bins=bins)
# plt.title('Predicted square of offshell W mass')
# plt.ylabel('Count')
# plt.xlabel(f'mass^2, [GeV^2]')
# plt.grid(True)
# plt.savefig(f"pred_mass_sqr.pdf", bbox_inches='tight')
# plt.clf()

# # offshellW_mass_sqr = offshellW_mass_sqr[offshellW_mass_sqr > 0.0]
# offshellW_mass = np.sqrt(np.abs(offshellW_mass_sqr))
# # offshellW_mass = offshellW_mass[offshellW_mass > 0.0]

# bins = np.linspace(0, 200, 100)
# plt.hist(offshellW_mass, bins=bins)
# plt.title('Predicted offshell W mass')
# plt.ylabel('Count')
# plt.xlabel(f'mass, [GeV]')
# # plt.figtext(0.75, 0.8, f"peak: {PredPeak(offshellW_mass):.2f}")
# # plt.figtext(0.75, 0.75, f"width: {PredWidth(offshellW_mass):.2f}")
# plt.grid(True)
# plt.savefig(f"pred_mass.pdf", bbox_inches='tight')
# plt.clf()

# print("true labels:")
# print(dw.test_labels.head())

# print("pred labels:")
# print(pred_df.head())

# print("pred mass:")
# print(offshellW_mass_sqr[0:5])

# tue_mass_sqr = np.array(dw.test_labels["offshellV_E"]**2
#                         - dw.test_labels["offshellV_px"]**2 
#                         - dw.test_labels["offshellV_py"]**2 
#                         - dw.test_labels["offshellV_pz"]**2)

# print("true mass:")
# print(tue_mass_sqr[0:5])