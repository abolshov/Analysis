import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from DataWrapper import DataWrapper
from NetUtils import TrainableSiLU, EpochLossUpdater
from LossUtils import Momentum3DLoss


target_names = {"genHbb_E": "E(H->bb)",
                "genHbb_px": "px(H->bb)",
                "genHbb_py": "py(H->bb)",
                "genHbb_pz": "pz(H->bb)",
                "genHVV_E": "E(H->VV)",
                "genHVV_px": "px(H->VV)",
                "genHVV_py": "py(H->VV)",
                "genHVV_pz": "pz(H->VV)" }


def PredWidth(pred_mass):
    q_84 = np.quantile(pred_mass, 0.84)
    q_16 = np.quantile(pred_mass, 0.16)
    width = q_84 - q_16
    return width 


def PredPeak(pred_mass):
    counts = np.bincount(pred_mass)
    peak = np.argmax(counts)
    return peak 


def PlotCompare2D(target, output, quantity, plotting_dir=None):
    plt.grid(False)
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist2d(target, output, bins=bins)
    var = quantity.split('_')[-1]
    plt.title(f'{var} comparison')
    plt.ylabel(f'predicted {target_names[quantity]}')
    plt.xlabel(f'true {target_names[quantity]}')
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f"cmp2d_{quantity}.pdf"), bbox_inches='tight')
    else:
        plt.savefig(f"cmp2d_{quantity}.pdf", bbox_inches='tight')
    plt.clf()


def Scheduler(epoch, lr):
    if epoch < 30:
        return lr
    else:
        if epoch % 2 == 0:
            return 0.9*lr
        return lr


def PlotMetric(history, model, metric, plotting_dir=None):
    plt.plot(history.history[metric], label=f'train_{metric}')
    plt.plot(history.history[f'val_{metric}'], label=f'val_{metric}')
    plt.title(f'{model} {metric}')
    plt.ylabel(metric)
    plt.xlabel('Epoch')
    plt.legend(loc='upper right')
    plt.grid(True)
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f"{metric}_{model}.pdf"), bbox_inches='tight')
    else:
        plt.savefig(f"{metric}_{model}.pdf", bbox_inches='tight')
    plt.clf()


input_file_names = []
with open("dl_train_files.txt", 'r') as file:
    input_file_names = [line.rstrip() for line in file.readlines()]

dw = DataWrapper()
dw.ReadFiles(input_file_names)
dw.AugmentDataset()
dw.FormTrainSet()

std = dw.label_scaler.scale_
mean = dw.label_scaler.mean_

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'
tf.random.set_seed(42)

model_name = "model_hh_dl_aug_train_silu_l2reg_mom3Dloss"
model_dir = model_name
os.makedirs(model_dir, exist_ok=True)

learning_rate_scheduler = tf.keras.callbacks.LearningRateScheduler(Scheduler)
early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor='val_loss',
    patience=3,
    min_delta=0.001
)

penalty_strength = tf.Variable(1e-6, dtype=tf.float32)
update_rate = 1.0001
loss_updater = EpochLossUpdater(penalty_strength, update_rate)

input_shape = dw.train_features.shape
num_units = 384
num_layers = 12
dropout_rate = 0.05
l2_strength = 0.01

model = tf.keras.Sequential()
for _ in range(num_layers):
    model.add(tf.keras.layers.Dense(num_units, 
                                    activation=TrainableSiLU(units=num_units), 
                                    kernel_initializer='random_normal', 
                                    bias_initializer='random_normal',
                                    kernel_regularizer=tf.keras.regularizers.L2(l2_strength)))
model.add(tf.keras.layers.Dense(6))

# model.compile(loss=CombinedLoss(penalty_strength), 
#               optimizer=tf.keras.optimizers.Adam(3e-4))
model.compile(loss=Momentum3DLoss(), 
              optimizer=tf.keras.optimizers.Adam(3e-4))
model.build(dw.train_features.shape)

# history = model.fit(dw.train_features,
#                     dw.train_labels,
#                     validation_split=0.2,
#                     verbose=1,
#                     batch_size=2048,
#                     epochs=100,
#                     callbacks=[learning_rate_scheduler, loss_updater])

history = model.fit(dw.train_features,
                    dw.train_labels,
                    validation_split=0.2,
                    verbose=1,
                    batch_size=2048,
                    epochs=100,
                    callbacks=[learning_rate_scheduler])

model.save(os.path.join(model_dir, f"{model_name}.keras"))
PlotMetric(history, model_name, "loss", plotting_dir=model_dir)

dw.Clear()
dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
# dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-800/nano_0.root")
dw.FormTestSet()

res = model.predict(dw.test_features)
# pred_dict = {label: res[:, i] for i, label in enumerate(dw.test_labels.columns)}
pred_dict = {label: res[:, i] for i, label in enumerate(["genHbb_px", "genHbb_py", "genHbb_pz", "genHVV_px", "genHVV_py", "genHVV_pz"])}
pred_df = pd.DataFrame.from_dict(pred_dict)
print(pred_df.describe())

Hbb_en = np.sqrt(125.0**2 + pred_df["genHbb_px"]**2 + pred_df["genHbb_py"]**2 + pred_df["genHbb_pz"]**2)
Hww_en = np.sqrt(125.0**2 + pred_df["genHVV_px"]**2 + pred_df["genHVV_py"]**2 + pred_df["genHVV_pz"]**2)
X_en = Hbb_en + Hww_en
pred_df["genHbb_E"] = Hbb_en
pred_df["genHVV_E"] = Hww_en

for col in pred_df.columns:
    PlotCompare2D(dw.test_labels[col], pred_df[col], col, plotting_dir=model_dir)

bins = np.linspace(0, 2000, 100)

X_E_pred = pred_df["genHbb_E"] + pred_df["genHVV_E"]
X_px_pred = pred_df["genHbb_px"] + pred_df["genHVV_px"]
X_py_pred = pred_df["genHbb_py"] + pred_df["genHVV_py"]
X_pz_pred = pred_df["genHbb_pz"] + pred_df["genHVV_pz"]

X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
print(f"positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
plt.hist(X_mass_pred, bins=bins)
plt.title('Predicted X->HH mass')
plt.ylabel('Count')
plt.xlabel(f'mass, [GeV]')
plt.figtext(0.75, 0.8, f"peak: {np.median(X_mass_pred):.2f}")
plt.figtext(0.75, 0.75, f"width: {PredWidth(X_mass_pred):.2f}")
plt.grid(True)
plt.savefig(os.path.join(model_dir, "pred_mass.pdf"), bbox_inches='tight')
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
plt.savefig(os.path.join(model_dir, "Hbb_mass.pdf"), bbox_inches='tight')
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
plt.savefig(os.path.join(model_dir, "Hvv_mass.pdf"), bbox_inches='tight')
plt.clf()
