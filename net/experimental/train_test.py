import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from DataWrapper import DataWrapper
from NetUtils import TrainableSiLU, CombinedLoss, EpochLossUpdater


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


def PlotCompare2D(target, output, quantity):
    plt.grid(False)
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist2d(target, output, bins=bins)
    var = quantity.split('_')[-1]
    plt.title(f'{var} comparison')
    plt.ylabel(f'predicted {target_names[quantity]}')
    plt.xlabel(f'true {target_names[quantity]}')
    plt.savefig(f"cmp2d_{quantity}.pdf", bbox_inches='tight')
    plt.clf()


def Scheduler(epoch, lr):
    if epoch < 20:
        return lr
    else:
        if epoch % 2 == 0:
            return 0.9*lr
        return lr


def PlotMetric(history, model, metric):
    plt.plot(history.history[metric], label=f'train_{metric}')
    plt.plot(history.history[f'val_{metric}'], label=f'val_{metric}')
    plt.title(f'{model} {metric}')
    plt.ylabel(metric)
    plt.xlabel('Epoch')
    plt.legend(loc='upper right')
    plt.grid(True)
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

model_name = "model_hh_dl_all_aug_trainable_silu_dynamic_comb_loss"

schedule_callback = tf.keras.callbacks.LearningRateScheduler(Scheduler)
early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor='val_loss',
    patience=3,
    min_delta=10e-3
)

penalty_strength = tf.Variable(1e-6, dtype=tf.float32)
update_rate = 1.0001
combined_loss = CombinedLoss(penalty_strength)
loss_updater = EpochLossUpdater(combined_loss.strength, update_rate)

input_shape = dw.train_features.shape

model = tf.keras.Sequential()
# model.add(tf.keras.layers.Input(shape=input_shape[1:], dtype="float32"))
# model.add(tf.keras.layers.Input(type_spec=tf.TensorSpec(shape=input_shape, dtype=tf.float32)))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
# model.add(tf.keras.layers.Dense(256, activation='silu', kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(256, activation=TrainableSiLU(units=256), kernel_initializer='random_normal', bias_initializer='random_normal'))
model.add(tf.keras.layers.Dense(8))

# model.compile(loss='logcosh', 
#               optimizer=tf.keras.optimizers.Adam(3e-4))
model.compile(loss=CombinedLoss(penalty_strength), 
              optimizer=tf.keras.optimizers.Adam(3e-4))

model.build(dw.train_features.shape)


history = model.fit(dw.train_features,
                    dw.train_labels,
                    validation_split=0.2,
                    verbose=1,
                    batch_size=2048,
                    epochs=100,
                    callbacks=[schedule_callback, loss_updater])

model.save(f"{model_name}.keras")
PlotMetric(history, model_name, "loss")

dw.Clear()
dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
# dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-800/nano_0.root")
dw.FormTestSet()

res = model.predict(dw.test_features)
pred_dict = {label: res[:, i] for i, label in enumerate(dw.test_labels.columns)}
pred_df = pd.DataFrame.from_dict(pred_dict)

# print("before inverse transform:")
print(pred_df.describe())

# pred_df *= std
# pred_df += mean

# print("after inverse transform:")
# print(pred_df.describe())

for col in pred_df.columns:
    PlotCompare2D(dw.test_labels[col], pred_df[col], col)

bins = np.linspace(0, 2000, 100)

X_E_pred = pred_df["genHbb_E"] + pred_df["genHVV_E"]
X_px_pred = pred_df["genHbb_px"] + pred_df["genHVV_px"]
X_py_pred = pred_df["genHbb_py"] + pred_df["genHVV_py"]
X_pz_pred = pred_df["genHbb_pz"] + pred_df["genHVV_pz"]

Hbb_en = np.sqrt(125.0**2 + pred_df["genHbb_px"]**2 + pred_df["genHbb_py"]**2 + pred_df["genHbb_pz"]**2)
Hww_en = np.sqrt(125.0**2 + pred_df["genHVV_px"]**2 + pred_df["genHVV_py"]**2 + pred_df["genHVV_pz"]**2)
X_en = Hbb_en + Hww_en

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
plt.savefig(f"pred_mass.pdf", bbox_inches='tight')
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

model.save(f"{model_name}.keras")