import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from Dataset import Dataset
from NetUtils import TrainableSiLU, EpochLossUpdater
from LossUtils import QuantileLoss


ground_truth_map = {"genHbb_E": "E(H->bb)",
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


def PlotCompare2D(target, output, quantity, quantile, plotting_dir=None):
    plt.grid(False)
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist2d(target, output, bins=bins)
    var = quantity.split('_')[1]
    plt.title(f'{var} comparison')
    # quantity = '_'.join(quantity.split('_')[:-1])
    plt.ylabel(f'predicted {ground_truth_map[quantity]}')
    plt.xlabel(f'true {ground_truth_map[quantity]}')
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f"cmp2d_{quantity}_{quantile}.pdf"), bbox_inches='tight')
    else:
        plt.savefig(f"cmp2d_{quantity}_{quantile}.pdf", bbox_inches='tight')
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


dataset = Dataset('dataset_config.yaml')
dataset.Load()
X_train, input_names, y_train, target_names = dataset.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'
tf.random.set_seed(42)

model_name = "model_hh_dl_train_silu_l2reg_qloss"
model_dir = model_name
os.makedirs(model_dir, exist_ok=True)

learning_rate_scheduler = tf.keras.callbacks.LearningRateScheduler(Scheduler)

input_shape = X_train.shape[1:]
num_units = 384
num_layers = 12
l2_strength = 0.01
quantiles = [0.16, 0.5, 0.84]
num_quantiles = len(quantiles)
epochs = 10
batch_size = 512

# prepare base part of the model
inputs = tf.keras.layers.Input(shape=input_shape)
x = tf.keras.layers.Dense(num_units, 
                          activation=TrainableSiLU(units=num_units), 
                          kernel_initializer='random_normal', 
                          bias_initializer='random_normal',
                          kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(inputs)
for _ in range(num_layers - 1):
    x = tf.keras.layers.Dense(num_units, 
                              activation=TrainableSiLU(units=num_units), 
                              kernel_initializer='random_normal', 
                              bias_initializer='random_normal',
                              kernel_regularizer=tf.keras.regularizers.L2(l2_strength))(x)

# prepeare output heads
num_targets = len(target_names)
losses = {}
outputs = []
ys_train = []
for idx, target_name in enumerate(target_names):
    target_array = y_train[:, idx]
    target_array = np.reshape(target_array, (-1, 1))
    losses[target_name] = QuantileLoss(quantiles=quantiles, name=f'loss_{target_name}')

    # repeat each target array num_quantiles times so that each output node has ground truth value
    target_array = np.repeat(target_array, repeats=num_quantiles, axis=-1)
    ys_train.append(target_array)

    output = tf.keras.layers.Dense(num_quantiles, name=target_name)(x)
    outputs.append(output)

model = tf.keras.Model(inputs=inputs, outputs=outputs)
model.compile(loss=losses, optimizer=tf.keras.optimizers.Adam(3e-4))

history = model.fit(X_train,
                    ys_train,
                    validation_split=0.2,
                    verbose=1,
                    batch_size=batch_size,
                    epochs=epochs,
                    callbacks=[learning_rate_scheduler])

model.save(os.path.join(model_dir, f"{model_name}.keras"))
PlotMetric(history, model_name, "loss", plotting_dir=model_dir)

# form testing dataset
def TestSelection(df, mod, parity, mass):
    return np.logical_and(df['event'] % mod == parity, df['X_mass'] == mass)

X_test, _, y_test, _ = dataset.Get(TestSelection, 2, 1, 800)
ys_pred = model.predict(X_test)

assert len(ys_pred) == len(target_names), f"mismatch between number of predicted values and ground truth values ({len(ys_pred)} vs {len(target_names)})"

# ys_pred is a list of np arrays of shape (num_events, num_quantiles)
pred_dict = {}
for i, name in enumerate(target_names):
    pred_array = ys_pred[i]
    for q_idx, quantile in enumerate(quantiles):
        q_pred_array = pred_array[:, q_idx]
        q = int(100*quantile)
        pred_dict[f'{name}_q{q}'] = q_pred_array

pred_df = pd.DataFrame.from_dict(pred_dict)
print(pred_df.head())
print(pred_df.columns)

bins = np.linspace(0, 2000, 100)
for quantile in quantiles:
    q = int(100*quantile)

    for i, col in enumerate(target_names):
        ground_truth = y_test[:, i]
        PlotCompare2D(ground_truth, pred_df[f'{col}_q{q}'], col, q, plotting_dir=model_dir)

    X_E_pred = pred_df[f"genHbb_E_q{q}"] + pred_df[f"genHVV_E_q{q}"]
    X_px_pred = pred_df[f"genHbb_px_q{q}"] + pred_df[f"genHVV_px_q{q}"]
    X_py_pred = pred_df[f"genHbb_py_q{q}"] + pred_df[f"genHVV_py_q{q}"]
    X_pz_pred = pred_df[f"genHbb_pz_q{q}"] + pred_df[f"genHVV_pz_q{q}"]

    X_mass_pred_sqr = X_E_pred**2 - X_px_pred**2 - X_py_pred**2 - X_pz_pred**2
    X_mass_sqr_pos = X_mass_pred_sqr[X_mass_pred_sqr > 0.0]
    print(f"quantile {q}: positive mass fraction: {len(X_mass_sqr_pos)}/{len(X_mass_pred_sqr)}={len(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])/len(X_mass_pred_sqr):.3f}")

    X_mass_pred = np.sqrt(X_mass_pred_sqr[X_mass_pred_sqr > 0.0])
    plt.hist(X_mass_pred, bins=bins)
    plt.title('Predicted X->HH mass')
    plt.ylabel('Count')
    plt.xlabel(f'mass, [GeV]')
    plt.figtext(0.75, 0.8, f"peak: {np.median(X_mass_pred):.2f}")
    plt.figtext(0.75, 0.75, f"width: {PredWidth(X_mass_pred):.2f}")
    plt.grid(True)
    plt.savefig(os.path.join(model_dir, f"pred_mass_q{q}.pdf"), bbox_inches='tight')
    plt.clf()
