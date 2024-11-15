import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# output: (n_events, 6) - produced by the net
# first three: px, py, pz of H->bb
# last three: px, py, pz of H->WW
# target: (n_events, 1) - true mass of X
@tf.function
def GetMXPred(output):
    H_mass = 125.0

    H_bb_p3 = output[:, 0:3]
    H_bb_E = tf.sqrt(H_mass*H_mass + tf.reduce_sum(H_bb_p3**2, axis=1, keepdims=True))

    H_WW_p3 = output[:, 3:6]
    H_WW_E = tf.sqrt(H_mass*H_mass + tf.reduce_sum(H_WW_p3**2, axis=1, keepdims=True))

    pred_mass = tf.sqrt((H_bb_E + H_WW_E)**2 - tf.reduce_sum((H_WW_p3 + H_bb_p3)**2, axis=1, keepdims=True))
    return pred_mass


# target: (n_events, 1) - true mass of X
# output: (n_events, 6) - produced by the net
# first three: px, py, pz of H->bb
# last three: px, py, pz of H->WW
@tf.function
def MXLossFunc(target, output):
    pred = GetMXPred(output)
    return (target - pred)**2/tf.cast(len(pred), tf.float32)


def PlotLoss(history, model):
    plt.plot(history.history['loss'], label='train_loss')
    plt.plot(history.history['val_loss'], label='val_loss')
    plt.xlabel('Epoch')
    plt.ylabel('Error')
    plt.title("Loss")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"Plots/loss_{model}.pdf", bbox_inches='tight')
    plt.clf()


def PlotPrediction(label_df, predicted_df, model):
    masspoints = [int(val) for val in label_df['X_mass'].unique()]
    X_mass_true = np.array(label_df['X_mass'])
    X_mass_pred = np.array(predicted_df['X_mass_pred'])
    mass_df = pd.DataFrame({"X_mass_true": X_mass_true, "X_mass_pred": X_mass_pred})

    for mp in masspoints:
        df = mass_df[mass_df['X_mass_true'] == mp]

        width = PredWidth(df['X_mass_pred'])
        peak = PredPeak(df['X_mass_pred'])

        bins = np.linspace(0, 2000, 100)
        plt.hist(df['X_mass_pred'], bins=bins)
        plt.title(f'{model} prediction')
        plt.xlabel('X mass [GeV]')
        plt.ylabel('Count')
        plt.figtext(0.75, 0.8, f"peak: {peak:.2f}")
        plt.figtext(0.75, 0.75, f"width: {width:.2f}")
        plt.grid(True)
        plt.savefig(f"Plots/X_mass_{mp}_GeV_{model}.pdf", bbox_inches='tight')
        plt.clf()


def PredWidth(pred_mass):
    q_84 = np.quantile(pred_mass, 0.84)
    q_16 = np.quantile(pred_mass, 0.16)
    width = q_84 - q_16
    return width 


def PredPeak(pred_mass):
    counts = np.bincount(pred_mass)
    peak = np.argmax(counts)
    return peak 