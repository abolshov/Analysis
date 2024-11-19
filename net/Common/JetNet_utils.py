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


# output: (n_events, 11) - produced by the net
# in output:
# first three: px, py, pz of H->bb
# next three: px, py, pz of H->WW
# next four: E, px, py, pz of W1
# W1 <-> genV1, W2 <-> genV2
@tf.function
def GetPredMasses(output):
    H_mass = 125.0

    H_bb_p3 = output[:, 0:3]
    H_bb_E = tf.sqrt(H_mass*H_mass + tf.reduce_sum(H_bb_p3**2, axis=1, keepdims=True))

    H_WW_p3 = output[:, 3:6]
    H_WW_E = tf.sqrt(H_mass*H_mass + tf.reduce_sum(H_WW_p3**2, axis=1, keepdims=True))

    W1_E = output[:, 6]
    W1_p3 = output[:, 7:10]

    W2_E = H_WW_E - W1_E
    W2_p3 = H_WW_p3 - W1_p3

    mW1_sqr = W1_E**2 - tf.reduce_sum(W1_p3**2, axis=1, keepdims=True)
    mW2_sqr = W2_E**2 - tf.reduce_sum(W2_p3**2, axis=1, keepdims=True)
    mX_sqr = (H_bb_E + H_WW_E)**2 - tf.reduce_sum((H_WW_p3 + H_bb_p3)**2, axis=1, keepdims=True)

    mW1 = tf.sign(mW1_sqr)*tf.sqrt(tf.abs(mW1_sqr))
    mW2 = tf.sign(mW2_sqr)*tf.sqrt(tf.abs(mW2_sqr))
    mX = tf.sign(mX_sqr)*tf.sqrt(tf.abs(mX_sqr))

    return tf.concat([mX, mW1, mW2], axis=1)


# target: (n_events, 1) - true mass of X
# output: (n_events, 6) - produced by the net
# in output:
# first three: px, py, pz of H->bb
# last three: px, py, pz of H->WW
@tf.function
def MXLossFunc(target, output):
    pred = GetMXPred(output)
    return (target - pred)**2/tf.cast(len(pred), tf.float32)


# target: (n_events, 3) - true mass of X, W1, W2
# output: (n_events, 11) - produced by the net
# in output:
# first three: px, py, pz of H->bb
# next three: px, py, pz of H->WW
# next four: E, px, py, pz of W1
@tf.function
def ThreeMassLossFunc(target, output):
    pred = GetPredMasses(output)
    mX_true = target[:, 0]
    mW1_true = target[:, 1]
    mW2_true = target[:, 2] 
    mX_pred = pred[:, 0]
    mW1_pred = pred[:, 1]
    mW2_pred = pred[:, 2]

    mX_loss = ((mX_true - mX_pred)/mX_true)**2
    mW1_loss = ((mW1_true - mW1_pred)/tf.maximum(mW1_true, 20.0))**2
    mW2_loss = ((mW2_true - mW2_pred)/tf.maximum(mW2_true, 20.0))**2
    
    return 1000*(mX_loss + mW1_loss + mW2_loss)/tf.cast(3*len(mX_pred), tf.float32)


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
