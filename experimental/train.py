import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from DataWrapper import DataWrapper


def Scheduler(epoch, lr):
    if epoch < 20:
        return lr
    else:
        if epoch % 2 == 0:
            return 0.9*lr
        return lr


# def Scheduler(epoch, lr):
#     if epoch % 2 == 0:
#         return 0.95*lr
#     return lr

mh = 125.0

def CustomLoss(y_true, y_pred):
    # base_loss = tf.keras.losses.log_cosh(y_true, y_pred)
    base_loss = tf.keras.losses.mse(y_true, y_pred)

    hbb_p3 = y_pred[:, 0:3]
    hbb_E = y_pred[:, 3]
    mh_bb_sqr = hbb_E*hbb_E - tf.reduce_sum(hbb_p3**2, axis=1, keepdims=True)

    hww_p3 = y_pred[:, 4:7]
    hww_E = y_pred[:, 7]
    mh_ww_sqr = hww_E*hww_E - tf.reduce_sum(hww_p3**2, axis=1, keepdims=True)

    mh_bb = tf.sign(mh_bb_sqr)*tf.sqrt(tf.abs(mh_bb_sqr))
    mh_ww = tf.sign(mh_ww_sqr)*tf.sqrt(tf.abs(mh_ww_sqr))

    # custom = 0.5*tf.keras.losses.log_cosh(mh_bb, mh) + 0.5*tf.keras.losses.log_cosh(mh_ww, mh)
    custom_loss = tf.reduce_mean(0.5*(mh_bb - mh)**2 + 0.5*(mh_ww - mh)**2)
    return base_loss + 0.1*custom_loss
    # return custom_loss


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

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ['TF_DETERMINISTIC_OPS'] = '1'
tf.random.set_seed(42)

model_name = "model_hh_custloss"
cpkt_path = f"{model_name}_cpkt.keras"

schedule_callback = tf.keras.callbacks.LearningRateScheduler(Scheduler)
checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(cpkt_path, 
                                                         monitor='val_loss', 
                                                         verbose=0, 
                                                         save_best_only=False, 
                                                         mode='min')

model = None
resume = False
if os.path.exists(cpkt_path):
    resume = True
    model = tf.keras.models.load_model(cpkt_path) 
else:
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    # model.add(tf.keras.layers.Dropout(0.2))
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    # model.add(tf.keras.layers.Dropout(0.2))
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    # model.add(tf.keras.layers.Dropout(0.2))
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    # model.add(tf.keras.layers.Dropout(0.2))
    model.add(tf.keras.layers.Dense(256, activation='relu', kernel_initializer='random_normal', bias_initializer='random_normal'))
    # model.add(tf.keras.layers.BatchNormalization())
    # model.add(tf.keras.layers.Dropout(0.2))
    model.add(tf.keras.layers.Dense(8))

if resume:
    model.compile(loss='log_cosh', 
                  optimizer=tf.keras.optimizers.Adam(1.29e-4))
else:
    # model.compile(loss='log_cosh', 
    #               optimizer=tf.keras.optimizers.Adam(3e-4))
    model.compile(loss=CustomLoss, 
                  optimizer=tf.keras.optimizers.Adam(3e-4))

model.build(dw.train_features.shape)

history = None
if resume:
    history = model.fit(dw.train_features,
                        dw.train_labels,
                        validation_split=0.2,
                        verbose=1,
                        batch_size=128,
                        epochs=30)
else:
    history = model.fit(dw.train_features,
                        dw.train_labels,
                        validation_split=0.2,
                        verbose=1,
                        batch_size=256,
                        epochs=80,
                        callbacks=[schedule_callback, checkpoint_callback])

model.save(f"{model_name}.keras")
PlotMetric(history, model_name, "loss")