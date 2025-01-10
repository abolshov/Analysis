import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt
from DataWrapper import DataWrapper
import math 
import tensorflow as tf
import seaborn as sns


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

def PlotCompare(hbb_value, hvv_value, feature, target):
    min_val = min(np.min(hbb_value), np.min(hvv_value))
    max_val = max(np.max(hbb_value), np.max(hvv_value))
    bins = np.linspace(min_val, max_val, 100)
    plt.hist(hbb_value, bins=bins, label="H->bb", alpha=0.8)
    plt.hist(hvv_value, bins=bins, label="H->VV", alpha=0.8)
    plt.title(f'Gradient of {target} wrt {feature}')
    plt.ylabel('Count')
    plt.xlabel(f'{feature}')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig(f"gradient_plots/cmp_grad_{target}_{feature}.pdf", bbox_inches='tight')
    plt.clf()


class PredictiveModel():
    """
    implements a vector valued function
                      f1(x1, ..., xm)
    f(x1, ..., xm) =        ...
                      fn(x1, ..., xm)
    """
    def __init__(self, model_name, input_names, target_names):
        self.input_names = input_names
        self.target_names = target_names
        self.model = tf.keras.models.load_model(model_name)


    def Predict(self, inputs):
        return self.model.predict(inputs, verbose=0)


    def Gradient(self, inputs):
        epsilon = 10e-5
        # epsilon = 0.0001

        n_examples, n_features = inputs.shape
        n_targets = len(self.target_names)

        grad = np.zeros((n_features, n_examples, n_targets))
        
        # loop over features and calculate partial derivatives of model wrt to this feature for ALL examples
        for i in range(n_features):
            print(f"Evaluating gradient wrt {self.input_names[i]}")
            
            inputs[:, i] -= epsilon
            fm1 = self.model.predict(inputs, verbose=0)
            inputs[:, i] -= epsilon
            fm2 = self.model.predict(inputs, verbose=0)

            inputs[:, i] += 3*epsilon
            fp1 = self.model.predict(inputs, verbose=0)
            inputs[:, i] += epsilon
            fp2 = self.model.predict(inputs, verbose=0)
            
            # matrix of gradients of each target wrt ith feature for all examples
            # shape: (n_examples, n_targets)
            res = (fm2 - 8*fm1 + 8*fp1 - fp2)/(12*epsilon)

            grad[i, :] = res
            inputs[:, i] -= 2*epsilon

        return grad


dw = DataWrapper()
dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
dw.FormTestSet()

feature_names = list(dw.test_features.columns)
target_names = list(dw.test_labels.columns)

model_path = "augmented_full/model_hh_dl.keras"
pm = PredictiveModel(model_path, feature_names, target_names)

# input_features = dw.test_features
# val = -0.1
# val = input_features.iloc[0, 0] - 0.5
# n_steps = 100
# step = 0.05
# for i in range(n_steps):
#     val += i*step
#     input_features.iloc[0, 0] = val
#     prediction = pm.Predict(input_features)
#     print(f"pz(H->VV)={prediction[0, 6]}, E(H->VV)={prediction[0, 7]}")
#     val += step

inputs = dw.test_features.values
grad = pm.Gradient(inputs)

HVV_pz_grad = pd.DataFrame(grad[:, :, 6].T, columns=feature_names)
HVV_E_grad = pd.DataFrame(grad[:, :, 7].T, columns=feature_names)

Hbb_pz_grad = pd.DataFrame(grad[:, :, 2].T, columns=feature_names)
Hbb_E_grad = pd.DataFrame(grad[:, :, 3].T, columns=feature_names)

for col in HVV_pz_grad.columns:
    hbb_val = Hbb_pz_grad[col]
    hvv_val = HVV_pz_grad[col]
    PlotCompare(hbb_val, hvv_val, col, "pz")

    hbb_val = Hbb_E_grad[col]
    hvv_val = HVV_E_grad[col]
    PlotCompare(hbb_val, hvv_val, col, "E")
