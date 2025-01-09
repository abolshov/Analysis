import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import matplotlib.pyplot as plt


from DataWrapper import DataWrapper


target_names = {"genHbb_E": "E(H->bb)",
                "genHbb_px": "px(H->bb)",
                "genHbb_py": "py(H->bb)",
                "genHbb_pz": "pz(H->bb)",
                "genHVV_E": "E(H->VV)",
                "genHVV_px": "px(H->VV)",
                "genHVV_py": "py(H->VV)",
                "genHVV_pz": "pz(H->VV)" }


def MakeBarPlot(x, y, target_name):
    plt.barh(y, x)
    plt.xlabel("corellation coefficient") 
    plt.title(f"Correlations of {target_names[target_name]} with input features")
    plt.savefig(f"correlation_plots/sl/corr_{target_name}.pdf", bbox_inches='tight')
    plt.clf()


dw = DataWrapper()
# dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2L2Nu_M-800/nano_0.root")
dw.ReadFile("/Users/artembolshov/Desktop/CMS/Di-Higgs/data/GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-800/nano_0.root")
df = dw.data
df = df.drop(columns=["event"])

labels_list = dw.labels 
features_list = dw.feature_list

correlations_all = df.corr()
target_feature_corr = correlations_all.filter(labels_list).drop(labels_list)

for target_name in labels_list:
    sorted_corr = target_feature_corr[target_name].sort_values(ascending=True, key=lambda x: abs(x))
    sorted_corr = sorted_corr[-10:]
    MakeBarPlot(sorted_corr.values, list(sorted_corr.index), target_name)
