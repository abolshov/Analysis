from JetNet import JetNet
from JetNetData import JetNetData
from JetNet_utils import PlotPrediction

import numpy as np

def main():
    features = ['genbjet1_px', 'genbjet1_py', 'genbjet1_pz', 'genbjet1_E', 'genbjet2_px', 'genbjet2_py', 'genbjet2_pz', 'genbjet2_E']
    labels = ['H_WW_px', 'H_WW_py', 'H_WW_pz', 'H_WW_E', 'X_mass']
    input_files = ["../JetNetTrain_M-1000.root"]

    data = JetNetData(features, labels)
    data.ReadFiles(input_files)

    net = JetNet(features, labels)
    net.LoadModel("./models/JetNet_v1.keras")

    pred = net.Predict(data.data[features])
    PlotPrediction(pred, data.data)

    
if __name__ == '__main__':
    main()