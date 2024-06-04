from JetNet import JetNet
from JetNetData import JetNetData
from JetNet_utils import PlotLoss, PlotPrediction


def main():
    features = ['genbjet1_px', 'genbjet1_py', 'genbjet1_pz', 'genbjet1_E', 'genbjet2_px', 'genbjet2_py', 'genbjet2_pz', 'genbjet2_E']
    labels = ['H_WW_px', 'H_WW_py', 'H_WW_pz', 'H_WW_E', 'X_mass']
    input_files = ["../JetNetTrain_M-500.root"]

    data = JetNetData(features, labels)
    data.ReadFiles(input_files)
    data.Shuffle()
    if len(input_files) == 1:
        data.TrainTestSplit()

    net = JetNet(features, labels)
    net.ConfigureModel(data.train_features.shape)
    history = net.Fit(data.train_features, data.train_labels)
    net.SaveModel("./models/")
    PlotLoss(history)
    if len(input_files) == 1:
        pred = net.Predict(data.test_features)
        PlotPrediction(pred, data.test_labels)

    
if __name__ == '__main__':
    main()