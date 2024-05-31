from JetNet import JetNet
from JetNet_utils import PlotLoss, PlotPrediction


def main():
    net = JetNet()
    input_file = "../JetNetTrain_M-500.root"
    net.ReadFile(input_file)
    history = net.Fit()
    PlotLoss(history)
    pred = net.Predict()
    PlotPrediction(pred, net.test_labels)

    
if __name__ == '__main__':
    main()