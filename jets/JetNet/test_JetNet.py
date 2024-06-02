from JetNet import JetNet
from JetNet_utils import PlotPrediction


def main():
    net = JetNet()
    input_file = "../JetNetTrain_M-700.root"
    net.LoadModel("./models/JetNet_v1.keras")
    net.ReadFile(input_file)
    pred = net.Predict()
    PlotPrediction(pred, net.test_labels)

    
if __name__ == '__main__':
    main()