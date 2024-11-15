from Common.JetNet import JetNet
from Common.DataWrapper import DataWrapper
from Common.JetNet_utils import PlotPrediction
import yaml 
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(prog='train_net', description='Trains Neural Net Estimator')
    parser.add_argument('config_file', type=str, help="File with neural net configuration")
    parser.add_argument('files', type=str, help="File with list of input files separated by newline character")
    parser.add_argument('model_path', type=str, help="Path where to save trained model")

    args = parser.parse_args()
    config = args.config_file
    path_to_model = args.model_path
    test_files = []
    with open(args.files, 'r') as file:
        test_files = [line[:-1] for line in file.readlines()]

    if not test_files:
        raise RuntimeError(f"file {args.files} contained empty list of input files")
    
    with open(config, 'r') as stream:
        cfg = yaml.safe_load(stream)

        dw = DataWrapper(cfg)
        dw.ReadFiles(test_files)

        dw.FormTestSet()  
    
        model_name = cfg['name']

        net = JetNet(cfg)
        net.LoadModel(f"{path_to_model}/{model_name}.keras")

        pred_df = net.Predict(dw.test_features)
        print(f"Test dataset contains {pred_df.shape[0]} events")
        print(pred_df.describe())
        w = np.quantile(pred_df['X_mass_pred'], 0.84) - np.quantile(pred_df['X_mass_pred'], 0.16)
        print(f"width   {w:.6f}")
        PlotPrediction(dw.test_labels, pred_df, model_name)
        np.savetxt(f"Output/mass_{model_name}.txt", pred_df.values, fmt='%10.5f')
        np.savetxt(f"Output/eventIds_{model_name}.txt", dw.test_events, fmt='%10.1i')
    
    
if __name__ == '__main__':
    main()
