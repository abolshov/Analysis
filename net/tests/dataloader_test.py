import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from experimental.Dataloader import Dataloader


def main():
    channel = 'DL'
    masspoint = 800
    suffix = '2B2Vto2B2JLNu' if channel == 'SL' else '2B2Vto2B2L2Nu'
    file = f'../train_data/Run3_2022/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    dataloader = Dataloader('../experimental/dataloader_config.yaml')
    dataloader.Load(file)
    X, input_names, y, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1)
    print(X.shape)
    print(y.shape)

    obj_data = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                     as_df=False,
                                     unravel=True,
                                     selector=lambda df: df['event'] % 2 == 1)
    num_features = 0
    for obj, (names, arr) in obj_data.items():
        print(f'{obj}: {arr.shape}')
        num_features += arr.shape[1]

    print(f'num_features: {num_features}, X.shape[1]: {X.shape[1]}')

if __name__ == '__main__':
    main()