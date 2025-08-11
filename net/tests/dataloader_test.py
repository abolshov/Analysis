from Dataloader import Dataloader


def main():
    file = 'nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)
    X, input_names, y, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1)
    print(X.shape)
    print(y.shape)

    obj_data = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                     as_df=False,
                                     selector=lambda df: df['event'] % 2 == 1)
    num_features = 0
    for obj, (names, arr) in obj_data.items():
        print(f'{obj}: {arr.shape}')
        num_features += arr.shape[1]

    print(f'num_features: {num_features}, X.shape[1]: {X.shape[1]}')

if __name__ == '__main__':
    main()