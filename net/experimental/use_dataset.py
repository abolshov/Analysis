from Dataset import Dataset


def main():
    dataset = Dataset("dataset_config.yaml")
    dataset.Load()
    X, input_names, y, target_names = dataset.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)
    print(X.shape) 
    print(y.shape) 

if __name__ == '__main__':
    main()