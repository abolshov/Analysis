import tensorflow as tf
import pandas as pd
import numpy as np
import uproot 
import awkward as ak 
import os
import yaml
import matplotlib.pyplot as plt

from Dataloader import Dataloader
from LossUtils import QuantileLoss
from Transformer import Transformer

from MiscUtils import ground_truth_map, Scheduler, pretty_vars, objects
from PlotUtils import PlotMetric, PlotHist, PlotCompare2D, PlotCovarMtrx

def main():
    params = {}

    # files = []
    # input_files = 'files_Run3_2022.txt'
    # with open(input_files, 'r') as file_cfg:
    #     files = [line[:-1] for line in file_cfg.readlines()]

    dataloader = Dataloader('dataloader_config.yaml')
    # dataloader.Load(files)
    dataloader.Load('../nano_0.root')
    inputs = dataloader.GetObjData(objects=['lep1', 'lep2', 'met', 'centralJet', 'SelectedFatJet'], 
                                   as_df=False,
                                   selector=lambda df: df['event'] % 2 == 0)

    print('Inputs:')
    for obj, (names, arr) in inputs.items():
        print(f'\t{obj}: {arr.shape}')

    targets = dataloader.GetObjData(objects=['genHbb', 'genHVV'], 
                                    as_df=False,
                                    selector=lambda df: df['event'] % 2 == 0)

    print('Targets:')
    for obj, (names, arr) in targets.items():
        print(f'\t{obj}: {arr.shape}')


if __name__ == '__main__':
    main()