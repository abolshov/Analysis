import tensorflow as tf
import pandas as pd
import numpy as np
import os
import yaml


from Dataloader import Dataloader
from NetUtils import *
from LossUtils import *
from LayerUtils import *


from MiscUtils import *
from PlotUtils import PlotHist
from ModelUtils import LoadModel


def main():
    # compare multiple models on multiple datasets (e.g. signal vs bkg)
    # 1. same model, different datasets
    # 2. same dataset, different models 

    # load datasets for model comparison
    file = 'nano_0.root'
    dataloader = Dataloader('dataloader_config.yaml')
    dataloader.Load(file)
    X, input_names, y, target_names = dataloader.Get()

    # load models
    model_configs = ['model_cfg.yaml']
    models = [LoadModel(cfg) for cfg in model_configs]


if __name__ == '__main__':
    main()