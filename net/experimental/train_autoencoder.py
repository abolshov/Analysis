import tensorflow as tf

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print(f"Enabled memory growth for {len(gpus)} GPU(s).")
    except RuntimeError as e:
        print(e)

import uproot
import pathlib
import awkward as ak
import numpy as np
import os
import re

from MiscUtils import MemoryMonitor

def main():
    print(f"Training autoencoder")
    mm = MemoryMonitor()
    mm.print_memory_usage()

if __name__ == "__main__":
    main()