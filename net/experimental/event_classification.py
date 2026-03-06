import tensorflow as tf

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print(f"Enabled memory growth for {len(gpus)} GPU(s).")
    except RuntimeError as e:
        # Memory growth must be set before GPUs have been initialized
        print(e)

import uproot
import awkward as ak
import numpy as np
import os
import re

def main():
    print("Training DNN for event classification")

    weight_file_pattern = re.compile(r"nParity(\d)_Merged_weight.root")
    event_file_pattern = re.compile(r"nParity(\d)_Merged.root")

    train_file_dir = "/home/artem/Desktop/CMS/data/DNN/SL/resolved/Run3_2022/Dataset"
    event_files = {}
    weight_files = {}
    for f in os.listdir(train_file_dir):
        weight_match = weight_file_pattern.search(f)
        event_match = event_file_pattern.search(f)
        abs_path = os.path.abspath(os.path.join(train_file_dir, f))
        if weight_match:
            parity = int(weight_match.group(1))
            weight_files[parity] = abs_path
        elif event_match:
            parity = int(event_match.group(1))
            event_files[parity] = abs_path
        else:
            continue

    parity_file_map = { p: (event_files[p], weight_files[p]) for p in event_files.keys()}
    for p, (event_file, weight_file) in parity_file_map.items():
        print(f"Parity {p}: event_file={event_file}, weight_file={weight_file}")

if __name__ == "__main__":
    main()