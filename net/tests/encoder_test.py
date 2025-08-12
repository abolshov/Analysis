import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import tensorflow as tf
from experimental.LayerUtils import Encoder

def main():
    encoder = Encoder(num_encoder_layers=6, d_model=128, num_heads=8, dff=512)
    lep_tensor = tf.random.normal(shape=(10, 8))
    met_tensor = tf.random.normal(shape=(10, 2))
    jet_tensor = tf.random.normal(shape=(10, 120))
    fatjet_tensor = tf.random.normal(shape=(10, 27))
    encoder_output = encoder([lep_tensor, met_tensor, jet_tensor, fatjet_tensor])
    print(f'Encoder output shape: {encoder_output.shape}')
    print(f'Encoder param count: {encoder.count_params()}')

if __name__ == '__main__':
    main()