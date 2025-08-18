import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import tensorflow as tf
from experimental.LayerUtils import Encoder

def main():
    input_mapping = [1]*2 + [2] + [3]*10 + [4]*3

    encoder = Encoder(num_encoder_layers=6, d_model=64, num_heads=8, dff=256, input_mapping=input_mapping)

    lep_inputs = [tf.random.normal(shape=(10, 4)) for _ in range(2)]
    met_inputs = [tf.random.normal(shape=(10, 2)) for _ in range(1)]
    jet_inputs = [tf.random.normal(shape=(10, 12)) for _ in range(10)]
    fatjet_inputs = [tf.random.normal(shape=(10, 9)) for _ in range(3)]
    inputs = lep_inputs + met_inputs + jet_inputs + fatjet_inputs

    output = encoder(inputs)
    print(f'Encoder output shape: {output.shape}')
    print(f'Encoder param count: {encoder.count_params()}')
    print(encoder.compute_output_shape(inputs))

if __name__ == '__main__':
    main()