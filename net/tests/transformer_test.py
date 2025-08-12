import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import tensorflow as tf
from experimental.Transformer import Transformer

def main():
    transformer = Transformer(num_encoder_layers=6, d_model=128, num_heads=8, dff=512)

    lep_tensor = tf.random.normal(shape=(10, 8))
    jet_tensor = tf.random.normal(shape=(10, 120))
    fatjet_tensor = tf.random.normal(shape=(10, 27))
    met_tensor = tf.random.normal(shape=(10, 2))

    inputs = [lep_tensor, met_tensor, jet_tensor, fatjet_tensor]
    transformer_output = transformer(inputs)

    print(f'Transformer output list length: {len(transformer_output)}')
    print(f'Transformer output shape: {transformer_output[0].shape}')
    print(f'Transformer param count: {transformer.count_params()}')


if __name__ == '__main__':
    main()