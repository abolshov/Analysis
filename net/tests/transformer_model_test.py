import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import tensorflow as tf
from experimental.Transformer import Transformer
from experimental.LossUtils import QuantileLoss

target_names = ['genHVV_px', 'genHVV_py', 'genHVV_pz', 'genHVV_E',
                'genHbb_px', 'genHbb_py', 'genHbb_pz', 'genHbb_E']

def main():
    regressor_cfg = {'num_quantiles': 3, 
                     'target_names': target_names, 
                     'use_energy_layer': True,
                     'use_quantile_ordering': False,
                     'means': [0.0 for _ in range(len(target_names))],
                     'scales': [1.0 for _ in range(len(target_names))]}
    transformer = Transformer(num_encoder_layers=6, d_model=128, num_heads=8, dff=512, regressor_cfg=regressor_cfg)

    lep_tensor = tf.random.normal(shape=(10, 8))
    jet_tensor = tf.random.normal(shape=(10, 120))
    fatjet_tensor = tf.random.normal(shape=(10, 27))
    met_tensor = tf.random.normal(shape=(10, 2))

    inputs = [lep_tensor, met_tensor, jet_tensor, fatjet_tensor]
    transformer_output = transformer(inputs)
    print(transformer.summary())

    losses = {}
    for target_name in target_names:
        losses[target_name] = QuantileLoss(quantiles=[0.16, 0.5, 0.84], 
                                           name=f'quantile_loss_{target_name}')
    print('Compiling done')

    transformer.compile(loss=losses,
                        optimizer=tf.keras.optimizers.Adam(3e-4))

    transformer.save('transformer_model.keras')
    loaded_transformer = tf.keras.models.load_model('transformer_model.keras')
    transformer_output = loaded_transformer(inputs)
    print(transformer.summary())

if __name__ == '__main__':
    main()