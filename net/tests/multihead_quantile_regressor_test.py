import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import tensorflow as tf
from experimental.LayerUtils import MultiheadQuantileRegressor

def main():
    target_names = ['genHVV_px', 'genHVV_py', 'genHVV_pz', 'genHVV_E',
                    'genHbb_px', 'genHbb_py', 'genHbb_pz', 'genHbb_E']

    regressor = MultiheadQuantileRegressor(num_quantiles=3, 
                                           target_names=target_names, 
                                           use_energy_layer=True,
                                           use_quantile_ordering=False)

    inputs = tf.random.normal(shape=(10, 512))
    regressor_output = regressor(inputs)
    num_regression_heads = len(regressor.heads)
    print(f'Num regressor heads: {num_regression_heads}')
    print(f'Regressor output list length: {len(regressor_output)}')
    print(f'Regressor output shape: {regressor_output[0].shape}')
    print(f'Regressor param count: {regressor.count_params()}')


if __name__ == '__main__':
    main()