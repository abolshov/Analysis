import tensorflow as tf

from typing import Callable, List
ActivationFunction = Callable[[tf.Tensor], tf.Tensor]

import numpy as np
from numpy.typing import NDArray

@tf.keras.utils.register_keras_serializable('DenseAutoencoder')
class DenseAutoencoder(tf.keras.models.Model):
    def __init__(self, 
                 *,
                 encoder_dims: List[int],
                 decoder_dims: List[int],
                 hidden_activation: str | ActivationFunction,
                 output_activation: str | ActivationFunction,
                 mean: NDArray[np.floating] | tf.Tensor = None,
                 variance: NDArray[np.floating] | tf.Tensor = None,
                 noise_std: float = 0,
                 **kwargs) -> None:
        
        super(DenseAutoencoder, self).__init__(**kwargs)
        self.encoder_dims = encoder_dims
        self.decoder_dims = decoder_dims
        self.hidden_activation = hidden_activation
        self.output_activation = output_activation
        self.mean = mean
        self.variance = variance
        self.noise_std = noise_std

        self.gaussian_noise = tf.keras.layers.GaussianNoise(noise_std)


        

    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded
