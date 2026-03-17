import tensorflow as tf

from typing import Callable, List
ActivationFunction = Callable[[tf.Tensor], tf.Tensor]

import numpy as np
from numpy.typing import NDArray

@tf.keras.utils.register_keras_serializable('AnomalyDetector')
class Autoencoder(tf.keras.models.Model):
    def __init__(self, 
                 *,
                 encoder_dims: List[int],
                 decoder_dims: List[int],
                 hidden_activation: str | ActivationFunction,
                 output_activation: str | ActivationFunction,
                 mean: NDArray[np.floating] | tf.Tensor = None,
                 variance: NDArray[np.floating] | tf.Tensor = None,
                 **kwargs) -> None:
        
        super(Autoencoder, self).__init__(**kwargs)

        self.hidden_activation = tf.keras.activations.get(hidden_activation)
        self.output_activation = tf.keras.activations.get(output_activation)

        self.encoder = tf.keras.Sequential()
        self.decoder = tf.keras.Sequential()

        use_norm = mean and variance
        if use_norm:
            self.decoder.add(tf.keras.layers.Normalization(mean=mean, variance=variance))
            
        self.encoder.add(
            [
                tf.keras.layers.Dense(dim, activation=hidden_activation) for dim in encoder_dims
            ]
        )

        self.decoder = tf.keras.Sequential(
            [
                tf.keras.layers.Dense(dim, activation=hidden_activation) for dim in decoder_dims[:-1]
            ]
        )
        self.decoder.add(tf.keras.layers.Dense(decoder_dims[-1], activation=output_activation))
        if use_norm:
            self.decoder.add(tf.keras.layers.Normalization(mean=mean, variance=variance, invert=True))

    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded
