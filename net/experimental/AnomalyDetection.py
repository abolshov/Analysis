import tensorflow as tf
from typing import Callable
ActivationFunction = Callable[[tf.Tensor], tf.Tensor]

import numpy as np
from numpy.typing import NDArray

@tf.keras.utils.register_keras_serializable('AnomalyDetector')
class Autoencoder(tf.keras.models.Model):
    def __init__(self, 
                 *,
                 input_dim: int,
                 latent_dim: int,
                 hidden_activation: str | ActivationFunction,
                 output_activation: str | ActivationFunction,
                 mean: NDArray | tf.Tensor = None,
                 variance: NDArray | tf.Tensor = None,
                 **kwargs) -> None:
        
        super(Autoencoder, self).__init__(**kwargs)

        self.hidden_activation = tf.keras.activations.get(hidden_activation)
        self.output_activation = tf.keras.activations.get(output_activation)

        self.norm = None
        if mean and variance:
            self.norm = tf.keras.layers.Normalization(mean=mean, variance=variance)
            # not sure reverse_norm is needed
            # self.reverse_norm = tf.keras.layers.Normalization(mean=mean, variance=variance, invert=True)

        self.encoder = None
        self.decoder = None

        # I want to build model more dynmically
        # self.encoder = tf.keras.Sequential(
        #     [
        #         tf.keras.layers.Dense(64, activation=hidden_activation),
        #         tf.keras.layers.Dense(32, activation=hidden_activation),
        #         tf.keras.layers.Dense(16, activation=hidden_activation),
        #         tf.keras.layers.Dense(8, activation=hidden_activation)
        #     ]
        # )

        # self.decoder = tf.keras.Sequential(
        #     [
        #         tf.keras.layers.Dense(16, activation=hidden_activation),
        #         tf.keras.layers.Dense(32, activation=hidden_activation),
        #         tf.keras.layers.Dense(64, activation=hidden_activation),
        #         tf.keras.layers.Dense(input_dim, activation=output_activation)
        #     ]
        # )

    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded
