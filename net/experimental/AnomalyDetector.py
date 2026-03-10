import tensorflow as tf

@tf.keras.utils.register_keras_serializable('AnomalyDetector')
class AnomalyDetector(tf.keras.models.Model):
    def __init__(self, **kwargs):
        super(AnomalyDetector, self).__init__(**kwargs)

        self.encoder = tf.keras.Sequential(
            [
                tf.keras.layers.Dense(64, activation="relu"),
                tf.keras.layers.Dense(32, activation="relu"),
                tf.keras.layers.Dense(16, activation="relu"),
                tf.keras.layers.Dense(8, activation="relu")
            ]
        )

        self.decoder = tf.keras.Sequential(
            [
                tf.keras.layers.Dense(16, activation="relu"),
                tf.keras.layers.Dense(32, activation="relu"),
                tf.keras.layers.Dense(64, activation="relu"),
                tf.keras.layers.Dense(98, activation="sigmoid")
            ]
        )

    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded



