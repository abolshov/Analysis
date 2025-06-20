import tensorflow as tf 

mh = tf.constant(125.0)

class TrainableSiLU(tf.keras.layers.Layer):
    def __init__(self, units=1, **kwargs):
        super(TrainableSiLU, self).__init__(**kwargs)
        self.units = units

    def build(self, input_shape):
        self.beta = self.add_weight(
            name='beta',
            shape=(self.units,),
            initializer='ones',  # Initialize beta to 1.0
            trainable=True,
            dtype="float32"
        )
        super(TrainableSiLU, self).build(input_shape)

    def call(self, inputs):
        return inputs * tf.keras.activations.sigmoid(self.beta * inputs)

    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super(TrainableSiLU, self).get_config()
        config.update({"beta": self.beta.numpy()})
        config.update({"units": self.units})
        return config

    @classmethod
    def from_config(cls, config):
        # Reconstruct the layer from the saved configuration
        return cls(**config)


class CombinedLoss(tf.keras.losses.Loss):
    def __init__(self, strength=0.01, name="combined_loss", **kwargs):
        super().__init__(name=name)
        self.strength = strength
        self.logcosh_loss = tf.keras.losses.LogCosh() 
        self.mse_loss = tf.keras.losses.MeanSquaredError() 

    def call(self, y_true, y_pred):
        logcosh = self.logcosh_loss(y_true, y_pred)
        mbb2 = tf.square(y_pred[:, 0]) - tf.reduce_sum(tf.square(y_pred[:, 1:4]), axis=1, keepdims=True)
        mww2 = tf.square(y_pred[:, 4]) - tf.reduce_sum(tf.square(y_pred[:, 5:8]), axis=1, keepdims=True)
        mbb = tf.sign(mbb2)*tf.sqrt(tf.abs(mbb2))
        mww = tf.sign(mww2)*tf.sqrt(tf.abs(mww2))
        # penalty = tf.sqrt(0.5*self.mse_loss(mh, mbb)) + 0.5*tf.sqrt(self.mse_loss(mh, mww))
        penalty = 0.5*self.mse_loss(mh, mbb) + 0.5*self.mse_loss(mh, mww)
        total_loss = logcosh + self.strength*penalty
        return total_loss 

    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super().get_config()
        config.update({"strength": self.strength})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)