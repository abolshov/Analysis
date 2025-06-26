import tensorflow as tf 

mh = tf.constant(125.0)

# decorators seem to be unavailable with my version of keras
# @tf.keras.saving.register_keras_serializable(package="TrainableSiLU")
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
        config.update({"name": "beta", 
                       "units": self.units})
        return config

    @classmethod
    def from_config(cls, config):
        # Reconstruct the layer from the saved configuration
        return cls(**config)

# @tf.keras.saving.register_keras_serializable(package="CombinedLoss")
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
        penalty = 0.5*self.mse_loss(mh, mbb) + 0.5*self.mse_loss(mh, mww)
        total_loss = logcosh + self.strength*penalty
        # print(f"penalty_frac={penalty.numpy()[0]/total_loss.numpy()[0]:.2f}, logcosh_frac={logcosh.numpy()[0]/total_loss.numpy()[0]:.2f}")
        # tf.print(self.strength*penalty/total_loss)
        # tf.print(logcosh/total_loss)
        return total_loss 

    def update_epoch(self):
        self.current_epoch.assign_add(1)

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


class MassLoss(tf.keras.losses.Loss):
    def __init__(self, strength=1.0, name="mass_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.mse_loss = tf.keras.losses.MeanSquaredError()
        self.strength = strength

    def call(self, y_true, y_pred):
        mbb2 = tf.square(y_pred[:, 3]) - tf.reduce_sum(tf.square(y_pred[:, 0:3]), axis=1, keepdims=True)
        mww2 = tf.square(y_pred[:, 7]) - tf.reduce_sum(tf.square(y_pred[:, 4:7]), axis=1, keepdims=True)
        mbb = tf.sign(mbb2)*tf.sqrt(tf.abs(mbb2))
        mww = tf.sign(mww2)*tf.sqrt(tf.abs(mww2))
        loss = 0.5*self.mse_loss(mh, mbb) + 0.5*self.mse_loss(mh, mww)
        return self.strength*loss 

    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super().get_config()
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)


class EpochLossUpdater(tf.keras.callbacks.Callback):
    def __init__(self, initial_value, update_rate):
        super().__init__()
        self.loss_parameter = initial_value
        self.update_rate = update_rate

    def on_epoch_end(self, epoch, logs=None):
        # Update the loss_param at the end of each epoch
        new_value = self.loss_parameter*self.update_rate  
        self.loss_parameter.assign(new_value)
        print(f"\nEpoch {epoch + 1}: Loss parameter updated to {self.loss_parameter.numpy():.6e}")


