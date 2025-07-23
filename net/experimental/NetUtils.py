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
            # name='beta',
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


