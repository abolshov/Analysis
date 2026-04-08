import tensorflow as tf 
import numpy as np

mh = tf.constant(125.0)
pi = tf.constant(np.pi)

@tf.keras.utils.register_keras_serializable('EnergyLayer')
class EnergyLayer(tf.keras.layers.Layer):
    def __init__(self, num_quantiles=1, normalize=False, means=None, scales=None, name='energy_layer', **kwargs):
        """
        assume that layout of means and scales is identical and is [px, py, pz, E]
        """
        super(EnergyLayer, self).__init__(name=name, **kwargs)
        self.num_quantiles = num_quantiles
        self.normalize = normalize
        
        self.means = means
        self.scales = scales
        if normalize: 
            if means is None or scales is None:
                raise RuntimeError(f'EnergyLayer: `means` or `scales` must not be `None` if `normalize` is `True`')

            # convert to tensor to make call happy and be able to automatically infer output shape
            self.means = tf.constant(means, dtype=tf.float32)
            self.scales = tf.constant(scales, dtype=tf.float32)

    def call(self, inputs):
        """
        assume that input shape is (num_events, num_quantiles*3)
        *3 because px, py, pz

        assume layout is [px_q16, px_q50, px_q84, ..., pz_q16, pz_q50, pz_q84]
        so that variables for the same quantile are within displacement 3 from each other

        assume that if normalize=True layer recieves (and outputs) z-scores 
        """

        # if normalize:
        #     assert self.means and self.scales "`normalize` parameter set to true but `means` and `scales` were not provided"

        energies = []
        for q in range(self.num_quantiles):
            p3 = inputs[:, q::self.num_quantiles]
            if self.normalize:
                # transform recieved input zs-scores to actual momenta
                p3 = p3*self.scales[:3] + self.means[:3]
           
            energy_q_sqr = tf.square(mh) + tf.reduce_sum(tf.square(p3), axis=1)           
            energy_q = tf.sign(energy_q_sqr)*tf.sqrt(energy_q_sqr)
            energy_q = tf.expand_dims(energy_q, axis=-1) 
            
            if self.normalize:
                # transform calculated energies to z scores
                energy_q = (energy_q - self.means[-1])/self.scales[-1]
            energies.append(energy_q)

        output = tf.concat(energies, axis=-1)
        return output

    def get_config(self):
        config = super().get_config()
        config.update({'num_quantiles': self.num_quantiles,
                       'normalize': self.normalize,
                       'means': tf.keras.utils.serialize_keras_object(self.means) if self.means is not None else None,
                       'scales': tf.keras.utils.serialize_keras_object(self.scales) if self.scales is not None else None})
        return config

    @classmethod
    def from_config(cls, config):
        """
        This method must be implemented if custom object takes any args other than builtin python types
        """
        means = config.pop('means')
        scales = config.pop('scales')
        if means:
            means = means['config']['value']
        if scales:
            scales = scales['config']['value']
        return cls(means=means, scales=scales, **config)
    

class QuantileOrderingLayer(tf.keras.layers.Layer):
    """
    https://www.kaggle.com/code/syerramilli/pfi-multi-quantile-neural-network-regression
    """
    def __init__(self, name='quantile_ordering_layer', **kwargs):
        super(QuantileOrderingLayer, self).__init__(name=name, **kwargs)
        
    def call(self, inputs):
        # Compute the left and right parts of the quantiles
        first_quants = inputs[..., 0:1]
        rem_quants = first_quants + tf.math.cumsum(tf.nn.softplus(inputs[..., 1:]), axis=-1)

        # Concatenate the parts along the last dimension
        out = tf.concat([first_quants, rem_quants], axis=-1)
        return out
    
@tf.keras.utils.register_keras_serializable('ResidualBlock')
class ResidualBlock(tf.keras.layers.Layer):
    def __init__(self, 
                 *,
                 units,
                 activation, 
                 **kwargs):
        
        super().__init__(**kwargs)
        
        self.add = tf.keras.layers.Add()
        self.activation = tf.keras.activations.get(activation)
        
        self.dense1 = None
        self.bn1 = None
        self.dense2 = None
        self.bn2 = None
        self.units = units

    def build(self, input_shape):
        input_dim = input_shape[-1]
        
        self.bn1 = tf.keras.layers.BatchNormalization()
        self.fc1 = tf.keras.layers.Dense(self.units)
        self.bn2 = tf.keras.layers.BatchNormalization()
        self.fc2 = tf.keras.layers.Dense(self.units)
        self.relu = tf.keras.layers.ReLU()
        
        # Projection for dimension mismatch
        if input_dim != self.units:
            self.projection = tf.keras.layers.Dense(self.units, use_bias=False)
        else:
            self.projection = None

        return super().build(input_shape)

    def call(self, 
             inputs, 
             *, 
             training=False):
        
        if self.projection is not None:
            residual = self.projection(inputs)
        else:
            residual = inputs
            
        x = self.fc1(self.relu(self.bn1(inputs, training=training)))
        x = self.fc2(self.relu(self.bn2(x, training=training)))
        x = self.add([x, residual])
        x = self.activation(x)
        return x
    
    def get_config(self):
        config = super().get_config()
        config.update({
            'units': self.units,
            'activation': tf.keras.activations.serialize(self.activation)
        })
        return config