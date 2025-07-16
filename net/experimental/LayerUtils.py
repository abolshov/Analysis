import tensorflow as tf 
import numpy as np

mh = tf.constant(125.0)

class EnergyLayer(tf.keras.layers.Layer):
    def __init__(self, num_quantiles=1, normalize=False, means=None, scales=None, name='energy_layer', **kwargs):
        """
        assume that layout of means and scales is identical and is [px, py, pz, E]
        """
        super(EnergyLayer, self).__init__(name=name, **kwargs)
        self.num_quantiles = num_quantiles
        self.normalize = normalize
        
        self.means = tf.constant(np.zeros(num_quantiles), dtype=tf.float32)
        self.scales = tf.constant(np.zeros(num_quantiles), dtype=tf.float32)
        if normalize: 
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
                       'means': self.means.numpy(),
                       'scales': self.scales.numpy()})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)
    

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


class PtLayer(tf.keras.layers.Layer):
    def __init__(self, name='pt_layer', **kwargs):
        super(PtLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # inputs: pt
        p2 = inputs[:, :2]
        return tf.sqrt(tf.square(tf.reduce_sum(p2, axis=1)))


class EtaLayer(tf.keras.layers.Layer):
    def __init__(self, name='eta_layer', **kwargs):
        super(EtaLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # inputs: eta
        p3 = inputs[:, :3]
        mod = tf.sqrt(tf.square(tf.reduce_sum(p3, axis=1)))
        return tf.atanh(p3[:, -1]/mod)


class PhiLayer(tf.keras.layers.Layer):
    def __init__(self, name='phi_layer', **kwargs):
        super(PhiLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # inputs: phi
        p3 = inputs[:, :3]
        return tf.atan2(p3[:, 1], p3[:, 0])


class EnErrLayer(tf.keras.layers.Layer):
    def __init__(self, name='energy_error_layer', **kwargs):
        super(EnErrLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # inputs: [E, \delta E]
        p3 = inputs[:, :3]
        dE = inputs[:, 3]
        dE = tf.expand_dims(dE, axis=-1) 
        energy_sqr = tf.square(mh) + tf.reduce_sum(tf.square(p3), axis=1)           
        energy = tf.sign(energy_sqr)*tf.sqrt(energy_sqr)
        energy = tf.expand_dims(energy, axis=-1) 
        output = tf.concat([energy, dE], axis=-1)
        return output