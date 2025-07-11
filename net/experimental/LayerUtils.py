import tensorflow as tf 

# if batch normalization is used, this value leads to errors
# it must be replaced to Higgs boson rest energy z-score of the train set distribution
mh = tf.constant(125.0)

class EnergyLayer(tf.keras.layers.Layer):
    def __init__(self, num_quantiles=1, name='energy_layer', **kwargs):
        super(EnergyLayer, self).__init__(name=name, **kwargs)
        self.num_quantiles = num_quantiles

    def call(self, inputs):
        """
        assume that input shape is (num_events, num_quantiles*3)
        *3 because px, py, pz

        assume layout is [px_q16, px_q50, px_q84, ..., pz_q16, pz_q50, pz_q84]
        so that variables for the same quantile are within displacement 3 from each other
        """

        energies = []
        for q in range(self.num_quantiles):
            p3 = inputs[:, q::self.num_quantiles]
            energy_q_sqr = tf.square(mh) + tf.reduce_sum(tf.square(p3), axis=1)           
            energy_q = tf.sign(energy_q_sqr)*tf.sqrt(energy_q_sqr)
            energy_q = tf.expand_dims(energy_q, axis=-1) 
            energies.append(energy_q)

        output = tf.concat(energies, axis=-1)
        return output

    def get_config(self):
        config = super().get_config()
        config.update({'num_quantiles': self.num_quantiles})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)
            