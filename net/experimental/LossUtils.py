import tensorflow as tf 
import numpy as np

mh = tf.constant(125.0)

# @tf.keras.saving.register_keras_serializable(package="CombinedLoss")
class CombinedLoss(tf.keras.losses.Loss):
    def __init__(self, strength=0.01, name="combined_loss", **kwargs):
        super().__init__(name=name)
        self.strength = strength
        self.logcosh_loss = tf.keras.losses.LogCosh() 
        self.mse_loss = tf.keras.losses.MeanSquaredError()

    def call(self, y_true, y_pred):
        logcosh = self.logcosh_loss(y_true, y_pred)
        mbb2 = tf.square(y_pred[:, 3]) - tf.reduce_sum(tf.square(y_pred[:, 0:3]), axis=1)
        mww2 = tf.square(y_pred[:, 7]) - tf.reduce_sum(tf.square(y_pred[:, 4:7]), axis=1)
        mbb = tf.sign(mbb2)*tf.sqrt(tf.abs(mbb2))
        mww = tf.sign(mww2)*tf.sqrt(tf.abs(mww2))
        penalty = 0.5*self.mse_loss(mh, mbb) + 0.5*self.mse_loss(mh, mww)
        total_loss = logcosh + self.strength*penalty
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
        mbb2 = tf.square(y_pred[:, 3]) - tf.reduce_sum(tf.square(y_pred[:, 0:3]), axis=1)
        mww2 = tf.square(y_pred[:, 7]) - tf.reduce_sum(tf.square(y_pred[:, 4:7]), axis=1)
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


class Momentum3DLoss(tf.keras.losses.Loss):
    def __init__(self, relative_weights={'energy': 0.225, 'momentum': 0.675, 'mass': 0.1}, name="momentum3D_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.loss_func = tf.keras.losses.LogCosh()
        self.energy_weight = relative_weights['energy']
        self.momentum_weight = relative_weights['momentum']
        self.mass_weight = relative_weights['mass']

    def call(self, y_true, y_pred):
        pred_hbb_energy_sqr = tf.square(mh) + tf.reduce_sum(tf.square(y_pred[:, 0:3]), axis=1)
        pred_hvv_energy_sqr = tf.square(mh) + tf.reduce_sum(tf.square(y_pred[:, 3:6]), axis=1)
        pred_hbb_energy = tf.sign(pred_hbb_energy_sqr)*tf.sqrt(tf.abs(pred_hbb_energy_sqr))
        pred_hvv_energy = tf.sign(pred_hvv_energy_sqr)*tf.sqrt(tf.abs(pred_hvv_energy_sqr))

        pred_x_energy = pred_hbb_energy + pred_hvv_energy
        pred_x_p3 = y_pred[:, 0:3] + y_pred[:, 3:6]
        pred_x_mass_sqr = tf.square(pred_x_energy) - tf.reduce_sum(tf.square(pred_x_p3), axis=1)
        pred_x_mass = tf.sign(pred_x_mass_sqr)*tf.sqrt(tf.abs(pred_x_mass_sqr))

        true_x_energy = y_true[:, 3] + y_true[:, 7]
        true_x_p3 = y_true[:, 0:3] + y_true[:, 4:7]
        true_x_mass = tf.sqrt(tf.square(true_x_energy) - tf.reduce_sum(tf.square(true_x_p3), axis=1))

        # make sure format of labels is (px, py, pz, E)
        energy_loss = (self.loss_func(y_true[:, 3], pred_hbb_energy) + self.loss_func(y_true[:, 7], pred_hvv_energy))/2
        momentum_loss = (self.loss_func(y_true[:, 0:3], y_pred[:, 0:3]) + self.loss_func(y_true[:, 4:7], y_pred[:, 3:6]))/2
        mass_loss = self.loss_func(true_x_mass, pred_x_mass)
        total_loss = self.energy_weight*energy_loss + self.momentum_weight*momentum_loss + self.mass_weight*mass_loss
        return total_loss

    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super().get_config()
        config.update({'energy_weight': self.energy_weight, 
                       'momentum_weight': self.momentum_weight,
                       'mass_weight': self.mass_weight})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)


class QuantileLoss(tf.keras.losses.Loss):
    """
    Loss for predicting array of quantiles for single variable
    """
    def __init__(self, quantiles=[0.5], noncrossing=True, name="quantile_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.quantiles = tf.constant(quantiles, dtype=tf.float32)
        self.noncrossing = noncrossing

    def call(self, y_true, y_pred):
        error = y_true - y_pred
        loss = tf.reduce_mean(tf.maximum(self.quantiles*error, error*(self.quantiles - 1)), axis=-1)
        
        if self.noncrossing:
            # # solution: https://github.com/scandido/quantile-regression-tf/blob/main/examples/Disallow_Crossing_Quantiles.ipynb
            # def mask(x):
            #     l = x.shape[-1]
            #     k = x.shape[-1] // 2
            #     m = tf.concat([-1*tf.ones(l - k - 1), tf.ones(l - k)], axis=0)
            #     return m * x

            # def sum_from_middle(x):
            #     d = len(x.shape) - 1
            #     knot = self.quantiles.shape[-1] // 2
            #     a = tf.cumsum(x[..., :knot + 1], axis=d, reverse=True)
            #     b = tf.cumsum(x[..., knot:], axis=d)
            #     return tf.concat([a[..., :-1], b], axis=d)

            # loss = mask(error)
            # loss = sum_from_middle(loss)

            quantile_gap = y_pred[:, 1:] - y_pred[:, :-1]
            quantile_order_penalty = tf.reduce_mean(0.5*(1 - tf.sign(quantile_gap))*tf.abs(quantile_gap), axis=-1)
            quantile_width_penalty = tf.abs(y_pred[:, 0] - y_pred[:, -1])
            penalty = quantile_order_penalty + 0.1*quantile_width_penalty
            return loss + penalty

        return loss

    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super().get_config()
        config.update({'quantiles': self.quantiles.numpy(),
                       'noncrossing': self.noncrossing})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)