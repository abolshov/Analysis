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
    def __init__(self, num_particles=2, quantiles=[0.5], name="quantile_loss", **kwargs):
        """
        num_particles: int - number of particles whoose 3D momentum trying to predict
        """
        super().__init__(name=name, **kwargs)
        self.num_quantiles = len(quantiles)
        # quantiles must have shape (num_quantiles, num_predictions)
        # each row repeats value of quantile as num_predictions times such that each output node computes every quantile
        # but lorentz vector has 4 components and we're predicting only 3
        # but will be comparing full (4D) lorentz vectors - so we'll add computed energy to preducted tensor
        self.num_predictions = 4*num_particles*self.num_quantiles
        self.num_particles = num_particles
        # quantile_array = np.array(quantiles)
        # quantile_array = np.reshape(quantile_array, (len(quantiles), 1))
        # quantile_array = np.repeat(quantile_array, 4*num_particles, axis=-1)
        # quantile_array = np.reshape(quantile_array, (-1, self.num_quantiles, 4*self.num_particles))
        quantile_array = np.repeat(np.array(quantiles)[:, None], 4*self.num_particles, axis=-1)
        quantile_array = quantile_array.reshape(-1, self.num_quantiles, 4*num_particles)
        self.quantiles = tf.constant(quantile_array, dtype=tf.float32)
        # self.quantiles = tf.reshape(self.quantiles, [-1, self.num_quantiles, 4*self.num_particles])
        tf.print(self.quantiles.shape)

    def call(self, y_true, y_pred):
        # y_pred has shape (batch_size, num_particles*num_quantiles)
        # second dimension is a flat array ordered as [px1, py1, pz1, px2, py2, pz2] * num_quantiles
        #                                              hbb  hbb  hbb  hvv  hvv  hvv 

        # need to repeat y_true num_quantiles times - each quantile has true output
        y_gt = tf.concat([y_true]*self.num_quantiles, axis=-1)
        # y_gt = tf.reshape(y_gt, [-1, self.num_quantiles, 4*self.num_particles])

        quantile_predictions = []
        for q in range(self.num_quantiles):
            start = q*self.num_quantiles
            end = start + 3*self.num_particles

            pred_hbb_p3 = y_pred[:, start:start + 3]
            pred_hvv_p3 = y_pred[:, start + 3:end]
            pred_hbb_energy_sqr = tf.square(mh) + tf.reduce_sum(tf.square(pred_hbb_p3), axis=-1)
            pred_hvv_energy_sqr = tf.square(mh) + tf.reduce_sum(tf.square(pred_hvv_p3), axis=-1)
            pred_hbb_energy = tf.sign(pred_hbb_energy_sqr)*tf.sqrt(tf.abs(pred_hbb_energy_sqr))
            pred_hvv_energy = tf.sign(pred_hvv_energy_sqr)*tf.sqrt(tf.abs(pred_hvv_energy_sqr))

            pred_hbb_energy = tf.reshape(pred_hbb_energy, [-1, 1])
            pred_hvv_energy = tf.reshape(pred_hvv_energy, [-1, 1])
            pred_hbb_p4 = tf.concat([pred_hbb_p3, pred_hbb_energy], axis=-1)
            pred_hvv_p4 = tf.concat([pred_hvv_p3, pred_hvv_energy], axis=-1)
            quantile_pred = tf.concat([pred_hbb_p4, pred_hvv_p4], axis=-1)
            quantile_predictions.append(quantile_pred)

        y = tf.concat(quantile_predictions, axis=-1)
        # y = tf.reshape(y, [-1, self.num_quantiles, 4*self.num_particles])
        assert y_gt.shape[-1] == self.quantiles.shape[-1]*self.num_quantiles, ("QuantileLoss: dimension of output vector must match number of predictions", 
                f"{y_gt.shape[-1]} vs {self.quantiles.shape[-1]*self.num_quantiles}")

        error = y_gt - y
        # error must have shape (batch_size, num_quantiles, total_num_components_pred_per_particle)
        # total_num_components_pred_per_particle is always self.num_particles * 4 (bc lorentz vector has 4 components - [px, py, pz, E])
        error = tf.reshape(error, [-1, self.num_quantiles, 4*self.num_particles])
        # reduce over all predictions and, all quantiles and all examples in the batch
        return tf.reduce_mean(tf.maximum(self.quantiles*error, error*(self.quantiles - 1)))


    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super().get_config()
        config.update({'quantiles': self.quantiles.numpy(),
                       'num_quantiles': self.num_quantiles,
                       'num_predictions': self.num_predictions,
                       'num_particles': self.num_particles})
        return config


    @classmethod
    def from_config(cls, config):
        return cls(**config)