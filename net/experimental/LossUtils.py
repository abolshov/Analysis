import tensorflow as tf 
import numpy as np

mh = tf.constant(125.0)
pi = tf.constant(np.pi)

# @tf.keras.saving.register_keras_serializable
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


@tf.keras.utils.register_keras_serializable('QuantileLoss')
class QuantileLoss(tf.keras.losses.Loss):
    """
    Loss for predicting array of quantiles for single variable
    """
    def __init__(self, quantiles=[0.5], order_penalty_rate=None, width_penalty_rate=None, name="quantile_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.quantiles = tf.constant(quantiles, dtype=tf.float32)
        self.order_penalty_rate = order_penalty_rate if order_penalty_rate and len(quantiles) > 1 else None
        self.width_penalty_rate = width_penalty_rate if width_penalty_rate and len(quantiles) > 1 else None

    def call(self, y_true, y_pred):
        error = y_true - y_pred
        loss = tf.reduce_mean(tf.maximum(self.quantiles*error, error*(self.quantiles - 1)), axis=-1)
        
        if self.order_penalty_rate:
            quantile_gap = y_pred[:, 1:] - y_pred[:, :-1]
            order_penalty = tf.reduce_mean(0.5*(1 - tf.sign(quantile_gap))*tf.abs(quantile_gap), axis=-1)
            return loss + self.order_penalty_rate*order_penalty

        if self.width_penalty_rate:
            width_penalty = self.width_penalty_rate*tf.abs(y_pred[:, -1] - y_pred[:, 0])
            return loss + width_penalty

        return loss

    def get_config(self):
        """
        Returns the serializable configuration of the layer.
        This is necessary for saving and loading the model with custom objects.
        """
        config = super().get_config()
        config.update({'quantiles': tf.keras.utils.serialize_keras_object(self.quantiles),
                       'order_penalty_rate': self.order_penalty_rate,
                       'width_penalty_rate': self.width_penalty_rate})
        return config

    @classmethod
    def from_config(cls, config):
        quantiles = config.pop('quantiles')
        return cls(quantiles=quantiles['config']['value'], **config)


class LogPtLoss(tf.keras.losses.Loss):
    def __init__(self, name="log_pt_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.log_cosh = tf.keras.losses.LogCosh()
        self.epsilon = tf.constant(0.1, dtype=tf.float32)

    def call(self, y_true, y_pred):
        return self.log_cosh(tf.math.log(y_true + self.epsilon), tf.math.log(y_pred + self.epsilon))



class DeltaPhiLoss(tf.keras.losses.Loss):
    def __init__(self, name="delta_phi_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.log_cosh = tf.keras.losses.LogCosh()

    def call(self, y_true, y_pred):
        dphi = y_true - y_pred
        dphi = tf.where(dphi > pi, dphi - 2*pi, dphi)
        dphi = tf.where(dphi <= -pi, dphi + 2*pi, dphi)
        loss = self.log_cosh(dphi, 0.0)
        return loss


class CombinedEnergyLoss(tf.keras.losses.Loss):
    def __init__(self, quantile=None, log=False, name="combined_energy_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.log_cosh = tf.keras.losses.LogCosh()
        self.quantile = quantile
        self.log = log

    def call(self, y_true, y_pred):
        # y_true - only contains ground truth for energy
        # there is no ground truth for energy error
        energy_loss = None
        if self.log:
            energy_loss = self.log_cosh(tf.math.log(y_true), tf.math.log(y_pred[:, 0]))
        else:
            energy_loss = self.log_cosh(y_true, y_pred[:, 0])
        
        if self.quantile:
            def Heaviside(x):
                return 0.5*(1 + tf.sign(x))

            def QR(q, x, mu):
                return (x - mu)*(q*Heaviside(x - mu) + (1 - q)*Heaviside(mu - x))

            dE_pred = y_pred[:, 1]
            qr1 = QR(0.5*(1 - self.quantile), tf.math.log(y_true), tf.math.log(y_pred[:, 0]) - tf.math.log(dE_pred))
            qr2 = QR(0.5*(1 + self.quantile), tf.math.log(y_true), tf.math.log(y_pred[:, 0]) + tf.math.log(dE_pred))
            energy_error_loss = qr1 + qr2
            return energy_loss + energy_error_loss
        else:
            return energy_loss
        

    def get_config(self):
        config = super().get_config()
        config.update({'quantile': self.quantile,
                       'log': self.log})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)


class PtEtaPhiELoss(tf.keras.losses.Loss):
    def __init__(self, quantile=None, log=False, name="PtEtaPhiELoss_loss", **kwargs):
        super().__init__(name=name, **kwargs)
        self.log_cosh = tf.keras.losses.LogCosh()
        self.quantile = quantile
        self.log = log

    def call(self, y_true, y_pred):
        # y_true: [pt, eta, phi, E]
        # y_pred: [px, py, pz, dE]
        
        energy_loss = None
        en_pred = tf.sqrt(tf.square(mh) + tf.reduce_sum(tf.square(y_pred[:, :3]), axis=1))
        if self.log:
            energy_loss = self.log_cosh(tf.math.log(y_true[:, 3] + 10), tf.math.log(en_pred + 10))
        else:
            energy_loss = self.log_cosh(y_true[:, 3], en_pred)
        
        if self.quantile:
            def Heaviside(x):
                return 0.5*(1 + tf.sign(x))

            def QR(q, x, mu):
                return (x - mu)*(q*Heaviside(x - mu) + (1 - q)*Heaviside(mu - x))

            dE_pred = tf.nn.softplus(y_pred[:, 3])
            qr1 = QR(0.5*(1 - self.quantile), tf.math.log(y_true[:, 3] + 10), tf.math.log(en_pred + 10) - tf.math.log(dE_pred + 10))
            qr2 = QR(0.5*(1 + self.quantile), tf.math.log(y_true[:, 3] + 10), tf.math.log(en_pred + 10) + tf.math.log(dE_pred + 10))
            energy_error_loss = qr1 + qr2
            energy_loss = energy_loss + energy_error_loss

        pt_loss = None
        p2_pred = y_pred[:, :2]
        pt_pred = tf.sqrt(tf.reduce_sum(tf.square(p2_pred), axis=1))
        # tf.print(pt_pred)
        pt_true = y_true[:, 0]
        # tf.print(pt_true)
        if self.log:
            pt_loss = self.log_cosh(tf.math.log(pt_true + 10), tf.math.log(pt_pred + 10))
        else:
            pt_loss = self.log_cosh(pt_true, pt_pred)

        # tf.print(pt_loss)

        eta_loss = None
        p3_pred = y_pred[:, :3]
        mod_pred = tf.sqrt(tf.reduce_sum(tf.square(p3_pred), axis=1))
        eta_pred = tf.atanh(p3_pred[:, -1]/mod_pred) 
        eta_true = y_true[:, 1]
        eta_loss = self.log_cosh(eta_true, eta_pred)

        # tf.print(eta_loss[0])

        phi_loss = None
        phi_pred = tf.atan2(p3_pred[:, 1], p3_pred[:, 0])
        dphi = y_true[:, 2] - phi_pred
        dphi = tf.where(dphi > pi, dphi - 2*pi, dphi)
        dphi = tf.where(dphi <= -pi, dphi + 2*pi, dphi)
        phi_loss = self.log_cosh(dphi, 0.0)

        # print('phi_loss:')
        # tf.print(phi_loss[0])

        # print('----------------------')

        total_loss = pt_loss + energy_loss + eta_loss + phi_loss
        return total_loss
        

    def get_config(self):
        config = super().get_config()
        config.update({'quantile': self.quantile,
                       'log': self.log})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)

@tf.keras.utils.register_keras_serializable('MultiheadLoss')
class MultiheadLoss(tf.keras.losses.Loss):
    def __init__(self, num_quantiles=1, scales=None, means=None, name='multihead_loss', **kwargs):
        super().__init__(name=name, **kwargs)
        self.mae = tf.keras.losses.MAE
        self.num_quantiles = num_quantiles
        self.scales = scales
        self.means = means

    def call(self, y_true, y_pred):
        # each head outputs num_quantiles predictions for 1 variable
        # y_pred: list of 8 tensors (batch_size, num_quantiles) each 
        # y_true: (batch_size, 8)

        # central quantile must be in the middle
        central_quantile = self.num_quantiles // 2
        central_pred = y_pred[:, central_quantile::self.num_quantiles]

        if self.means is not None and self.scales is not None:
            central_pred = central_pred*self.scales + self.means
            y_true = y_true*self.scales + self.means

        true_hbb = y_true[:, 4:]
        true_hvv = y_true[:, :4]
        true_x = true_hbb + true_hvv
        true_x_mass = tf.sqrt(tf.square(true_x[:, 3]) - tf.reduce_sum(tf.square(true_x[:, :3]), axis=1))

        pred_hbb = central_pred[:, 4:]
        pred_hvv = central_pred[:, :4]
        pred_x = pred_hbb + pred_hvv
        pred_x_mass_sqr = tf.square(pred_x[:, 3]) - tf.reduce_sum(tf.square(pred_x[:, :3]), axis=1)
        pred_x_mass = tf.sign(pred_x_mass_sqr)*(tf.sqrt(tf.abs(pred_x_mass_sqr)))
        
        loss = self.mae(true_x_mass, pred_x_mass)
        return loss

    def get_config(self):
        config = super().get_config()
        config.update({'num_quantiles': self.num_quantiles,
                       'scales': self.scales,
                       'means': self.means})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)

        
