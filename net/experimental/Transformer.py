import tensorflow as tf

from experimental.LayerUtils import Encoder, MultiheadQuantileRegressor


@tf.keras.utils.register_keras_serializable('Transformer')
class Transformer(tf.keras.Model):
    """
    Transformer mass estimator model. Takes arguments of encoder block and dictionary with arguments for regressor.
    To enable flexibility for regressor, its args must be passed in dictionary
    Args:
        num_encoder_layers:
        d_model: dimension of embedding, keys and queries
        num_heads: number of attention heads
        dff: dimension of feed-forward network in encoder
        dropout_rate: dropout proba
        regressor_cfg: dictionary with arguments of regressor being used
    """
    def __init__(self, *, 
                 num_encoder_layers, 
                 d_model, 
                 num_heads, 
                 dff, 
                 dropout_rate=0.01, 
                 regressor_cfg={},
                 **kwargs):
        super().__init__()

        # save arguments to put them in config
        self.num_encoder_layers = num_encoder_layers
        self.d_model = d_model
        self.num_heads = num_heads
        self.dff = dff 
        self.dropout_rate = dropout_rate
        self.regressor_cfg = regressor_cfg

        self.encoder = Encoder(num_encoder_layers=num_encoder_layers, 
                               d_model=d_model, 
                               num_heads=num_heads, 
                               dff=dff)
        self.regressor = MultiheadQuantileRegressor(**regressor_cfg)

    def call(self, inputs):
        # inputs: list of tensors corresponding to features of different objects
        # embedded featrues: (batch_size, num_objects, embedding_dim) = (batch_size, num_objects, d_model)
        encoded_features = self.encoder(inputs) # (batch_size, num_objects, d_model)
        flattened_features = tf.keras.layers.Flatten()(encoded_features) # (batch_size, num_objects * d_model)
        predictions = self.regressor(flattened_features) # list of (batch_size, num_quantiles)
        return predictions

    def get_config(self):
        config = super().get_config()
        config.update({'num_encoder_layers': self.num_encoder_layers,
                       'd_model': self.d_model,
                       'num_heads': self.num_heads,
                       'dff': self.dff,
                       'dropout_rate': self.dropout_rate,
                       'regressor_cfg': self.regressor_cfg})
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)