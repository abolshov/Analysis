import tensorflow as tf

from experimental.LayerUtils import Encoder, Embedding, MultiheadQuantileRegressor

target_names = ['genHVV_px', 'genHVV_py', 'genHVV_pz', 'genHVV_E',
                'genHbb_px', 'genHbb_py', 'genHbb_pz', 'genHbb_E']

class Transformer(tf.keras.Model):
    def __init__(self, *, num_encoder_layers, d_model, num_heads, dff, dropout_rate=0.01, **kwargs):
        super().__init__()

        self.encoder = Encoder(num_encoder_layers=num_encoder_layers, d_model=d_model, num_heads=num_heads, dff=dff)
        self.regressor = MultiheadQuantileRegressor(num_quantiles=3, 
                                                    target_names=target_names, 
                                                    use_energy_layer=True,
                                                    use_quantile_ordering=False)

    def call(self, inputs):
        # inputs: list of tensors corresponding to features of different objects
        # embedded featrues: (batch_size, num_objects, embedding_dim) = (batch_size, num_objects, d_model)
        encoded_features = self.encoder(inputs) # (batch_size, num_objects, d_model)
        flattened_features = tf.keras.layers.Flatten()(encoded_features) # (batch_size, num_objects * d_model)
        predictions = self.regressor(flattened_features) # list of (batch_size, num_quantiles)
        return predictions
