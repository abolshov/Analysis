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


class PtLayer(tf.keras.layers.Layer):
    def __init__(self, name='pt_layer', **kwargs):
        super(PtLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # inputs: pt
        p2 = inputs[:, :2]
        return tf.sqrt(tf.reduce_sum(tf.square(p2), axis=1))


class EtaLayer(tf.keras.layers.Layer):
    def __init__(self, name='eta_layer', **kwargs):
        super(EtaLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # outputs: eta
        p3 = inputs[:, :3]
        mod = tf.sqrt(tf.reduce_sum(tf.square(p3), axis=1))
        output = tf.atanh(p3[:, -1]/mod) 
        return output


class PhiLayer(tf.keras.layers.Layer):
    def __init__(self, name='phi_layer', **kwargs):
        super(PhiLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # outputs: phi
        p3 = inputs[:, :3]
        output = tf.atan2(p3[:, 1], p3[:, 0])
        output = pi*tf.nn.tanh(output)
        return output


class EnergyErrorLayer(tf.keras.layers.Layer):
    def __init__(self, name='energy_error_layer', **kwargs):
        super(EnergyErrorLayer, self).__init__(name=name, **kwargs)

    def call(self, inputs):
        # inputs: px, py, pz, \delta E
        # outputs: [E, \delta E]
        p3 = inputs[:, :3]
        dE = inputs[:, 3]
        dE = tf.nn.softplus(dE)
        dE = tf.expand_dims(dE, axis=-1) 
        energy_sqr = tf.square(mh) + tf.reduce_sum(tf.square(p3), axis=1)           
        energy = tf.sign(energy_sqr)*tf.sqrt(energy_sqr)
        energy = tf.expand_dims(energy, axis=-1) 
        output = tf.concat([energy, dE], axis=-1)
        return output


@tf.keras.utils.register_keras_serializable('Attention')
class Attention(tf.keras.layers.Layer):
    def __init__(self, **kwargs):
        super().__init__()
        self.mha = tf.keras.layers.MultiHeadAttention(**kwargs)
        self.layernorm = tf.keras.layers.LayerNormalization()
        self.add = tf.keras.layers.Add()

    def call(self, x):
        attn_output = self.mha(query=x, value=x, key=x)
        # print(attn_output.shape) # Expected: (batch_size, seq_length, key_dim)
        x = self.add([x, attn_output])
        x = self.layernorm(x)
        # print(x.shape) # Expected: (batch_size, seq_length, key_dim)
        return x


@tf.keras.utils.register_keras_serializable('FeedForward')
class FeedForward(tf.keras.layers.Layer):
    def __init__(self, *, d_model, dff, dropout_rate, **kwargs):
        super().__init__()
        self.seq = tf.keras.Sequential([
            tf.keras.layers.Dense(dff, activation='gelu'),
            tf.keras.layers.Dense(d_model),
            tf.keras.layers.Dropout(dropout_rate)
        ])
        self.add = tf.keras.layers.Add()
        self.layer_norm = tf.keras.layers.LayerNormalization()

    def call(self, x):
        x = self.add([x, self.seq(x)])
        x = self.layer_norm(x) 
        return x


@tf.keras.utils.register_keras_serializable('EncoderLayer')
class EncoderLayer(tf.keras.layers.Layer):
    def __init__(self, *, d_model, num_heads, dff, dropout_rate, **kwargs):
        super().__init__()

        self.self_attention = Attention(num_heads=num_heads, key_dim=d_model, dropout=dropout_rate)
        self.feedforward = FeedForward(d_model=d_model, dff=dff, dropout_rate=dropout_rate)

    def call(self, x):
        x = self.self_attention(x) # Expected x.shape = (batch_size, seq_length, key_dim)
        x = self.feedforward(x)
        return x


@tf.keras.utils.register_keras_serializable('EmbeddingLayer')
class EmbeddingLayer(tf.keras.layers.Layer):
    def __init__(self, *, dim_embedding, dropout_rate, **kwargs):
        """
        Args:
            num_layers: number of layers in embedding perceptron
            num_units: number of units in layers
            dim_embedding: dimension of the embdedding vector
        """
        super().__init__()
        
        self.seq = tf.keras.Sequential([tf.keras.layers.Dense(dim_embedding, activation='gelu')])
        self.seq.add(tf.keras.layers.LayerNormalization())
        if dropout_rate and dropout_rate > 0.0:
            self.seq.add(tf.keras.layers.Dropout(dropout_rate))

    def call(self, x):
        out = self.seq(x)
        return out


@tf.keras.utils.register_keras_serializable('Embedding')
class Embedding(tf.keras.layers.Layer):
    def __init__(self, *, dim_embedding, dropout_rate, **kwargs):
        """
        Args:
            dim_embedding: dimension of the hidden respresentation vector
            dropout_rate: probability of dropout
        """
        super().__init__()

        self.lep_embedding = EmbeddingLayer(dim_embedding=dim_embedding, dropout_rate=dropout_rate)
        self.met_embedding = EmbeddingLayer(dim_embedding=dim_embedding, dropout_rate=dropout_rate)
        self.jet_embedding = EmbeddingLayer(dim_embedding=dim_embedding, dropout_rate=dropout_rate)
        self.fatjet_embedding = EmbeddingLayer(dim_embedding=dim_embedding, dropout_rate=dropout_rate)

    def call(self, inputs):
        lep_input, met_input, jet_input, fatjet_input = inputs
        lep_out = self.lep_embedding(lep_input)
        met_out = self.met_embedding(met_input)
        jet_out = self.jet_embedding(jet_input)
        fatjet_out = self.fatjet_embedding(fatjet_input)
        out = tf.stack([lep_out, met_out, jet_out, fatjet_out], axis=1)
        return out


@tf.keras.utils.register_keras_serializable('Encoder')
class Encoder(tf.keras.layers.Layer):
    def __init__(self, *, 
                 num_encoder_layers, 
                 d_model,
                 num_heads, 
                 dff,
                 dropout_rate=0.01, 
                 name='encoder', 
                 **kwargs): 
        """
        Implements encoder part of transformer. Consists of num_encoder_layers identical layers. Each layer contains
        MultihedAttention layer and feedforward network (2 dense layers)

        Args:
            d_model: dimension of keys and queries
            num_heads: number of attention heads
            dff: number of units in the first dense layer of the feed-froward network following attention heads
            dropout: dropout rate applied in feedforward network
        """
        super().__init__()
        self.num_encoder_layers = num_encoder_layers

        # project input to "vector representation" (batch_size, n_features) -> (batch_size, seq_len, d_model)
        # for us seq_len = 4 [leptons, met, ak4_jets, ak8_jets]
        self.embedding = Embedding(dim_embedding=d_model, dropout_rate=dropout_rate)
        
        self.encoder_layers = [EncoderLayer(d_model=d_model, 
                                            num_heads=num_heads, 
                                            dff=dff, 
                                            dropout_rate=dropout_rate) for _ in range(num_encoder_layers)]
        self.dropout = tf.keras.layers.Dropout(dropout_rate)

    def call(self, x):
        x = self.embedding(x)
        x = self.dropout(x)
        for i in range(self.num_encoder_layers):
            x = self.encoder_layers[i](x)
        return x


@tf.keras.utils.register_keras_serializable('MultiheadQuantileRegressor')
class MultiheadQuantileRegressor(tf.keras.layers.Layer):
    def __init__(self, *, 
                 num_quantiles, 
                 target_names, 
                 use_energy_layer, 
                 use_quantile_ordering,
                 means=None,
                 scales=None,
                 **kwargs):
        super().__init__()

        self.heads = []
        for idx, target_name in enumerate(target_names):
            if use_energy_layer and 'E' in target_name:
                continue

            if use_quantile_ordering:
                self.heads.append(tf.keras.Sequential([tf.keras.layers.Dense(num_quantiles),
                                                       QuantileOrderingLayer(name=target_name)]))
            else:
                self.heads.append(tf.keras.layers.Dense(num_quantiles, name=target_name))

        self.hbb_energy_layer = None
        self.hvv_energy_layer = None
        self.hvv_en_idx = None
        self.hbb_en_idx = None

        if use_energy_layer:
            self.hvv_en_idx = target_names.index('genHVV_E')
            self.hbb_en_idx = target_names.index('genHbb_E')

            normalize = means is not None and scales is not None
            self.hbb_energy_layer = EnergyLayer(num_quantiles=num_quantiles, 
                                                normalize=normalize, 
                                                means=means[4:] if normalize else None,
                                                scales=scales[4:] if normalize else None,
                                                name='genHbb_E')
            
            self.hvv_energy_layer = EnergyLayer(num_quantiles=num_quantiles, 
                                                normalize=normalize, 
                                                means=means[:4] if normalize else None,
                                                scales=scales[:4] if normalize else None,
                                                name='genHVV_E')

    def call(self, inputs): 
        # inputs: list of (batch_size, num_objects * d_model), len = len(self.heads)
        # outputs: list of (batch_size, num_quantiles)
        outputs = []
        for i in range(len(self.heads)):
            out = self.heads[i](inputs)
            outputs.append(out)

        if self.hvv_en_idx and self.hbb_en_idx:
            hbb_input = tf.concat(outputs[3:], axis=-1)
            hvv_input = tf.concat(outputs[:3], axis=-1)

            outputs.insert(self.hvv_en_idx, self.hvv_energy_layer(hvv_input))
            outputs.insert(self.hbb_en_idx, self.hbb_energy_layer(hbb_input))

        return outputs
