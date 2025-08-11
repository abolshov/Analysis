import tensorflow as tf

from LayerUtils import Encoder, Embedding 

class Transformer(tf.keras.Model):
    def __init__(self, *, num_encoder_layers, d_model, num_heads, dff, dropout_rate=0.01):
        super().__init__()

        self.embedding = Embedding(dim_embedding=d_model, dropout_rate=dropout_rate)
        self.encoder = Encoder(num_encoder_layers=num_encoder_layers, d_model=d_model, num_heads=num_heads, dff=dff)
        self.regressor = 

    def call(self, inputs):
        # inputs: list of tensors corresponding to features of different objects
        # embedded featrues: (batch_size, num_objects, embedding_dim) = (batch_size, num_objects, d_model)
        embedded_features = self.embdedding(inputs)
        encoded_features = self.encoder(embedded_features)
        predictions = self.regressor(encoded_features)
        return predictions
