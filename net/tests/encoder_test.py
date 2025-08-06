import tensorflow as tf
from LayerUtils import Encoder, Embedding

def encoder_test():
    print('Encoder test')
    t = tf.random.normal(shape=(10, 157))
    encoder = Encoder(num_encoder_layers=6, 
                      num_proj_layers=2,
                      dim_proj_layer=512,
                      d_model=512,
                      num_heads=8, 
                      dff=2048,
                      dropout_rate=0.02)
    out = encoder(t)
    print(f'Encoder output shape: {out.shape}')

def embedding_test():
    print('Embedding test')
    t = tf.random.normal(shape=(10, 157))
    embedding = Embedding(num_layers=3, 
                          dim_layer=512,
                          dim_embedding=512)
    out = embedding(t)
    print(f'Embedding output shape: {out.shape}')

def main():
    encoder_test()
    embedding_test()

if __name__ == '__main__':
    main()