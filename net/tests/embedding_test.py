import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import tensorflow as tf
from experimental.LayerUtils import Embedding, EmbeddingLayer

def main():
    lep_tensor = tf.random.normal(shape=(10, 8))
    lep_embedding = EmbeddingLayer(dim_embedding=512, dropout_rate=0.01)
    lep_out = lep_embedding(lep_tensor)
    print(f'Lepton embedding output shape: {lep_out.shape}')
    print(f'Lepton embedding param count: {lep_embedding.count_params()}')

    jet_tensor = tf.random.normal(shape=(10, 120))
    jet_embedding = EmbeddingLayer(dim_embedding=512, dropout_rate=0.01)
    jet_out = jet_embedding(jet_tensor)
    print(f'Jet embedding output shape: {jet_out.shape}')
    print(f'Jet embedding param count: {jet_embedding.count_params()}')

    fatjet_tensor = tf.random.normal(shape=(10, 27))
    fatjet_embedding = EmbeddingLayer(dim_embedding=512, dropout_rate=0.01)
    fatjet_out = fatjet_embedding(fatjet_tensor)
    print(f'Fatjet embedding output shape: {fatjet_out.shape}')
    print(f'Fatjet embedding param count: {fatjet_embedding.count_params()}')

    met_tensor = tf.random.normal(shape=(10, 2))
    met_embedding = EmbeddingLayer(dim_embedding=512, dropout_rate=0.01)
    met_out = met_embedding(met_tensor)
    print(f'MET embedding output shape: {met_out.shape}')
    print(f'MET embedding param count: {met_embedding.count_params()}')

    particle_embedding = Embedding(dim_embedding=512, dropout_rate=0.01)
    particle_out = particle_embedding([lep_tensor, met_tensor, jet_tensor, fatjet_tensor])
    print(f'Particle embedding output shape: {particle_out.shape}')
    print(f'Particle embedding param count: {particle_embedding.count_params()}')

if __name__ == '__main__':
    main()