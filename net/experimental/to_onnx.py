import tensorflow as tf
import tf2onnx
import onnx

from ModelUtils import LoadModel, ConvertToONNX

def main():
    keras_model, training_params = LoadModel('model_cfg.yaml', allow_missing_train_cfg=False)
    print(keras_model.summary())

    num_features = training_params['input_shape'][0]
    input_signature = [tf.TensorSpec([None, num_features], tf.float32, name='input')]
    model_dir = training_params['model_dir']
    onnx_model = ConvertToONNX(keras_model, input_signature, model_dir)

    onnx.checker.check_model(onnx_model)
    print(onnx.helper.printable_graph(onnx_model.graph))


if __name__ == '__main__':
    main()