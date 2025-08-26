import tensorflow as tf
import tf2onnx
import onnx

from ModelUtils import LoadModel, ConvertToONNX

def main():
    keras_model, training_params = LoadModel('model_cfg.yaml', allow_missing_train_cfg=False)
    print(keras_model.summary())

    input_signature = [tf.TensorSpec([None, 157], tf.float32, name='input')]
    onnx_model = ConvertToONNX(keras_model, input_signature, out_path='predict_quantiles3D_DL_v4')

    onnx.checker.check_model(onnx_model)
    print(onnx.helper.printable_graph(onnx_model.graph))


if __name__ == '__main__':
    main()