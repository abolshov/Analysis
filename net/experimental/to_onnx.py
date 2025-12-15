import tensorflow as tf
import tf2onnx
import onnx
import yaml

from ModelUtils import LoadModel, ConvertToONNX

def main():
    cfg = None
    with open('model_cfg.yaml', 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    keras_model, training_params = LoadModel(cfg, allow_missing_train_cfg=False)
    print(keras_model.summary())

    num_features = training_params['input_shape'][0]
    input_signature = [tf.TensorSpec([None, num_features], tf.float32, name='input')]
    model_dir = training_params['model_dir']
    onnx_name = cfg['onnx_name']
    onnx_model = ConvertToONNX(keras_model, input_signature, save_to=model_dir, onnx_name=onnx_name)

    onnx.checker.check_model(onnx_model)
    print(onnx.helper.printable_graph(onnx_model.graph))


if __name__ == '__main__':
    main()