import tensorflow as tf
import yaml
import importlib
import os
import onnx
import tf2onnx


def LoadModel(cfg_file_name):
    cfg = None
    with open(cfg_file_name, 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    model = None
    model_dir = cfg['directory']
    model_name = cfg['name']
    model_fmt = cfg['format']
    custom_objects = {}
    if 'custom_objects' in cfg and cfg['custom_objects']:
        module_mapping = cfg['custom_objects']
        for object_name, module_name in module_mapping.items():  
            module = importlib.import_module(module_name)
            custom_objects[object_name] = getattr(module, object_name)

    path_to_model = os.path.join(model_dir, f'{model_name}.{model_fmt}')
    model = tf.keras.models.load_model(path_to_model, custom_objects=custom_objects)

    # load training params if available
    training_params = None
    with open(os.path.join(model_dir, 'params.yaml'), 'r') as train_cfg_file:
        training_params = yaml.safe_load(train_cfg_file)
    return model, training_params

def ConvertToONNX(model, input_signature, out_path=None, op_set=13):
    """
    input_signature ~ [tf.TensorSpec([None, 157], tf.float32, name='input')]
    """

    onnx_model, _ = tf2onnx.convert.from_keras(model,
                                               input_signature=input_signature,
                                               opset=op_set)
    if out_path:
        onnx.save(onnx_model, os.path.join(out_path, f'{model.name}.onnx'))
    return onnx_model
