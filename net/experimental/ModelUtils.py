import tensorflow as tf
import yaml
import importlib
import os
import onnx
import tf2onnx
import sys


def LoadModel(cfg_file_name, allow_missing_train_cfg=True):
    """
        function for loading .keras model
        returns: tuple (model, training_parameters)
        if allow_missing_train_cfg is True, training_parameters will be None if .yaml with is not found, otherwise exception will be re-thrown
    """
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
    training_cfg_name = f'params_{model_name}.yaml'
    traing_cfg_path = os.path.join(model_dir, training_cfg_name)
    try:
        with open(traing_cfg_path, 'r') as train_cfg_file:
            training_params = yaml.safe_load(train_cfg_file)
    except FileNotFoundError:
        print(f'Model training config file {training_cfg_name} not found in {model_dir}')
        if not allow_missing_train_cfg:
            raise
    return model, training_params

def ConvertToONNX(model, input_signature, save_to=None, op_set=13):
    """
    input_signature ~ [tf.TensorSpec([None, 157], tf.float32, name='input')]
    if save_to is None, model is not saved
    """

    onnx_model, _ = tf2onnx.convert.from_keras(model,
                                               input_signature=input_signature,
                                               opset=op_set)
    if save_to:
        onnx.save(onnx_model, os.path.join(save_to, f'{model.name}.onnx'))
    return onnx_model
