import tensorflow as tf
import yaml
import importlib
import os


def LoadModel(cfg_file_name):
    cfg = None
    with open (cfg_file_name, 'r') as cfg_file:
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
            custom_objects[object_name] =  getattr(module, object_name)

    path_to_model = os.path.join(model_dir, f'{model_name}.{model_fmt}')
    model = tf.keras.models.load_model(path_to_model, custom_objects=custom_objects)
    return model