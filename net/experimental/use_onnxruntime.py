import onnxruntime as ort
import numpy as np
import yaml
import os

from Dataloader import Dataloader
from PlotUtils import PlotHist
from PrepUtils import ComputeMass

def StartSession(model_dir=None, event_parity=None):
    """
    Returns: (session, training_params : dict)
    """
    loc = locals()
    invalid_args = [key for key in loc if loc[key] is None]
    if invalid_args:
        raise ValueError(f'Arguments {invalid_args} have illegal value `None`')

    if event_parity not in ['even', 'odd']:
        raise ValueError(f'Illegal event parity `{event_parity}`, only `even` or `odd` are allowed.')

    session = ort.InferenceSession(os.path.join(model_dir, f'{model_dir}_{event_parity}.onnx'))
    params = {}
    with open(os.path.join(model_dir, f'params_{model_dir}_{event_parity}.yaml'), 'r') as cfg:
        params = yaml.safe_load(cfg)

    return session, params

def main():
    # evaluation params
    model_dir = 'predict_quantiles3D_DL_woMP550'
    channel = 'DL'
    masspoint = 550
    suffix = '2B2Vto2B2JLNu' if channel == 'SL' else '2B2Vto2B2L2Nu'
    file = f'../train_data/Run3_2022/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    sample = [t for t in file.split('/') if suffix in t][0]
    plot_name = f'plot_{channel}_{sample}'

    session_even, params_even = StartSession(model_dir=model_dir, event_parity='even')
    session_odd, params_odd = StartSession(model_dir=model_dir, event_parity='odd')

    input_name_even = session_even.get_inputs()[0].name
    output_names_even = [out.name for out in session_even.get_outputs()]
    input_name_odd = session_odd.get_inputs()[0].name
    output_names_odd = [out.name for out in session_odd.get_outputs()]

    assert input_name_even == input_name_odd and output_names_even == output_names_odd, 'Mismatch between even/odd input/output names'
    input_name = input_name_even
    output_names = output_names_even

    dataloader = Dataloader('../experimental/dataloader_config.yaml')
    dataloader.Load(file)
    X_even, input_names, y_even, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)
    X_odd, _, y_odd, _ = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1)

    assert params_even['standardize'] == params_odd['standardize'], '`standardize` value mismatch between even/odd config set'
    standardize = params_even['standardize']
    if standardize:
        input_means_even = params_even['input_train_means']
        input_scales_even = params_even['input_train_scales']
        X_even -= input_means_even
        X_even /= input_scales_even

        input_means_odd = params_odd['input_train_means']
        input_scales_odd = params_odd['input_train_scales']
        X_odd -= input_means_odd
        X_odd /= input_scales_odd

    # apply model trained on even to odd and vice versa
    outputs_odd = session_even.run(output_names, {input_name: X_odd.astype(np.float32)})
    outputs_even = session_odd.run(output_names, {input_name: X_even.astype(np.float32)})
    central_odd = np.array([out[:, 1] for out in outputs_odd[:-1]]).T
    central_even = np.array([out[:, 1] for out in outputs_even[:-1]]).T
    
    target_scales_odd = np.array(params_odd['target_train_scales'])
    target_means_odd = np.array(params_odd['target_train_means'])
    target_scales_even = np.array(params_even['target_train_scales'])
    target_means_even = np.array(params_even['target_train_means'])

    central_odd *= target_scales_odd
    central_odd += target_means_odd
    central_even *= target_scales_even
    central_even += target_means_even
    central = np.concatenate((central_even, central_odd))

    mass = ComputeMass(central)
    PlotHist(data=mass, 
             bins=np.linspace(0, 1500, 100),
             title="Predicted X->HH mass",
             ylabel='Count',
             xlabel='Mass, [GeV]',
             plotting_dir=model_dir,
             file_name=plot_name,
             peak=True,
             width=True)

if __name__ == '__main__':
    main()