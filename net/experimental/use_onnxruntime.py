import onnxruntime as ort
import numpy as np
import yaml
import os
import uproot

from Dataloader import Dataloader
from PlotUtils import PlotHist
from PrepUtils import ComputeMass, ComputeVisMass
from ErrorProp import ErrorPropagator

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
    model_dir = 'predict_quantiles3D_DL_v8'
    channel = 'DL'
    masspoint = 1000
    suffix = '2B2Vto2B2JLNu' if channel == 'SL' else '2B2Vto2B2L2Nu'
    file = f'../train_data/Run3_2022/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    sample = [t for t in file.split('/') if suffix in t][0]
    base_plot_name = f'plot_{channel}_{sample}'

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

    # for visible mass computation
    df_even = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0, df=True)
    df_odd = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1, df=True)
    event_even = df_even['event']
    event_odd = df_odd['event']
    event = np.concatenate((event_even, event_odd))

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

    # srip off output of combined layer
    outputs_odd = np.array([out for out in outputs_odd[:-1]])
    outputs_even = np.array([out for out in outputs_even[:-1]])
    outputs_odd = np.transpose(outputs_odd, (1, 0, 2))
    outputs_even = np.transpose(outputs_even, (1, 0, 2))

    if standardize:
        target_scales_odd = np.array(params_odd['target_train_scales'])
        target_means_odd = np.array(params_odd['target_train_means'])
        target_scales_even = np.array(params_even['target_train_scales'])
        target_means_even = np.array(params_even['target_train_means'])

        target_scales_odd = np.repeat(target_scales_odd[:, None], 3, axis=1)
        target_means_odd = np.repeat(target_means_odd[:, None], 3, axis=1)
        target_scales_even = np.repeat(target_scales_even[:, None], 3, axis=1)
        target_means_even = np.repeat(target_means_even[:, None], 3, axis=1)

        outputs_odd *= target_scales_odd
        outputs_odd += target_means_odd
        outputs_even *= target_scales_even
        outputs_even += target_means_even

    central_odd = outputs_odd[:, :, 1]
    central_even = outputs_even[:, :, 1]
    central = np.concatenate((central_even, central_odd))

    errors_odd = outputs_odd[:, :, 2] - outputs_odd[:, :, 0] 
    errors_even = outputs_even[:, :, 2] - outputs_even[:, :, 0] 

    assert params_even['use_energy_layer'] == params_odd['use_energy_layer'], '`use_energy_layer` mismatch between even/odd models training configs'
    use_energy_layer = params_even['use_energy_layer']
    if use_energy_layer:
        # remove energy errors if energy was computed from mass constraint 
        errors_odd = np.delete(errors_odd, [3, 7], axis=1)
        errors_even = np.delete(errors_even, [3, 7], axis=1)

    global_corr_mtrx_even = np.array(params_even['global_corr_mtrx'])
    global_corr_mtrx_odd = np.array(params_odd['global_corr_mtrx'])
    # ErrorPropagator should not take momenta of H->bb and H->WW in the con structor
    # move it to Propagate method arguments
    ep_even = ErrorPropagator(global_corr_mtrx_even, central_even[:, :3], central_even[:, 4:7])
    ep_odd = ErrorPropagator(global_corr_mtrx_odd, central_odd[:, :3], central_odd[:, 4:7])
    errors_odd = ep_odd.Propagate(errors_odd)
    errors_even = ep_even.Propagate(errors_even)
    errors = np.concatenate((errors_even, errors_odd))

    mass = ComputeMass(central)
    PlotHist(data=mass, 
             bins=np.linspace(0, 2500, 100),
             title="Predicted X->HH mass",
             ylabel='Count',
             xlabel='Mass, [GeV]',
             plotting_dir=model_dir,
             file_name=f'{base_plot_name}_mass',
             peak=True,
             width=True,
             count=True)

    PlotHist(data=errors, 
             bins=np.linspace(0, 500, 50),
             title="Predicted $M_X$ errors",
             ylabel='Count',
             xlabel='$M_X$ error, [GeV]',
             plotting_dir=model_dir,
             file_name=f'{base_plot_name}_error',
             peak=True,
             width=True,
             count=True)

    # compute vis_mass
    vis_mass_even = ComputeVisMass(df_even, channel)
    vis_mass_odd = ComputeVisMass(df_odd, channel)
    vis_mass = np.concatenate((vis_mass_even, vis_mass_odd))

    # write stuff to file
    with uproot.recreate(os.path.join(model_dir, f'nn_out_{channel}_M{masspoint}.root')) as output_file:
        data = {'nn_mass': mass,
                'nn_mass_error': errors,
                'vis_mass': vis_mass,   
                'event': event}

        output_file['Events'] = data

if __name__ == '__main__':
    main()