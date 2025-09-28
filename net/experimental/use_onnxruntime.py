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

def RunInference(file=None, model_dir=None):
    # file = f'../train_data/Run3_2022EE/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    tokens = file.split('/')
    sample_type = tokens[-2]
    masspoint = int(sample_type.split('_')[-1])
    sample = sample_type.split('_')[0]
    channel = 'DL' if '2B2Vto2B2L2Nu' in sample else 'SL'
    era = tokens[2]

    print(file)

    output_base_dir = os.path.join(model_dir, era, channel)
    os.makedirs(output_base_dir, exist_ok=True)
    base_plot_name = f'plot_{channel}_{sample}_M_{masspoint}'

    session_even, params_even = StartSession(model_dir=model_dir, event_parity='even')
    session_odd, params_odd = StartSession(model_dir=model_dir, event_parity='odd')

    input_name_even = session_even.get_inputs()[0].name
    output_names_even = [out.name for out in session_even.get_outputs()]
    input_name_odd = session_odd.get_inputs()[0].name
    output_names_odd = [out.name for out in session_odd.get_outputs()]

    assert input_name_even == input_name_odd and output_names_even == output_names_odd, 'Mismatch between even/odd input/output names'
    input_name = input_name_even
    output_names = output_names_even

    dataloader = Dataloader(f'../experimental/dataloader_cfg_{channel}.yaml')
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

    global_corr_mtrx_odd = np.array(params_even['global_corr_mtrx'])
    global_corr_mtrx_even = np.array(params_odd['global_corr_mtrx'])
    # ErrorPropagator should not take momenta of H->bb and H->WW in the con structor
    # move it to Propagate method arguments
    ep_even = ErrorPropagator(global_corr_mtrx_even, central_even[:, :3], central_even[:, 4:7])
    ep_odd = ErrorPropagator(global_corr_mtrx_odd, central_odd[:, :3], central_odd[:, 4:7])
    errors_odd = ep_odd.Propagate(errors_odd)
    errors_even = ep_even.Propagate(errors_even)
    errors = np.concatenate((errors_even, errors_odd))

    mass = ComputeMass(central)
    PlotHist(data=mass, 
             bins=np.linspace(0, masspoint*2, 100),
             title="Predicted X->HH mass",
             ylabel='Count',
             xlabel='Mass, [GeV]',
             plotting_dir=output_base_dir,
             file_name=f'{base_plot_name}_mass',
             peak=True,
             width=True,
             count=True)

    # PlotHist(data=errors, 
    #          bins=np.linspace(0, 250, 50),
    #          title="Predicted $M_X$ errors",
    #          ylabel='Count',
    #          xlabel='$M_X$ error, [GeV]',
    #          plotting_dir=output_base_dir,
    #          file_name=f'{base_plot_name}_error',
    #          peak=True,
    #          width=True,
    #          count=True)

    # compute vis_mass
    vis_mass_even = ComputeVisMass(df_even, channel)
    vis_mass_odd = ComputeVisMass(df_odd, channel)
    vis_mass = np.concatenate((vis_mass_even, vis_mass_odd))

    # write stuff to file
    output_file_name = f'nn_out_{channel}_M{masspoint}.root' if masspoint is not None else f'nn_out_{channel}_{sample}.root'
    with uproot.recreate(os.path.join(output_base_dir, output_file_name)) as output_file:
        data = {'nn_mass': mass,
                'nn_mass_error': errors,
                'vis_mass': vis_mass,   
                'event': event}

        output_file['Events'] = data

def main():
    # # evaluation params
    # model_dir = 'predict_quantiles3D_DL_v8'
    # channel = 'DL'
    # masspoint = 250
    # # masspoint = None # TTbar or nonres
    # suffix = '2B2Vto2B2JLNu' if channel == 'SL' else '2B2Vto2B2L2Nu'
    # file = f'../train_data/Run3_2022EE/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    # # file = '../train_data/Run3_2022/GluGlutoHHto2B2Vto2L2Nu_kl_1p00_kt_1p00_c2_0p00/nano_0.root'
    # # file = 'TTbar.root'
    # sample = [t for t in file.split('/') if suffix in t][0] if masspoint is not None else 'bkg'
    # base_plot_name = f'plot_{channel}_{sample}'

    # session_even, params_even = StartSession(model_dir=model_dir, event_parity='even')
    # session_odd, params_odd = StartSession(model_dir=model_dir, event_parity='odd')

    # input_name_even = session_even.get_inputs()[0].name
    # output_names_even = [out.name for out in session_even.get_outputs()]
    # input_name_odd = session_odd.get_inputs()[0].name
    # output_names_odd = [out.name for out in session_odd.get_outputs()]

    # assert input_name_even == input_name_odd and output_names_even == output_names_odd, 'Mismatch between even/odd input/output names'
    # input_name = input_name_even
    # output_names = output_names_even

    # dataloader = Dataloader(f'../experimental/dataloader_cfg_{channel}.yaml')
    # dataloader.Load(file)
    # X_even, input_names, y_even, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0)
    # X_odd, _, y_odd, _ = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1)

    # # for visible mass computation
    # df_even = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 0, df=True)
    # df_odd = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1, df=True)
    # event_even = df_even['event']
    # event_odd = df_odd['event']
    # event = np.concatenate((event_even, event_odd))

    # assert params_even['standardize'] == params_odd['standardize'], '`standardize` value mismatch between even/odd config set'
    # standardize = params_even['standardize']
    # if standardize:
    #     input_means_even = params_even['input_train_means']
    #     input_scales_even = params_even['input_train_scales']
    #     X_even -= input_means_even
    #     X_even /= input_scales_even

    #     input_means_odd = params_odd['input_train_means']
    #     input_scales_odd = params_odd['input_train_scales']
    #     X_odd -= input_means_odd
    #     X_odd /= input_scales_odd

    # # apply model trained on even to odd and vice versa
    # outputs_odd = session_even.run(output_names, {input_name: X_odd.astype(np.float32)})
    # outputs_even = session_odd.run(output_names, {input_name: X_even.astype(np.float32)})

    # # srip off output of combined layer
    # outputs_odd = np.array([out for out in outputs_odd[:-1]])
    # outputs_even = np.array([out for out in outputs_even[:-1]])
    # outputs_odd = np.transpose(outputs_odd, (1, 0, 2))
    # outputs_even = np.transpose(outputs_even, (1, 0, 2))

    # if standardize:
    #     target_scales_odd = np.array(params_odd['target_train_scales'])
    #     target_means_odd = np.array(params_odd['target_train_means'])
    #     target_scales_even = np.array(params_even['target_train_scales'])
    #     target_means_even = np.array(params_even['target_train_means'])

    #     target_scales_odd = np.repeat(target_scales_odd[:, None], 3, axis=1)
    #     target_means_odd = np.repeat(target_means_odd[:, None], 3, axis=1)
    #     target_scales_even = np.repeat(target_scales_even[:, None], 3, axis=1)
    #     target_means_even = np.repeat(target_means_even[:, None], 3, axis=1)

    #     outputs_odd *= target_scales_odd
    #     outputs_odd += target_means_odd
    #     outputs_even *= target_scales_even
    #     outputs_even += target_means_even

    # central_odd = outputs_odd[:, :, 1]
    # central_even = outputs_even[:, :, 1]
    # central = np.concatenate((central_even, central_odd))

    # errors_odd = outputs_odd[:, :, 2] - outputs_odd[:, :, 0] 
    # errors_even = outputs_even[:, :, 2] - outputs_even[:, :, 0] 

    # assert params_even['use_energy_layer'] == params_odd['use_energy_layer'], '`use_energy_layer` mismatch between even/odd models training configs'
    # use_energy_layer = params_even['use_energy_layer']
    # if use_energy_layer:
    #     # remove energy errors if energy was computed from mass constraint 
    #     errors_odd = np.delete(errors_odd, [3, 7], axis=1)
    #     errors_even = np.delete(errors_even, [3, 7], axis=1)

    # global_corr_mtrx_odd = np.array(params_even['global_corr_mtrx'])
    # global_corr_mtrx_even = np.array(params_odd['global_corr_mtrx'])
    # # ErrorPropagator should not take momenta of H->bb and H->WW in the con structor
    # # move it to Propagate method arguments
    # ep_even = ErrorPropagator(global_corr_mtrx_even, central_even[:, :3], central_even[:, 4:7])
    # ep_odd = ErrorPropagator(global_corr_mtrx_odd, central_odd[:, :3], central_odd[:, 4:7])
    # errors_odd = ep_odd.Propagate(errors_odd)
    # errors_even = ep_even.Propagate(errors_even)
    # errors = np.concatenate((errors_even, errors_odd))

    # mass = ComputeMass(central)
    # PlotHist(data=mass, 
    #          bins=np.linspace(0, masspoint*2, 100),
    #          title="Predicted X->HH mass",
    #          ylabel='Count',
    #          xlabel='Mass, [GeV]',
    #          plotting_dir=model_dir,
    #          file_name=f'{base_plot_name}_mass',
    #          peak=True,
    #          width=True,
    #          count=True)

    # PlotHist(data=errors, 
    #          bins=np.linspace(0, 250, 50),
    #          title="Predicted $M_X$ errors",
    #          ylabel='Count',
    #          xlabel='$M_X$ error, [GeV]',
    #          plotting_dir=model_dir,
    #          file_name=f'{base_plot_name}_error',
    #          peak=True,
    #          width=True,
    #          count=True)

    # # compute vis_mass
    # vis_mass_even = ComputeVisMass(df_even, channel)
    # vis_mass_odd = ComputeVisMass(df_odd, channel)
    # vis_mass = np.concatenate((vis_mass_even, vis_mass_odd))

    # # write stuff to file
    # output_file_name = f'nn_out_{channel}_M{masspoint}.root' if masspoint is not None else f'nn_out_{channel}_{sample}.root'
    # with uproot.recreate(os.path.join(model_dir, output_file_name)) as output_file:
    #     data = {'nn_mass': mass,
    #             'nn_mass_error': errors,
    #             'vis_mass': vis_mass,   
    #             'event': event}

    #     output_file['Events'] = data

    # model_dir = 'predict_quantiles3D_SL_v3'
    model_dir = 'predict_quantiles3D_DL_v8'

    files = []
    input_files = 'dl_train_files.txt'
    # input_files = 'sl_train_files.txt'
    with open(input_files, 'r') as file_cfg:
        files = [line[:-1] for line in file_cfg.readlines() if '2022EE' not in line and 'Graviton' not in line]

    for file in files:
        RunInference(file=file, model_dir=model_dir)

if __name__ == '__main__':
    main()