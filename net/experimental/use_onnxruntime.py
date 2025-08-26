import onnxruntime as ort
import numpy as np
import yaml
import os

from Dataloader import Dataloader
from PlotUtils import PlotHist
from PrepUtils import ComputeMass

def main():
    session = ort.InferenceSession('predict_quantiles3D_DL_v4/model.onnx')
    training_params = {}
    with open('predict_quantiles3D_DL_v4/params_model.yaml', 'r') as train_cfg_file:
        training_params = yaml.safe_load(train_cfg_file)

    input_name = session.get_inputs()[0].name
    output_names = [out.name for out in session.get_outputs()]

    channel = 'DL'
    masspoint = 800
    suffix = '2B2Vto2B2JLNu' if channel == 'SL' else '2B2Vto2B2L2Nu'
    file = f'../train_data/Run3_2022/GluGlutoRadiontoHHto{suffix}_M_{masspoint}/nano_0.root'
    dataloader = Dataloader('../experimental/dataloader_config.yaml')
    dataloader.Load(file)
    X, input_names, y, target_names = dataloader.Get(lambda df, mod, parity: df['event'] % mod == parity, 2, 1)

    if training_params['standardize']:
        input_means = training_params['input_train_means']
        input_scales = training_params['input_train_scales']
        X -= input_means
        X /= input_scales

    outputs = session.run(output_names, {input_name: X.astype(np.float32)})
    central = np.array([out[:, 1] for out in outputs[:-1]]).T
    
    target_scales = np.array(training_params['target_train_scales'])
    target_means = np.array(training_params['target_train_means'])

    central *= target_scales
    central += target_means

    mass = ComputeMass(central)
    PlotHist(data=mass, 
             bins=np.linspace(0, 2500, 100),
             title="Predicted X->HH mass",
             ylabel='Count',
             xlabel='Mass, [GeV]',
             plotting_dir=os.path.join(training_params['model_dir'], 'plots'),
             file_name='mass_onnx',
             peak=True,
             width=True)



if __name__ == '__main__':
    main()