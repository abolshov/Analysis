import numpy as np
import vector
import pandas as pd

def ComputeMass(central):
    hvv_en = central[:, 3]
    hbb_en = central[:, -1]

    hvv_p3 = central[:, :3]
    hbb_p3 = central[:, 4:-1]

    x_en = hvv_en + hbb_en
    x_p3 = hvv_p3 + hbb_p3
    x_mass_sqr = np.square(x_en) - np.sum(np.square(x_p3), axis=1)
    neg_mass = x_mass_sqr <= 0.0
    x_mass = np.sqrt(np.abs(x_mass_sqr))
    x_mass = np.where(neg_mass, -1.0, x_mass)
    return x_mass

def PreparePredictions(model, training_params, X_test):
    ys_pred = None
    if training_params['standardize']:
        input_means = training_params['input_train_means']
        input_scales = training_params['input_train_scales']
        X_test -= input_means
        X_test /= input_scales

        print(X_test.shape)

        ys_pred = model.predict(X_test)

        if training_params['add_mass_loss']:
            # drop dimension containing concatenation of all elements before it
            ys_pred = ys_pred[:-1]

        # now ys_pred has shape (n_heads, n_events, n_quantiles) 
        ys_pred = np.array(ys_pred)
        target_scales = np.array(training_params['target_train_scales'])
        target_means = np.array(training_params['target_train_means'])

        # broadcast arrays with means from shape (n_heads,) to shape of ys_pred by inserting new axes
        ys_pred *= target_scales[:, None, None]
        ys_pred += target_means[:, None, None]
        # swap heads and events dimension so events is first
        # (n_heads, n_events, n_quantiles) → (n_events, n_heads, n_quantiles)
        ys_pred = ys_pred.transpose(1, 0, 2)
    else:
        ys_pred = model.predict(X)
        if training_params['add_mass_loss']:
            ys_pred = ys_pred[:-1]
        ys_pred = np.array(ys_pred)
        ys_pred = ys_pred.transpose(1, 0, 2)

    return ys_pred

def ComputeVisMass(df, channel):
    b1 = vector.zip({'px': df['centralJet_1_px'],
                     'py': df['centralJet_1_py'],
                     'pz': df['centralJet_1_pz'],
                     'E': df['centralJet_1_E']})

    b2 = vector.zip({'px': df['centralJet_2_px'],
                     'py': df['centralJet_2_py'],
                     'pz': df['centralJet_2_pz'],
                     'E': df['centralJet_2_E']})

    hbb = b1 + b2
    
    lep1 = vector.zip({'px': df['lep1_px'], 'py': df['lep1_py'], 'pz': df['lep1_pz'], 'E': df['lep1_E']})
    met_E = np.sqrt(df['met_px']**2 + df['met_py']**2)
    met = vector.zip({'px': df['met_px'], 'py': df['met_py'], 'pz': 0.0, 'E': met_E})

    hvv = None
    match channel:
        case 'DL':
            lep2 = vector.zip({'px': df['lep2_px'], 'py': df['lep2_py'], 'pz': df['lep2_pz'], 'E': df['lep2_E']})
            hvv = lep1 + lep2 + met
        case 'SL':
            light_scores = df[[f'centralJet_{i + 1}_btagPNetQvG' for i in range(2, 10)]]
            light_scores = light_scores.values
            light_indices = np.argsort(light_scores, axis=1)
            light_jet_1_idx = light_indices[:, 0] + 2
            light_jet_2_idx = light_indices[:, 1] + 2

            momentum_comp = ['px', 'py', 'pz', 'E']
            jet1_dict = {}
            jet2_dict = {}
            for comp in momentum_comp:
                names1 = [f'centralJet_{i}_{comp}' for i in light_jet_1_idx]
                names2 = [f'centralJet_{i}_{comp}' for i in light_jet_2_idx]

                jet1_dict[comp] = df.to_numpy()[np.arange(len(df)), df.columns.get_indexer(names1)]
                jet2_dict[comp] = df.to_numpy()[np.arange(len(df)), df.columns.get_indexer(names2)]

            light_jet_1 = vector.zip(jet1_dict)
            light_jet_2 = vector.zip(jet2_dict)
            
            hvv = light_jet_1 + light_jet_2 + lep1 + met

    return (hbb + hvv).mass