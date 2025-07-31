import numpy as np
import vector

ground_truth_map = {"genHbb_E": "E(H->bb)",
                    "genHbb_px": r"$P_x$(H->bb)",
                    "genHbb_py": r"$P_y$(H->bb)",
                    "genHbb_pz": r"$P_z$(H->bb)",
                    "genHbb_pt": r"$P_T$(H->bb)",
                    "genHbb_eta": r"$\eta$(H->bb)",
                    "genHbb_phi": r"$\phi$(H->bb)",
                    "genHVV_E": "E(H->VV)",
                    "genHVV_px": r"$P_x$(H->VV)",
                    "genHVV_py": r"$P_y$(H->VV)",
                    "genHVV_pz": r"$P_z$(H->VV)",
                    "genHVV_pt": r"$P_T$(H->VV)",
                    "genHVV_eta": r"$\eta$(H->VV)",
                    "genHVV_phi": r"$\phi$(H->VV)" }

pretty_vars = {'pt': r'$P_T$',
               'px': r'$P_x$',
               'py': r'$P_y$',
               'pz': r'$P_z$',
               'eta': r'$\eta$',
               'phi': r'$\phi$'}

objects = {'genHbb': 'H->bb',
           'genHVV': 'H->VV'}


def PredWidth(arr):
    q_84 = np.quantile(arr, 0.84)
    q_16 = np.quantile(arr, 0.16)
    width = q_84 - q_16
    return width 


def PredPeak(arr, bins='auto'):
    counts, edges = np.histogram(arr, bins=bins)
    binmax = np.argmax(counts)
    peak = (edges[binmax] + edges[binmax + 1])/2
    return peak 


def Scheduler(epoch, lr):
    if epoch < 30:
        return lr
    else:
        if epoch % 2 == 0:
            return 0.9*lr
        return lr


def ToPtEtaPhiE(data):
    """
    converts PxPyPzE vector to PtEtaPhiE vector
    """

    assert data.shape[-1] == 4, f'Wrong dimension {data.shape[-1]} of the input vector'

    p4 = vector.zip({'px': data[:, 0], 'py': data[:, 1], 'pz': data[:, 2], 'E': data[:, 3]})
    return np.stack([p4.pt.to_numpy(), p4.eta.to_numpy(), p4.phi.to_numpy(), p4.E.to_numpy()], axis=1)