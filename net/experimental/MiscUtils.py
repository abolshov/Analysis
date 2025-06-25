import numpy as np

target_names = {"genHbb_E": "E(H->bb)",
                "genHbb_px": "px(H->bb)",
                "genHbb_py": "py(H->bb)",
                "genHbb_pz": "pz(H->bb)",
                "genHVV_E": "E(H->VV)",
                "genHVV_px": "px(H->VV)",
                "genHVV_py": "py(H->VV)",
                "genHVV_pz": "pz(H->VV)" }


def PredWidth(pred_mass):
    q_84 = np.quantile(pred_mass, 0.84)
    q_16 = np.quantile(pred_mass, 0.16)
    width = q_84 - q_16
    return width 


def PredPeak(pred_mass):
    counts = np.bincount(pred_mass)
    peak = np.argmax(counts)
    return peak 


def Scheduler(epoch, lr):
    if epoch < 30:
        return lr
    else:
        if epoch % 2 == 0:
            return 0.9*lr
        return lr