import numpy as np

ground_truth_map = {"genHbb_E": "E(H->bb)",
                    "genHbb_px": "px(H->bb)",
                    "genHbb_py": "py(H->bb)",
                    "genHbb_pz": "pz(H->bb)",
                    "genHbb_pt": "Pt(H->bb)",
                    "genHbb_eta": "Eta(H->bb)",
                    "genHbb_phi": "Phi(H->bb)",
                    "genHVV_E": "E(H->VV)",
                    "genHVV_px": "px(H->VV)",
                    "genHVV_py": "py(H->VV)",
                    "genHVV_pz": "pz(H->VV)",
                    "genHVV_pt": "Pt(H->VV)",
                    "genHVV_eta": "Eta(H->VV)",
                    "genHVV_phi": "Phi(H->VV)" }


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