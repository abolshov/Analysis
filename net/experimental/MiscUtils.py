import numpy as np
import vector
import psutil
import os
import uproot
import pathlib
import math
import awkward as ak
import re

from typing import List, Dict, Union, Tuple
from numpy.typing import NDArray

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

class MemoryMonitor:
    def __init__(self) -> None:
        self.process = psutil.Process(os.getpid())

    def print_memory_usage(self, 
                           *,
                           msg: str = None) -> None:
        memory_bytes = self.process.memory_info().rss
        memory_mb = memory_bytes / (1024 * 1024)
        if msg:
            print(msg)
        print(f"Memory usage {memory_mb:.2f} MB")

def to_numpy(*,
             ak_array: ak.Array, 
             dtype=np.float32) -> np.ndarray:
    """
    Pre-allocate numpy array for better memory efficiency
    """
    fields = ak.fields(ak_array)
    M = len(ak_array)
    N = len(fields)
    
    result = np.empty((M, N), dtype=dtype)
    
    for i, field in enumerate(fields):
        result[:, i] = ak.to_numpy(ak_array[field])
    
    return result

def load_file(*,
              tree_name: str,
              file_path: str | os.PathLike | pathlib.Path,
              list_of_branches: List[str],
              convert_to_numpy: bool):

    file = uproot.open(file_path)
    tree = file[tree_name]
    branches = tree.arrays(list_of_branches)

    if convert_to_numpy:
        return to_numpy(ak_array=branches)
    
    return branches

def nearest_pow2(n: int) -> int:
    if n <= 0:
        return 1
    
    p_high = math.ceil(math.log2(n))
    v_high = 2**p_high
    return v_high

def map_input_files(*,
                    directory: str | os.PathLike | pathlib.Path,
                    event_file_pattern: re.Pattern[str],
                    weight_file_pattern: re.Pattern[str]) -> Dict[int, Union[str, os.PathLike, pathlib.Path]]:
    
    event_files = {}
    weight_files = {}
    for f in os.listdir(directory):
        weight_match = weight_file_pattern.search(f)
        event_match = event_file_pattern.search(f)
        abs_path = os.path.abspath(os.path.join(directory, f))
        if weight_match:
            parity = int(weight_match.group(1))
            weight_files[parity] = abs_path
        elif event_match:
            parity = int(event_match.group(1))
            event_files[parity] = abs_path
        else:
            continue

    parity_file_map = { p: (event_files[p], weight_files[p]) for p in event_files.keys()}
    return parity_file_map

def clean_extreme_values(*,
                         data: NDArray[np.floating],
                         pos_thrsh: float = 2500.0,
                         neg_thrsh: float = -1000.0) -> Tuple[NDArray[np.bool_], NDArray[np.floating]]:
    mask = ~np.any((data > pos_thrsh) | (data < neg_thrsh), axis=1)
    return mask, data[mask]