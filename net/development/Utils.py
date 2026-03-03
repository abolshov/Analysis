import numpy as np
import awkward as ak
import vector
import re
from typing import Dict, Union, List, Literal

def pad_chunk(chunk: ak.Array,
              branch_paddings: Dict[str, int],
              to_cartesian: bool) -> Dict[str, np.ndarray]:
    """
    Takes an Awkward Array chunk and converts branches into a dictionary 
    of flat NumPy arrays, padding jagged branches with 0.
    Also takes an optional argumnet `branch_padings` specifying maximum padding length of branches.
    """
    pat = re.compile(r'^(?P<particle>.+)_(?P<variable>pt|eta|phi|mass)$')
    particles = {}

    res = {}
    
    for name in chunk.fields:
        data = chunk[name]
        # Check if the data is jagged (multi-dimensional)
        if data.ndim > 1:
            # Calculate the maximum length in this specific chunk
            if branch_paddings is None or branch_paddings.get(name) is None:
                raise RuntimeError(f'Cannot consistently pad branch `{name}` across all chunks because padding length was not provided.')
            
            max_len = branch_paddings[name]
            
            # Pad to the max length and fill empty slots with 0
            padded = ak.pad_none(data, max_len, clip=True)
            filled = ak.fill_none(padded, 0)
            
            # Convert the regular Awkward structure to a NumPy ndarray
            for i in range(max_len):
                res[f'{name}_{i+1}'] = ak.to_numpy(filled[:, i])
        else:
            # For 1D (scalar) branches, convert directly
            res[name] = ak.to_numpy(data)
            
    return res

def mask_indices(indices: Union[ak.Array, np.ndarray, List[int]], 
                 array_to_slice: ak.Array) -> ak.Array:
    """
    https://stackoverflow.com/questions/78657390/accessing-elements-of-an-awkward-array-that-are-not-a-passed-in-index
    Generates a boolean mask for an Awkward Array based on a list of target indices.
    
    This function compares the local indices of `array_to_slice` against a provided 
    set of `indices`. It is particularly useful for ragged arrays where you need 
    to filter elements by their position within a sub-list.

    Args:
        indices: The indices to be masked. Can be a flat list, a NumPy array, 
            or an Awkward Array.
        array_to_slice: The target Awkward Array from which to generate the mask.

    Returns:
        A boolean Awkward Array with the same structure as `array_to_slice`, 
        where True indicates the index was present in the `indices` input.

    Example:
        >>> data = ak.Array([[10, 20, 30], [40, 50]])
        >>> idx = [[0, 2]]
        >>> mask_indices(idx, data)
        <Array [[True, False, True], [True, False]] >
    """
    # Generate local indices (0, 1, 2...) for each sub-list
    # Create a Cartesian product to compare every local index with every target index
    whole_set, in_set = ak.unzip(
        ak.cartesian([ak.local_index(array_to_slice), indices], nested=True)
    )
    
    # Return True if any comparison in the nested dimension matches
    return ak.any(whole_set == in_set, axis=-1)

def make_p4(data: ak.Array,
            prefix: str,
            coordinates: List[str],
            default_init_missing: bool) -> ak.Array:
    """
    Construct a 4-momentum array using the `vector` library from an Awkward Array.

    This function reads kinematic components from `data` using a common prefix
    and assembles them into a vector-compatible 4-momentum object.

    Parameters
    ----------
    data : ak.Array
        Awkward Array containing kinematic fields. Expected field names depend
        on the chosen coordinate system and are constructed as
        ``f"{prefix}_<component>"``.
    prefix : str
        Prefix used to locate kinematic fields in `data`
        (e.g. ``"jet"``, ``"muon"``).
    coordinates : {['px', 'py', 'pz', 'E'], ['pt', 'eta', 'phi', 'mass']}
        Coordinate system used to build the 4-momentum:
    default_init_missing : bool
        Flag for default-initializing missing branches to 0.
        Useful for MET: MET only has 2 features (e.g. pt and phi), but if 4D vector of MET is desired, True should be passed.

    Returns
    -------
    ak.Array
        An Awkward Array with vector behaviors enabled, representing
        Lorentz 4-vectors compatible with the `vector` library.
    """
    if default_init_missing:
        p4_dict = vector.zip({var: data[f"{prefix}_{var}"] if f"{prefix}_{var}" in data.fields else 0 for var in coordinates})
    else:
        p4_dict = {}
        for var in coordinates:
            if var not in data.fields:
                continue
            p4_dict[var] = data[f"{prefix}_{var}"]
    p4 = vector.zip(p4_dict)
    return p4

def dr2(data: ak.Array, 
        prefix1: str, 
        prefix2: str):
    """
    Calculates the squared angular separation (Delta R^2) between two particles.
    Using R^2 is more memory-efficient as it avoids the expensive sqrt() operation.
    """
    deta = data[f'{prefix1}_eta'] - data[f'{prefix2}_eta']
    dphi = ak.fill_none(data[f'{prefix1}_phi'], 0).delta_phi(ak.fill_none(data[f'{prefix2}_phi'], 0))
    return deta**2 + dphi**2

def stack_columns(arr_dict: Dict[str, int],
                  col_names: List[str]) -> np.ndarray:
        arrays_to_stack = []
        for name in col_names:
            arr = arr_dict[name]
            
            # If the array is 2D (n_events, max_objects), 
            # we flatten the objects into the feature width.
            if arr.ndim == 2:
                # Reshape from (N, M) to (N, M) - already compatible for hstack
                # If it was a jagged-padded array, it's now N rows of M features
                arrays_to_stack.append(arr)
            else:
                # If it's 1D (n_events,), reshape to (n_events, 1) to stack it
                arrays_to_stack.append(arr.reshape(-1, 1))
        
        # Horizontal stack: glues arrays along the feature axis (axis 1)
        return np.hstack(arrays_to_stack)