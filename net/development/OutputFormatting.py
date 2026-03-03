import numpy as np
import awkward as ak
import re

from typing import List, Dict, Tuple, List
from Utils import stack_columns, make_p4

class Tabulizer:
    """
    Class for converting ak.Array to tuple of numpy arrays.
    Args:
        particle_paddings: provides padding length per particle.
    """
    def __init__(self,
                 feature_branches: List[str],
                 target_branches: List[str],
                 particle_paddings: Dict[str, int],
                 input_kinematic_variables: List[str],
                 output_kinematic_variables: List[str],
                 met4D: bool) -> None:
        
        if set(input_kinematic_variables) != {'pt', 'eta', 'phi', 'mass'} or set(input_kinematic_variables) != {'px', 'py', 'pz', 'E'}:
            raise ValueError(f'Illegal format of kinematics {input_kinematic_variables} for input branches.')
        
        self.feature_branches = feature_branches
        self.target_branches = target_branches
        self.particle_paddings = particle_paddings
        self.input_kinematic_variables = input_kinematic_variables
        self.output_kinematic_variables = output_kinematic_variables
        self.input_kinematic_pattern = re.compile(rf'^(?P<particle>.+)_(?P<variable>{'|'.join(map(re.escape, input_kinematic_variables))})$')
        self.met4D = met4D

    def format(self,
               chunk: ak.Array) -> Tuple[np.ndarray, np.ndarray]:
        
        chunk_dict = self._pad_and_flatten_chunk(chunk)
        features = stack_columns(arr_dict=chunk_dict, col_names=self.feature_branches)
        targets = stack_columns(arr_dict=chunk_dict, col_names=self.target_branches)
        return features, targets
    
    def _pad_and_flatten_chunk(self, 
                               chunk: ak.Array) -> Dict[str, np.ndarray]:
        res = {}

        # find 'particles' (e.g. lep1, met, centralJet, etc)
        particles = []
        for name in chunk.fields:
            match = self.input_kinematic_pattern.search(name)
            if match:
                particle = match.group('particle')
                particles.append(particle)
        
        # all 'particle' branches follow this pattern
        particle_pattern = re.compile(rf'({'|'.join(particles)})_*')

        # free branches do not have particle pattern in them
        free_branches = [name for name in chunk.fields if not particle_pattern.search(name)]

        for fb in free_branches:
            data = chunk[fb]
            if data.ndim != 1:
                raise RuntimeError(f'Branch {fb} has illegal dimension {data.ndim}.')
            res[name] = ak.to_numpy(data)

        for p in particles:
            if self.particle_paddings is None or self.particle_paddings.get(p) is None:
                raise RuntimeError(f'Cannot consistently pad branch `{name}` across all chunks because padding length was not provided.')
            max_len = self.particle_paddings[p]

            # get branches such as btag, etc and save them
            particle_feature_branches = [pb for pb in chunk.fields if particle_pattern.search(pb) and not self.input_kinematic_pattern.search(pb)]
            for pb in particle_feature_branches:
                tokens = pb.split('_')
                variable = tokens[-1]
                prefix = '_'.join(tokens[:-1])
                data = chunk[pb]
                if max_len > 1:
                    padded = ak.pad_none(data, max_len, clip=True)
                    filled = ak.fill_none(padded, 0)
                    for i in range(max_len):
                        res[f'{prefix}_{i+1}_{variable}'] = ak.to_numpy(filled[:, i])
                else:
                    res[f'{prefix}_{variable}'] = ak.to_numpy(data)
            
            # now make p4 to save momentum data
            p4 = make_p4(data=chunk, 
                         prefix=p, 
                         coordinates=self.input_kinematic_variables,
                         default_init_missing=self.met4D)
            
            for var in self.output_kinematic_variables:
                if not hasattr(p4, var):
                    # that will save me for the case of MET!
                    # make_p4 will create met as 2D vector and it does not have 4D attributes (pz or Es)
                    continue
                
                arr = getattr(p4, var)
                if max_len > 1:
                    padded_arr = ak.pad_none(arr, max_len, clip=True)
                    filled_arr = ak.fill_none(padded_arr, 0)
                    for i in range(max_len):
                        res[f'{p}_{i+1}_{var}'] = ak.to_numpy(filled_arr[:, i])
                else:
                    res[f'{p}_{var}'] = ak.to_numpy(arr)

        return res
