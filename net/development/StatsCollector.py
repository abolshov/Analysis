import numpy as np
import awkward as ak
import pathlib
import uproot

class StatsCollector:
    def __init__(self,
                 file_list: list[pathlib.Path],
                 branches_to_load: list[str],
                 features: list[str],
                 targets: list[str],
                 step_size: str | int,
                 verbose=False):
        
        self.verbose = verbose
        if self.verbose:
            print(f'Initializing StatsCollector with branches_to_load={branches_to_load}, features={features}, targets={targets}, file_lsit={file_list}')

        self.file_list = file_list
        self.branches_to_load = branches_to_load
        self.features = features
        self.targets = targets
        self.step_size = step_size

        self.n = 0
        self.sum_x_features = np.zeros(len(features))
        self.sum_x2_features = np.zeros(len(features))
        self.sum_x_targets = np.zeros(len(targets))
        self.sum_x2_targets = np.zeros(len(targets))

    def collect(self):
        for file_path in self.file_list:
            if self.verbose:
                print(f'Processing {file_path}')

            for chunk in uproot.iterate(file_path,
                                        "Events", 
                                        self.branches_to_load, 
                                        step_size=self.step_size, 
                                        library="np"):
                # select events

                
                # calc stats
                self.n += ...
                self.sum_x_features += ...
                self.sum_x2_features += ...

    def get_stats(self):
        pass
