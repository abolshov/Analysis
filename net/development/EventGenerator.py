import numpy as np
import awkward as ak
import pathlib
import uproot

from typing import Any

# from .OutputFormatter import OutputFormatter
# from .Augmentor import Augmentor
# from .EventSelection import Selector

class EventGenerator:
    """
    Class implementing a generator of events for a given mass point.
    It iterates over list of files corresponding to selected masspoint and
    yields selected events (selection logic will be done in a separate class).
    Additionally will be able to pass selected events to augmentation algorithm
    to increase the dataset size (augmentaton will be habled by a separate class).
    Will be able to yield outputs according to format provided by user-defined
    formatter class (to be implemented).
    """
    def __init__(self,
                 masspoint: int,
                 file_list: list[pathlib.Path],
                 branches_to_load: list[str],
                 features: list[str],
                 targets: list[str],
                 step_size: str | int,
                 yield_size: int=1, 
                 selector: Any=None):
                #  output_formatter: OutputFormatter=None,
                #  augmentor: Augmentor=None):
        
        self.masspoint = masspoint
        self.file_list = file_list
        self.branches_to_load = branches_to_load
        self.features = features
        self.targets = targets
        self.step_size = step_size
        self.yield_size = yield_size

        # these need to be implemented
        # selector object is the same for all mass points, should be created in the flow before mass generator and passed to it
        # it only reduces tree along the `rows` dimension (events), but not `columns` dimension (event objects)
        self.selector = selector 

        # self.augmentor = augmentor

        # formatter object is also the same for all mass points, should be created in the flow before mass generator and passed to it
        # it reduces tree along the `columns` dimension (event objects) and transforms and computes new variables and removes those that are not needed
        # it also changes makes shape of the output of what will be yielded by generator equal what it is supposed to be 
        # self.output_formatter = output_formatter
    
    def get_generator(self):
        """
        Method provides generator (required by tensorflow)
        """

        def generator():
            for chunk in uproot.iterate(
                self.file_list,
                "Events",
                self.branches_to_load,
                self.step_size,
            ):
                # select events
                if self.selector:
                    chunk = self.selector.apply_selections(chunk)

                # prepare features, labels according to desired format
                features, labels = ...

                # perform augmentation if needed

                # yield results
                yield features, labels
        
        return generator()
