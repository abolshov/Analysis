import numpy as np
import uproot
import os
import awkward as ak
import numpy as np
import vector

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from MiscUtils import MemoryMonitor

from typing import Tuple
import numpy.typing as npt

class Transformer(nn.Module):
    def __init__(self,
                 *,
                 input_dim: int = 4, 
                 hidden_dim: int = 64, 
                 num_heads: int = 4, 
                 num_layers: int = 2, 
                 output_dim: int = 4) -> None:
        super(Transformer, self).__init__()
        
        # Input projection
        self.input_proj = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.1)
        )
        
        # Transformer encoder layers
        # this layer implements a transformer as laid out in the paper Attention Is All You Need.
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=hidden_dim,
            nhead=num_heads,
            dim_feedforward=hidden_dim * 4,
            batch_first=True,
            dropout=0.1,
            activation='gelu'  # Using GELU as in the paper
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Global pooling and regression
        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, output_dim)
        )
    
    def forward(self, 
                x: torch.Tensor) -> torch.Tensor:
        # x shape: [batch_size, num_part, 4]
        
        # Project input
        x = self.input_proj(x)  # Shape: [batch_size, num_part, hidden_dim]
        
        # Apply transformer
        x = self.transformer(x)  # Shape: [batch_size, num_part, hidden_dim]
        
        # Global pooling (mean over particles)
        x = x.mean(dim=1)  # Shape: [batch_size, hidden_dim]
        
        # regression head
        x = self.regressor(x).squeeze()  # Shape: [batch_size, output_dim]
        return x
    
def prepare_input_vectors(*,
                          input_path: str | os.PathLike,
                          tree_name: str,
                          channel: str,
                          num_jets: int,
                          num_fatjets: int) -> Tuple[npt.NDArray, npt.NDArray]:
    file = uproot.open(input_path)
    tree = file[tree_name]
    branches = tree.arrays(
        [
            "lep1_pt",
            "lep1_eta",
            "lep1_phi",
            "lep1_mass",
            "lep2_pt",
            "lep2_eta",
            "lep2_phi",
            "lep2_mass",
            "PuppiMET_pt",
            "PuppiMET_phi",
            "centralJet_pt",
            "centralJet_eta",
            "centralJet_phi",
            "centralJet_mass",
            "SelectedFatJet_pt",
            "SelectedFatJet_eta",
            "SelectedFatJet_phi",
            "SelectedFatJet_mass",
            # "X_mass",
            "genHbb_pt",
            "genHbb_eta",
            "genHbb_phi",
            "genHbb_mass",
            "genHVV_pt",
            "genHVV_eta",
            "genHVV_phi",
            "genHVV_mass",
        ]
    )

    lepton_selection = (branches["lep1_pt"] > 0.0) & (branches["lep2_pt"] > 0.0)
    branches = branches[lepton_selection]

    # +1 for MET
    tot_num_particles = num_jets + num_fatjets + 1
    if channel == 'SL':
        tot_num_particles += 1
    elif channel == 'DL':
        tot_num_particles += 2
    else:
        raise RuntimeError(f"Illegal channel {channel}.")

    # output shape: [events, tot_num_particles, 4]
    N_events = len(branches)

    # particles layout:
    # jet1|...|jetN|fatjet1|...|fatjetM|lep1|lep2|MET|
    # particle feautre layout: E|px|py|pz
    features = np.empty((N_events, tot_num_particles, 4))
    labels = np.empty((N_events, 1, 4)) # only p4 of X

    jets = ak.pad_none(
        branches[["centralJet_pt",
                  "centralJet_eta",
                  "centralJet_phi",
                  "centralJet_mass"]], 
        num_jets, 
        clip=True
    )
    jets = ak.fill_none(jets, 0.0)
    jets_p4 = vector.zip({'pt': jets['centralJet_pt'], 
                          'eta': jets['centralJet_eta'],
                          'phi': jets['centralJet_phi'],
                          'mass': jets['centralJet_mass']})

    for jet_idx in range(num_jets):
        for comp_idx, comp_name in enumerate(['E', 'px', 'py', 'pz']):
            features[:, jet_idx, comp_idx] = getattr(jets_p4[:, jet_idx], comp_name)

    fatjets = ak.pad_none(
        branches[["SelectedFatJet_pt",
                  "SelectedFatJet_eta",
                  "SelectedFatJet_phi",
                  "SelectedFatJet_mass"]], 
        num_fatjets, 
        clip=True
    )
    fatjets = ak.fill_none(fatjets, 0.0)
    fatjets_p4 = vector.zip({'pt': fatjets['SelectedFatJet_pt'], 
                             'eta': fatjets['SelectedFatJet_eta'],
                             'phi': fatjets['SelectedFatJet_phi'],
                             'mass': fatjets['SelectedFatJet_mass']})

    for fatjet_idx in range(num_fatjets):
        for comp_idx, comp_name in enumerate(['E', 'px', 'py', 'pz']):
            # need to take into account offset: fatjets start at num_jets
            features[:, num_jets + fatjet_idx, comp_idx] = getattr(fatjets_p4[:, fatjet_idx], comp_name)

    lep1_p4 = vector.zip({'pt': branches['lep1_pt'], 
                          'eta': branches['lep1_eta'],
                          'phi': branches['lep1_phi'],
                          'mass': branches['lep1_mass']})
    offset = num_jets + num_fatjets
    for comp_idx, comp_name in enumerate(['E', 'px', 'py', 'pz']):
        # lep1 starts at num_jets + num_fatjets
        features[:, offset, comp_idx] = getattr(lep1_p4, comp_name)

    if channel == 'DL':
        lep2_p4 = vector.zip({'pt': branches['lep2_pt'], 
                              'eta': branches['lep2_eta'],
                              'phi': branches['lep2_phi'],
                              'mass': branches['lep2_mass']})
        offset += 1
        for comp_idx, comp_name in enumerate(['E', 'px', 'py', 'pz']):
            # lep2 starts at num_jets + num_fatjets + 1 (after lep1)
            features[:, offset, comp_idx] = getattr(lep2_p4, comp_name)

    met_p4 = vector.zip({'pt': branches['PuppiMET_pt'], 
                         'eta': 0.0,
                         'phi': branches['PuppiMET_phi'],
                         'mass': 0.0})
    offset += 1
    for comp_idx, comp_name in enumerate(['E', 'px', 'py', 'pz']):
        # met starts at num_jets + num_fatjets + 2
        features[:, offset, comp_idx] = getattr(met_p4, comp_name)

    # now deal with labels - compute X_p4 = Hbb_p4 + Hww_p4
    Hbb_p4 = vector.zip({'pt': branches['genHbb_pt'], 
                         'eta': branches['genHbb_eta'],
                         'phi': branches['genHbb_phi'],
                         'mass': branches['genHbb_mass']})
    
    HVV_p4 = vector.zip({'pt': branches['genHVV_pt'], 
                         'eta': branches['genHVV_eta'],
                         'phi': branches['genHVV_phi'],
                         'mass': branches['genHVV_mass']})

    X_p4 = Hbb_p4 + HVV_p4

    for comp_idx, comp_name in enumerate(['E', 'px', 'py', 'pz']):
        labels[:, 0, comp_idx] = getattr(X_p4, comp_name)

    return features, labels

def main():
    mm = MemoryMonitor()
    mm.print_memory_usage(msg="Training transformer")

    input_path = "/home/artem/Desktop/CMS/data/DeepHME/big_file_DL_M700.root"
    features, labels = prepare_input_vectors(input_path=input_path,
                                             tree_name="Events",
                                             channel='DL',
                                             num_jets=10,
                                             num_fatjets=2)
    mm.print_memory_usage(msg=f"After loading root file, shapes: {features.shape}, {labels.shape}")

if __name__ == "__main__":
    main()