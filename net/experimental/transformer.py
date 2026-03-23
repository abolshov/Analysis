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
        # x shape: [batch_size, 4, input_dim]
        
        # Project input
        x = self.input_proj(x)  # Shape: [batch_size, 4, hidden_dim]
        
        # Apply transformer
        x = self.transformer(x)  # Shape: [batch_size, 4, hidden_dim]
        
        # Global pooling (mean over particles)
        x = x.mean(dim=1)  # Shape: [batch_size, hidden_dim]
        
        # regression head
        x = self.regressor(x).squeeze()  # Shape: [batch_size, output_dim]
        return x

def prepare_input_vectors(*,
                          input_path: str | os.PathLike,
                          tree_name: str,
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
            "X_mass",
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

    print(f"\tlen(branches)={len(branches)}")
    lepton_selection = (branches["lep1_pt"] > 0.0) & (branches["lep2_pt"] > 0.0)
    print(f"\tnum_valid_events={np.sum(lepton_selection)}")
    branches = branches[lepton_selection]
    print(f"\tlen(branches)={len(branches)}")

    # todo
    features = np.empty()
    labels = np.empty()

    pass

def main():
    mm = MemoryMonitor()
    mm.print_memory_usage(msg="Training transformer")

    input_path = "/home/artem/Desktop/CMS/data/DeepHME/big_file_DL_M700.root"
    _ = prepare_input_vectors(input_path=input_path,
                              tree_name="Events",
                              num_jets=10,
                              num_fatjets=2)

    mm.print_memory_usage(msg="After loading root file")

if __name__ == "__main__":
    main()