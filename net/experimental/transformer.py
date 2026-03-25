import numpy as np
import uproot
import os
import awkward as ak
import numpy as np
import vector
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from MiscUtils import MemoryMonitor

from typing import Tuple, Dict, List
import numpy.typing as npt

from dataclasses import dataclass
from typing import Optional

@dataclass
class NormalizationParams:
    feature_mean: npt.NDArray
    feature_std: npt.NDArray
    label_mean: npt.NDArray
    label_std: npt.NDArray

class Transformer(nn.Module):
    def __init__(self,
                 *,
                 input_dim: int = 4, 
                 hidden_dim: int = 64, 
                 num_heads: int = 4, 
                 num_layers: int = 2, 
                 output_dim: int = 4,
                 dropout : float = 0.1) -> None:
        super(Transformer, self).__init__()
        
        # Input projection
        self.input_proj = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout)
        )
        
        # Transformer encoder layers
        # this layer implements a transformer as laid out in the paper Attention Is All You Need.
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=hidden_dim,
            nhead=num_heads,
            dim_feedforward=hidden_dim * 4,
            batch_first=True,
            dropout=dropout,
            activation='gelu'  # Using GELU as in the paper
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Global pooling and regression
        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
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
        x = torch.unsqueeze(x, 1)
        return x
    
class DeepHMEDataset(Dataset):
    def __init__(self,
                 *,
                 file_path: str | os.PathLike,
                 tree_name: str,
                 channel: str,
                 num_jets: int,
                 num_fatjets: int,
                 normalize: bool = False,
                 device: torch.device | None = None,
                 norm_params: Optional[NormalizationParams] = None):
        
        features, labels, computed_norm_params = prepare_input_vectors(
            input_path=file_path,
            tree_name=tree_name,
            channel=channel,
            num_jets=num_jets,
            num_fatjets=num_fatjets,
            normalize=normalize,
            norm_params=norm_params
        )
        
        self.features = torch.as_tensor(features, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.float32)
        self.norm_params = computed_norm_params

        if device is not None:
            self.features = self.features.to(device)
            self.labels = self.labels.to(device)
    
    def __len__(self) -> int:
        return self.features.shape[0]
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.features[idx], self.labels[idx]
    
    def denormalize_labels(self, normalized_labels: torch.Tensor) -> torch.Tensor:
        """Convert model predictions back to physical units."""
        if self.norm_params is None:
            return normalized_labels
        mean = torch.as_tensor(self.norm_params.label_mean, dtype=torch.float32)
        std = torch.as_tensor(self.norm_params.label_std, dtype=torch.float32)
        return normalized_labels * std + mean
    
    @property
    def num_particles(self) -> int:
        return self.features.shape[1]
    
    @property
    def feature_dim(self) -> int:
        return self.features.shape[2]

def prepare_input_vectors(*,
                          input_path: str | os.PathLike,
                          tree_name: str,
                          channel: str,
                          num_jets: int,
                          num_fatjets: int,
                          normalize: bool,
                          norm_params: Optional[NormalizationParams] = None) -> Tuple[npt.NDArray, npt.NDArray, Optional[NormalizationParams]]:
    
    def make_p4(pt, eta, phi, mass):
        """Create a 4-vector from pt, eta, phi, mass."""
        return vector.zip({'pt': pt, 'eta': eta, 'phi': phi, 'mass': mass})
    
    def make_p4_from_branches(branches, prefix):
        """Create a 4-vector from branches with a common prefix."""
        return make_p4(
            branches[f'{prefix}_pt'],
            branches[f'{prefix}_eta'],
            branches[f'{prefix}_phi'],
            branches[f'{prefix}_mass']
        )
    
    def extract_cartesian(p4):
        """Extract [E, px, py, pz] as a stacked numpy array."""
        return np.stack([p4.E, p4.px, p4.py, p4.pz], axis=-1)
    
    def process_collection(branches, prefix, num_objects):
        """Process a variable-length collection (jets/fatjets) into fixed-size array."""
        keys = [f'{prefix}_{var}' for var in ('pt', 'eta', 'phi', 'mass')]
        padded = ak.pad_none(branches[keys], num_objects, clip=True)
        padded = ak.fill_none(padded, 0.0)
        p4 = make_p4(padded[keys[0]], padded[keys[1]], padded[keys[2]], padded[keys[3]])
        # Extract all objects at once: shape [N_events, num_objects, 4]
        return np.stack([
            np.asarray(getattr(p4, comp)) for comp in ('E', 'px', 'py', 'pz')
        ], axis=-1)
    
    # Validate channel
    if channel not in ('SL', 'DL'):
        raise RuntimeError(f"Illegal channel {channel}.")
    
    # Build branch list dynamically to load only what's needed
    base_branches = [
        "lep1_pt", "lep1_eta", "lep1_phi", "lep1_mass",
        "PuppiMET_pt", "PuppiMET_phi",
        "centralJet_pt", "centralJet_eta", "centralJet_phi", "centralJet_mass",
        "SelectedFatJet_pt", "SelectedFatJet_eta", "SelectedFatJet_phi", "SelectedFatJet_mass",
        "genHbb_pt", "genHbb_eta", "genHbb_phi", "genHbb_mass",
        "genHVV_pt", "genHVV_eta", "genHVV_phi", "genHVV_mass",
    ]
    if channel == 'DL':
        base_branches.extend(["lep2_pt", "lep2_eta", "lep2_phi", "lep2_mass"])
    
    # Load data
    with uproot.open(input_path) as file:
        branches = file[tree_name].arrays(base_branches)
    
    # Apply selection
    lepton_selection = branches["lep1_pt"] > 0.0
    if channel == 'DL':
        lepton_selection = lepton_selection & (branches["lep2_pt"] > 0.0)
    branches = branches[lepton_selection]
    
    # Calculate dimensions
    num_leptons = 2 if channel == 'DL' else 1
    tot_num_particles = num_jets + num_fatjets + num_leptons + 1  # +1 for MET
    N_events = len(branches)
    
    # Pre-allocate output arrays
    features = np.empty((N_events, tot_num_particles, 4), dtype=np.float32)
    labels = np.empty((N_events, 1, 4), dtype=np.float32)
    
    # Fill features using slicing (avoid per-index loops)
    offset = 0
    
    # Jets
    features[:, offset:offset + num_jets, :] = process_collection(
        branches, 'centralJet', num_jets
    )
    offset += num_jets
    
    # Fat jets
    features[:, offset:offset + num_fatjets, :] = process_collection(
        branches, 'SelectedFatJet', num_fatjets
    )
    offset += num_fatjets
    
    # Lepton 1
    features[:, offset, :] = extract_cartesian(make_p4_from_branches(branches, 'lep1'))
    offset += 1
    
    # Lepton 2 (DL only)
    if channel == 'DL':
        features[:, offset, :] = extract_cartesian(make_p4_from_branches(branches, 'lep2'))
        offset += 1
    
    # MET
    met_p4 = make_p4(branches['PuppiMET_pt'], 0.0, branches['PuppiMET_phi'], 0.0)
    features[:, offset, :] = extract_cartesian(met_p4)
    
    # Labels: X = Hbb + HVV
    X_p4 = make_p4_from_branches(branches, 'genHbb') + make_p4_from_branches(branches, 'genHVV')
    labels[:, 0, :] = extract_cartesian(X_p4)
    
    if normalize:
        n_events, n_particles, n_components = features.shape
        features_flat = features.reshape(n_events, -1)
        labels_flat = labels.reshape(n_events, -1)
        
        if norm_params is None:
            # Training: compute normalization parameters
            feature_mean = features_flat.mean(axis=0)
            feature_std = features_flat.std(axis=0)
            feature_std[feature_std == 0] = 1.0  # Avoid division by zero
            
            label_mean = labels_flat.mean(axis=0)
            label_std = labels_flat.std(axis=0)
            label_std[label_std == 0] = 1.0
            
            norm_params = NormalizationParams(
                feature_mean=feature_mean,
                feature_std=feature_std,
                label_mean=label_mean,
                label_std=label_std
            )
        
        # Apply normalization (works for both train and inference)
        features_flat = (features_flat - norm_params.feature_mean) / norm_params.feature_std
        labels_flat = (labels_flat - norm_params.label_mean) / norm_params.label_std
        
        features = features_flat.reshape(n_events, n_particles, n_components)
        labels = labels_flat.reshape(n_events, 1, n_components)
    
    return features, labels, norm_params if normalize else None

def get_loader(*,
               file_path: str | os.PathLike,
               tree_name: str,
               channel: str,
               num_jets: int,
               num_fatjets: int,
               batch_size: int,
               normalize: bool,
               num_workers: int,
               device: torch.device | None = None,
               norm_params: Optional[NormalizationParams] = None) -> Tuple[DataLoader, Optional[NormalizationParams]]:
    
    ds = DeepHMEDataset(
        file_path=file_path,
        tree_name=tree_name,
        channel=channel,
        num_jets=num_jets,
        num_fatjets=num_fatjets,
        normalize=normalize,
        norm_params=norm_params,
        device=device
    )

    loader = DataLoader(
        ds,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True
    )

    return loader, ds.norm_params

def plot_training_history(*,
                          history: Dict[str, npt.NDArray], 
                          plot_name: str,
                          plotting_dir: str | os.PathLike,
                          metrics: List[str] = ['loss']):
    """
    Plot training history for neural networks
    """
    # _, axes = plt.subplots(1, len(metrics), figsize=(15, 5))
    _, axes = plt.subplots(1, len(metrics))
    if len(metrics) == 1:
        axes = [axes]

    try: 
        # this is the syntax for keras 
        hist = history.history
    except: 
        hist = history
    for ax, metric in zip(axes, metrics):
        ax.plot(hist[metric], label=f'Training {metric}')
        ax.plot(hist[f'val_{metric}'], label=f'Validation {metric}')
        ax.set_xlabel('Epoch')
        ax.set_ylabel(metric.capitalize())
        ax.legend()
        ax.grid(True)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plotting_dir, f"{plot_name}.pdf"), format="pdf")
    plt.clf()
    plt.close()

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    mm = MemoryMonitor()
    mm.print_memory_usage(msg=f"Training transformer on {device}")

    num_jets = 6
    num_fatjets = 1
    batch_size = 128
    num_workers = 2
    num_epochs = 10

    train_file_path = "/home/artem/Desktop/CMS/data/DeepHME/big_file_DL_M700.root"
    train_loader, train_norm = get_loader(
        file_path=train_file_path,
        tree_name="Events",
        channel="DL",
        num_jets=num_jets,
        num_fatjets=num_fatjets,
        batch_size=batch_size,
        normalize=True,
        norm_params=None,
        num_workers=num_workers,
    )
    mm.print_memory_usage(msg=f"After creating train loader with {len(train_loader)} batches")

    val_file_path = "/home/artem/Desktop/CMS/data/DeepHME/Run3_2023BPix/XtoYHto2B2Wto2B2L2Nu_MX_700_MY_125/anaTuple_0.root"
    val_loader, _ = get_loader(
        file_path=val_file_path,
        tree_name="Events",
        channel="DL",
        num_jets=num_jets,
        num_fatjets=num_fatjets,
        batch_size=batch_size,
        normalize=True,
        norm_params=train_norm,
        num_workers=num_workers
    )
    mm.print_memory_usage(msg=f"After creating validation loader with {len(val_loader)} batches")

    model = Transformer().to(device)
    criterion = nn.MSELoss()
    optimizer = optim.AdamW(model.parameters(), 
                            lr=0.001,
                            weight_decay=0.01)

    # num_steps = len(train_loader)*num_epochs
    # cosine = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=10, T_mult=2)
    # n_warmup_steps = 2000
    # warmup = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lambda step: min((step + 1)/n_warmup_steps, 1.0))
    # scheduler = torch.optim.lr_scheduler.SequentialLR(optimizer, 
    #                                                   schedulers=[warmup, cosine], 
    #                                                   milestones=[n_warmup_steps])

    history = {
        'loss': [],
        'val_loss': []
    }

    # Training loop
    for epoch in range(num_epochs):
        model.train()
        train_loss = 0
        
        for features, labels in train_loader:
            features, labels = features.to(device), labels.to(device)
            
            optimizer.zero_grad()
            outputs = model(features)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            # scheduler.step()
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0
        
        with torch.no_grad():
            for features, labels in val_loader:
                features, labels = features.to(device), labels.to(device)
                outputs = model(features)
                loss = criterion(outputs, labels)
                val_loss += loss.item()
        
        # Calculate epoch metrics
        train_loss /= len(train_loader)
        val_loss /= len(val_loader)
        
        # Append to history
        history['loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        
        print(f'Epoch {epoch + 1}/{num_epochs}:')
        print(f'\ttrain_loss: {train_loss:.4f}')
        print(f'\tval_loss: {val_loss:.4f}')

    model_dir = os.path.join("/home/artem/Desktop/CMS/models/Transformer/")
    os.makedirs(model_dir, exist_ok=True)
    plot_training_history(history=history, 
                          plot_name="loss",
                          plotting_dir=model_dir)    

    torch.save(model, os.path.join(model_dir, "transformer.pth"))

if __name__ == "__main__":
    main()