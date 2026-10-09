"""The sequence model, defined and smoke-tested before it is ever fitted (D-026).

A GRU over the H history steps, each step carrying the CGM value and - for the primary
model - the three covariates the protocol allows (basal, bolus, carbs; D-010). The head
predicts the glucose change from the last history slot, so a zero output is persistence:
the model starts at the floor every later number is measured from, and has to earn any
departure from it. Input scaling is a fixed affine map, not a fitted scaler.

Nothing here touches data. Saturday's training script is src/train.py.
"""
from __future__ import annotations

import torch
from torch import nn

CGM_SCALE = (140.0, 60.0)        # (centre, scale) mg/dL - fixed constants, not fitted
COV_SCALE = {"basal": 1.0, "bolus": 1.0, "carbs": 20.0}   # units per feature, fixed


class GlucoseGRU(nn.Module):
    def __init__(self, n_features: int = 4, hidden: int = 64, layers: int = 1, dropout: float = 0.0):
        super().__init__()
        self.n_features = n_features
        self.gru = nn.GRU(n_features, hidden, num_layers=layers, batch_first=True,
                          dropout=dropout if layers > 1 else 0.0)
        self.head = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, 1))
        nn.init.zeros_(self.head[-1].weight)
        nn.init.zeros_(self.head[-1].bias)           # starts exactly at persistence

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x: [B, H, F] scaled features, channel 0 the scaled CGM. Returns [B] mg/dL."""
        if x.ndim != 3 or x.shape[-1] != self.n_features:
            raise ValueError(f"expected [B, H, {self.n_features}], got {tuple(x.shape)}")
        out, _ = self.gru(x)
        delta = self.head(out[:, -1]).squeeze(-1)
        last_cgm = x[:, -1, 0] * CGM_SCALE[1] + CGM_SCALE[0]
        return last_cgm + delta * CGM_SCALE[1]


def features(history, basal=None, bolus=None, carbs=None) -> torch.Tensor:
    """Stack a window set's arrays into [B, H, F] with the fixed scaling. Covariates are
    optional so the CGM-only variant uses the same path."""
    cols = [(torch.as_tensor(history, dtype=torch.float32) - CGM_SCALE[0]) / CGM_SCALE[1]]
    for name, arr in (("basal", basal), ("bolus", bolus), ("carbs", carbs)):
        if arr is not None:
            cols.append(torch.as_tensor(arr, dtype=torch.float32) / COV_SCALE[name])
    return torch.stack(cols, dim=-1)
