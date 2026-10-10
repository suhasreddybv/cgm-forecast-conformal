"""Sequence models, defined and smoke-tested before any was fitted (D-026), extended on
10 October for the staircase (D-032): recurrent cell choice, dropout, a multi-step head for
the full-horizon loss, and a dilated causal CNN as the second family.

Every model predicts the CHANGE from the last history slot, scaled, with a zero-initialised
output layer: the untrained model is exactly persistence, the floor every later number is
measured from. Input scaling is a fixed affine map, not a fitted scaler. The multi-step head
returns [B, n_out]; the challenge scores column -1, the horizon's final step.
"""
from __future__ import annotations

import torch
from torch import nn

CGM_SCALE = (140.0, 60.0)        # (centre, scale) mg/dL - fixed constants, not fitted
COV_SCALE = {"basal": 1.0, "bolus": 1.0, "carbs": 20.0}   # units per feature, fixed


def _head(hidden: int, n_out: int, dropout: float) -> nn.Sequential:
    head = nn.Sequential(nn.Dropout(dropout), nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, n_out))
    nn.init.zeros_(head[-1].weight)
    nn.init.zeros_(head[-1].bias)           # starts exactly at persistence
    return head


def _finish(x: torch.Tensor, delta: torch.Tensor) -> torch.Tensor:
    last_cgm = x[:, -1, 0] * CGM_SCALE[1] + CGM_SCALE[0]
    return last_cgm[:, None] + delta * CGM_SCALE[1]


class GlucoseRNN(nn.Module):
    """GRU or LSTM over the H history steps. `n_out` = 1 for the final step only, h for the
    full-horizon loss."""

    def __init__(self, n_features: int = 4, hidden: int = 64, layers: int = 1, dropout: float = 0.0,
                 cell: str = "gru", n_out: int = 1):
        super().__init__()
        self.n_features, self.n_out, self.cell = n_features, n_out, cell
        rnn = {"gru": nn.GRU, "lstm": nn.LSTM}[cell]
        self.rnn = rnn(n_features, hidden, num_layers=layers, batch_first=True,
                       dropout=dropout if layers > 1 else 0.0)
        self.head = _head(hidden, n_out, dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 3 or x.shape[-1] != self.n_features:
            raise ValueError(f"expected [B, H, {self.n_features}], got {tuple(x.shape)}")
        out, _ = self.rnn(x)
        return _finish(x, self.head(out[:, -1]))


class GlucoseGRU(GlucoseRNN):
    """The Step 1 model and the name the Friday smoke test uses."""

    def __init__(self, n_features: int = 4, hidden: int = 64, layers: int = 1, dropout: float = 0.0, n_out: int = 1):
        super().__init__(n_features, hidden, layers, dropout, "gru", n_out)

    def forward(self, x):
        y = super().forward(x)
        return y[:, -1] if self.n_out == 1 else y


class _CausalBlock(nn.Module):
    def __init__(self, ch_in, ch_out, kernel, dilation, dropout):
        super().__init__()
        self.pad = (kernel - 1) * dilation
        self.conv1 = nn.Conv1d(ch_in, ch_out, kernel, dilation=dilation)
        self.conv2 = nn.Conv1d(ch_out, ch_out, kernel, dilation=dilation)
        self.drop = nn.Dropout(dropout)
        self.skip = nn.Conv1d(ch_in, ch_out, 1) if ch_in != ch_out else nn.Identity()

    def forward(self, x):                    # x: [B, C, T]
        y = torch.relu(self.conv1(nn.functional.pad(x, (self.pad, 0))))
        y = self.drop(y)
        y = torch.relu(self.conv2(nn.functional.pad(y, (self.pad, 0))))
        return torch.relu(self.drop(y) + self.skip(x))


class GlucoseTCN(nn.Module):
    """Dilated causal 1-D CNN (Zhu's lineage): `levels` blocks with dilations 1, 2, 4, ...,
    the last time step into the same delta head."""

    def __init__(self, n_features: int = 4, hidden: int = 64, levels: int = 3, kernel: int = 3,
                 dropout: float = 0.0, n_out: int = 1):
        super().__init__()
        self.n_features, self.n_out = n_features, n_out
        blocks, ch = [], n_features
        for i in range(levels):
            blocks.append(_CausalBlock(ch, hidden, kernel, 2 ** i, dropout))
            ch = hidden
        self.blocks = nn.Sequential(*blocks)
        self.head = _head(hidden, n_out, dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 3 or x.shape[-1] != self.n_features:
            raise ValueError(f"expected [B, H, {self.n_features}], got {tuple(x.shape)}")
        z = self.blocks(x.transpose(1, 2))[:, :, -1]
        return _finish(x, self.head(z))


def build_model(family: str, n_features: int, n_out: int = 1, **cfg) -> nn.Module:
    if family in ("gru", "lstm"):
        return GlucoseRNN(n_features, cfg.get("hidden", 64), cfg.get("layers", 1), cfg.get("dropout", 0.0), family, n_out)
    if family == "tcn":
        return GlucoseTCN(n_features, cfg.get("hidden", 64), cfg.get("levels", 3), cfg.get("kernel", 3),
                          cfg.get("dropout", 0.0), n_out)
    raise ValueError(f"unknown family {family!r}")


def features(history, basal=None, bolus=None, carbs=None) -> torch.Tensor:
    """Stack a window set's arrays into [B, H, F] with the fixed scaling. Covariates are
    optional so the CGM-only variant uses the same path."""
    cols = [(torch.as_tensor(history, dtype=torch.float32) - CGM_SCALE[0]) / CGM_SCALE[1]]
    for name, arr in (("basal", basal), ("bolus", bolus), ("carbs", carbs)):
        if arr is not None:
            cols.append(torch.as_tensor(arr, dtype=torch.float32) / COV_SCALE[name])
    return torch.stack(cols, dim=-1)
