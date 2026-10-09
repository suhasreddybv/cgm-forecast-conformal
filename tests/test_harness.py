"""The model harness (D-025, D-026): built and exercised on synthetic tensors only.
No fitting on real data happens here or anywhere before Saturday."""
from __future__ import annotations

import numpy as np
import pytest
import torch

from src.data.windows import WindowSet
from src.models.sequence import CGM_SCALE, GlucoseGRU, features
from src.train import SEEDS, VAL_FRACTION, set_seed, temporal_split, tensors, train_one


def test_forward_pass_shapes_and_persistence_initialisation():
    set_seed(0)
    m = GlucoseGRU(n_features=4)
    x = torch.randn(8, 24, 4)
    y = m(x)
    assert y.shape == (8,)
    # the head is zero-initialised, so the untrained model IS persistence
    last = x[:, -1, 0] * CGM_SCALE[1] + CGM_SCALE[0]
    assert torch.allclose(y, last)
    with pytest.raises(ValueError):
        m(torch.randn(8, 24, 3))


def test_features_stack_with_fixed_scaling():
    h = np.full((5, 6), 140.0)
    x = features(h, np.ones((5, 6)), np.zeros((5, 6)), np.full((5, 6), 20.0))
    assert x.shape == (5, 6, 4)
    assert torch.allclose(x[..., 0], torch.zeros(5, 6)) and float(x[0, 0, 3]) == 1.0
    assert features(h).shape == (5, 6, 1)


def _ws(n=100, H=6, nan_rows=()):
    rng = np.random.default_rng(0)
    hist = rng.uniform(60, 300, (n, H))
    basal = np.ones((n, H))
    for r in nan_rows:
        basal[r, 0] = np.nan
    ts = np.datetime64("2021-01-01T00:00:00", "s") + (np.arange(n) * 300).astype("timedelta64[s]")
    return WindowSet("000", "train", 6, H, target_ts=ts, target=hist[:, -1] + rng.normal(0, 5, n),
                     history=hist, max_bracket_s=np.zeros(n, int), basal=basal,
                     bolus=np.zeros((n, H)), carbs=np.zeros((n, H)),
                     n_readings=n, n_candidates=n, n_dropped_long_gap=0, n_dropped_no_history=0)


def test_temporal_split_takes_the_tail_and_never_overlaps():
    ws = _ws(100)
    fit, val = temporal_split(ws)
    assert len(val) == 20 and len(fit) == 80 and VAL_FRACTION == 0.2
    assert val.min() == 80 and ws.target_ts[val].min() > ws.target_ts[fit].max()
    assert not set(fit) & set(val)


def test_nan_covariate_windows_are_dropped_only_with_covariates():
    ws = _ws(50, nan_rows=(3, 7))
    x, y, keep = tensors(ws, np.arange(50), covariates=True)
    assert len(keep) == 48 and 3 not in keep and 7 not in keep and not torch.isnan(x).any()
    x2, y2, keep2 = tensors(ws, np.arange(50), covariates=False)
    assert len(keep2) == 50 and x2.shape[-1] == 1


def test_training_loop_runs_on_synthetic_tensors_and_is_seed_reproducible():
    rng = np.random.default_rng(1)
    n, H = 400, 6
    hist = rng.uniform(80, 250, (n, H))
    target = hist[:, -1] + 0.5 * (hist[:, -1] - hist[:, -2]) + rng.normal(0, 3, n)
    x = features(hist)
    y = torch.as_tensor(target, dtype=torch.float32)
    a, info_a = train_one(x[:320], y[:320], x[320:], y[320:], seed=0, hidden=16, max_epochs=5, patience=3)
    b, info_b = train_one(x[:320], y[:320], x[320:], y[320:], seed=0, hidden=16, max_epochs=5, patience=3)
    assert info_a["val_curve"] == info_b["val_curve"], "same seed, same curve"
    assert all(torch.equal(p, q) for p, q in zip(a.state_dict().values(), b.state_dict().values()))
    c, info_c = train_one(x[:320], y[:320], x[320:], y[320:], seed=1, hidden=16, max_epochs=5, patience=3)
    assert info_c["val_curve"] != info_a["val_curve"], "different seed, different curve"
    assert info_a["epochs_run"] <= 5 and info_a["seed"] == 0
    assert len(SEEDS) >= 3
