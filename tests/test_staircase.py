"""The staircase's building blocks (D-032): every family starts at persistence, the
multi-step head and masked loss, full-horizon targets, job keys. Synthetic only."""
from __future__ import annotations

import numpy as np
import pytest
import torch

from src.models.sequence import CGM_SCALE, GlucoseGRU, build_model, features
from src.train import _loss, set_seed, temporal_split, train_one
from src.staircase import GRID, STEP1, job_key, make_job


@pytest.mark.parametrize("family,cfg", [("gru", {}), ("lstm", dict(layers=2)), ("tcn", dict(levels=3, kernel=3))])
@pytest.mark.parametrize("n_out", [1, 6])
def test_every_family_starts_at_persistence_for_every_output_step(family, cfg, n_out):
    set_seed(0)
    m = build_model(family, 4, n_out, hidden=32, dropout=0.2, **cfg)
    m.eval()
    x = torch.randn(5, 24, 4)
    y = m(x)
    assert y.shape == (5, n_out)
    last = x[:, -1, 0] * CGM_SCALE[1] + CGM_SCALE[0]
    assert torch.allclose(y, last[:, None].expand(5, n_out))
    with pytest.raises(ValueError):
        m(torch.randn(5, 24, 3))


def test_step1_class_keeps_the_friday_contract():
    m = GlucoseGRU(n_features=4)
    assert m(torch.randn(3, 24, 4)).shape == (3,)


def test_masked_loss_ignores_steps_without_a_reading():
    pred = torch.tensor([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
    y = torch.tensor([[1.0, float("nan"), 5.0], [float("nan"), float("nan"), 3.0]])
    assert float(_loss(pred, y)) == pytest.approx((0 + 4 + 0) / 3)
    assert float(_loss(pred, torch.tensor([3.0, 3.0]))) == 0.0        # 1-D target: final step


def test_full_horizon_training_reaches_the_final_step():
    rng = np.random.default_rng(0)
    hist = rng.uniform(80, 250, (300, 6))
    y = np.stack([hist[:, -1] + k * 2.0 for k in range(1, 7)], axis=1)
    y[::7, 2] = np.nan                                              # a missing intermediate step
    x = features(hist)
    m, info = train_one(x[:240], torch.as_tensor(y[:240], dtype=torch.float32),
                        x[240:], torch.as_tensor(y[240:], dtype=torch.float32), seed=0, hidden=16, max_epochs=3, patience=3)
    assert m.n_out == 6 and info["epochs_run"] <= 3
    assert m(x[:4]).shape == (4, 6)


def test_job_keys_are_deterministic_and_distinct():
    a = make_job("step1", "per-patient", "559", 6, STEP1, 0)
    b = make_job("step1", "per-patient", "559", 6, STEP1, 0)
    c = make_job("step1", "per-patient", "559", 6, dict(STEP1, hidden=128), 0)
    assert a["key"] == b["key"] == job_key(a) and a["key"] != c["key"]
    assert "H24" in a["key"] and "s0" in a["key"]
    assert [ax for ax, _ in GRID] == ["H", "hidden", "family", "layers", "dropout", "lr"]
