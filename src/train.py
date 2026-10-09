"""Training harness for the sequence model - written and logged on 9 October, not run
until Saturday 10 October (D-025, D-026).

    python -m src.train --patient 559 --horizon 6 --history 24 --seed 0
    python -m src.train --population --held-out 559 --horizon 6 --history 24 --seed 0

Rules it enforces, all from docs/protocol.md:
  - inputs are the Day 2 windows unchanged; primary models see CGM, basal, bolus, carbs;
  - the validation hold-out is the last VAL_FRACTION of each training file's windows by
    time (protocol clarification of 9 Oct); early stopping and checkpoint choice use it
    and nothing else; the test file is loaded only for the final score;
  - population models exclude the held-out patient's whole training file;
  - windows with a NaN covariate are dropped from fitting and from scoring (D-020);
  - the seed is set for torch, numpy and python before anything random happens, and is
    written into the result row.
"""
from __future__ import annotations

import argparse
import json
import random
import time

import numpy as np
import torch

from src.data.loader import PATIENTS, REPO_ROOT
from src.data.windows import WindowSet, windows_for
from src.eval.baselines import covariate_ok
from src.eval.metrics import mae, mape, rmse
from src.models.sequence import GlucoseGRU, features

VAL_FRACTION = 0.2          # protocol clarification, 9 Oct 2026
SEEDS = (0, 1, 2)           # the minimum three; more may be added, none dropped


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def temporal_split(ws: WindowSet, val_fraction: float = VAL_FRACTION):
    """Indices of (fit, hold-out) windows: the hold-out is the last val_fraction of the
    windows in time order. Windows are already time-ordered by construction; asserted."""
    order = np.argsort(ws.target_ts, kind="stable")
    assert (order == np.arange(len(ws))).all(), "windows are not in time order"
    n_val = int(round(len(ws) * val_fraction))
    cut = len(ws) - n_val
    fit, val = np.arange(cut), np.arange(cut, len(ws))
    if n_val and cut:
        assert ws.target_ts[val].min() > ws.target_ts[fit].max()
    return fit, val


def tensors(ws: WindowSet, idx: np.ndarray, covariates: bool):
    keep = idx[covariate_ok(ws)[idx]] if covariates else idx
    x = features(ws.history[keep], *((ws.basal[keep], ws.bolus[keep], ws.carbs[keep]) if covariates else ()))
    y = torch.as_tensor(ws.target[keep], dtype=torch.float32)
    return x, y, keep


def fit_sets(patient: str | None, population: bool, horizon: int, history: int, covariates: bool):
    """(fit tensors, hold-out tensors) under the protocol's split rules."""
    pids = [p for p in PATIENTS if p != patient] if population else [patient]
    xf, yf, xv, yv = [], [], [], []
    for pid in pids:
        ws = windows_for(pid, "train", horizon, history)
        fit, val = temporal_split(ws)
        a, b, _ = tensors(ws, fit, covariates)
        c, d, _ = tensors(ws, val, covariates)
        xf.append(a); yf.append(b); xv.append(c); yv.append(d)
    return torch.cat(xf), torch.cat(yf), torch.cat(xv), torch.cat(yv)


def train_one(x_fit, y_fit, x_val, y_val, seed: int, hidden: int = 64, layers: int = 1,
              lr: float = 1e-3, batch: int = 256, max_epochs: int = 100, patience: int = 10,
              device: str = "cpu") -> tuple[GlucoseGRU, dict]:
    set_seed(seed)
    model = GlucoseGRU(n_features=x_fit.shape[-1], hidden=hidden, layers=layers).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    best, best_state, bad, history = float("inf"), None, 0, []
    g = torch.Generator().manual_seed(seed)
    for epoch in range(max_epochs):
        model.train()
        perm = torch.randperm(len(x_fit), generator=g)
        for i in range(0, len(perm), batch):
            j = perm[i:i + batch]
            opt.zero_grad()
            loss = torch.mean((model(x_fit[j].to(device)) - y_fit[j].to(device)) ** 2)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            val = float(torch.sqrt(torch.mean((model(x_val.to(device)) - y_val.to(device)) ** 2)))
        history.append(val)
        if val < best - 1e-3:
            best, bad = val, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= patience:
                break
    model.load_state_dict(best_state)
    return model, dict(seed=seed, epochs_run=len(history), best_epoch=int(np.argmin(history)) + 1,
                       val_rmse=round(best, 3), val_curve=[round(v, 3) for v in history])


def score_test(model: GlucoseGRU, patient: str, horizon: int, history: int, covariates: bool,
               device: str = "cpu") -> dict:
    """The only place the test file is read, after fitting and selection are over."""
    ws = windows_for(patient, "test", horizon, history)
    x, y, keep = tensors(ws, np.arange(len(ws)), covariates)
    model.eval()
    with torch.no_grad():
        p = model(x.to(device)).cpu().numpy()
    yv = y.numpy()
    return dict(n_scored=len(yv), rmse=round(rmse(yv, p), 3), mae=round(mae(yv, p), 3),
                mape=round(mape(yv, p), 3), n_dropped_nan_covariate=len(ws) - len(keep))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", help="patient-specific model, or the held-out patient with --population")
    ap.add_argument("--population", action="store_true")
    ap.add_argument("--horizon", type=int, required=True, choices=(6, 12))
    ap.add_argument("--history", type=int, required=True, choices=(6, 12, 24))
    ap.add_argument("--no-covariates", action="store_true", help="CGM-only variant (secondary)")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--hidden", type=int, default=64)
    ap.add_argument("--layers", type=int, default=1)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--max-epochs", type=int, default=100)
    ap.add_argument("--patience", type=int, default=10)
    ap.add_argument("--out", default=str(REPO_ROOT / "results" / "sequence_runs.jsonl"))
    a = ap.parse_args(argv)
    if not a.patient:
        ap.error("--patient is required (the held-out patient when --population)")
    cov = not a.no_covariates
    t0 = time.time()
    xf, yf, xv, yv = fit_sets(a.patient, a.population, a.horizon, a.history, cov)
    model, fit_info = train_one(xf, yf, xv, yv, a.seed, a.hidden, a.layers, a.lr,
                                max_epochs=a.max_epochs, patience=a.patience)
    test = score_test(model, a.patient, a.horizon, a.history, cov)
    row = dict(patient=a.patient, fit="population (LOPO)" if a.population else "per-patient",
               horizon_steps=a.horizon, history_steps=a.history, covariates=cov,
               hidden=a.hidden, layers=a.layers, lr=a.lr, val_fraction=VAL_FRACTION,
               n_fit=len(yf), n_val=len(yv), seconds=round(time.time() - t0, 1), **fit_info, **test)
    with open(a.out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps({k: v for k, v in row.items() if k != "val_curve"}))


if __name__ == "__main__":
    main()
