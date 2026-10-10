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
    written into the result row;
  - every run is scored under BOTH evaluation variants (D-027): the primary protocol and
    the conforming `bglp` variant, from the first run. `--fit-variant` says which windows
    the model was fitted on (primary by default); both test scores are written.
"""
from __future__ import annotations

import argparse
import json
import random
import time

import numpy as np
import torch

from src.data.loader import PATIENTS, REPO_ROOT
from src.data.windows import VARIANTS, WindowSet, windows_for
from src.eval.baselines import covariate_ok
from src.eval.metrics import mae, mape, rmse
from src.models.sequence import GlucoseGRU, build_model, features

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


def tensors(ws: WindowSet, idx: np.ndarray, covariates: bool, multi: np.ndarray | None = None):
    """multi: [n, h] intermediate-and-final targets for the full-horizon loss (NaN where no
    real reading exists at that step); when given, y is [n, h], else [n]."""
    keep = idx[covariate_ok(ws)[idx]] if covariates else idx
    x = features(ws.history[keep], *((ws.basal[keep], ws.bolus[keep], ws.carbs[keep]) if covariates else ()))
    y = torch.as_tensor((multi if multi is not None else ws.target)[keep], dtype=torch.float32)
    return x, y, keep


def multi_step_targets(ws: WindowSet, tol_s: int = 60) -> np.ndarray:
    """[n, h] targets at T - (h-k)*300 s for k = 1..h, taken from the patient's own series
    only where a real reading sits within tol_s of the step (NaN otherwise; the loss masks
    them). Column -1 is the scored target. Training signal only - never scored."""
    from src.data.loader import load_patient
    from src.data.windows import STEP_S, _series_for
    cgm, _, _ = _series_for(load_patient(ws.patient_id, ws.split))
    t = cgm.ts.astype("datetime64[s]").astype(np.int64)
    T = ws.target_ts.astype("datetime64[s]").astype(np.int64)
    steps = T[:, None] - (ws.horizon - np.arange(1, ws.horizon + 1))[None, :] * STEP_S
    r = np.minimum(np.searchsorted(t, steps, side="left"), len(t) - 1)
    lft = np.maximum(r - 1, 0)
    near = np.where(np.abs(t[r] - steps) <= np.abs(t[lft] - steps), r, lft)
    ok = np.abs(t[near] - steps) <= tol_s
    out = np.where(ok, cgm.values[near], np.nan)
    assert np.allclose(out[:, -1], ws.target), "the final step must be the scored target"
    return out


def fit_sets(patient: str | None, population: bool, horizon: int, history: int, covariates: bool,
             variant: str = "primary", full_horizon: bool = False, patients=None, extra_2018_test: bool = False):
    """(fit tensors, hold-out tensors) under the protocol's split rules. `patients` overrides
    the patient list (pre-training on all twelve); `extra_2018_test` adds the 2018 cohort's
    test files to the fitting set, which the challenge rules permit for pre-training (D-035).
    The hold-out is always each file's temporal tail; a 2018 test file contributes its tail
    too. Full-horizon targets are [n, h] with NaN where no real reading exists."""
    if patients is None:
        patients = [p for p in PATIENTS if p != patient] if population else [patient]
    files = [(p, "train") for p in patients]
    if extra_2018_test:
        from src.data.loader import COHORT_2018
        files += [(p, "test") for p in COHORT_2018 if p in patients]
    xf, yf, xv, yv = [], [], [], []
    for pid, split in files:
        ws = windows_for(pid, split, horizon, history, variant=variant)
        multi = multi_step_targets(ws) if full_horizon else None
        fit, val = temporal_split(ws)
        a, b, _ = tensors(ws, fit, covariates, multi)
        c, d, _ = tensors(ws, val, covariates, multi)
        xf.append(a); yf.append(b); xv.append(c); yv.append(d)
    return torch.cat(xf), torch.cat(yf), torch.cat(xv), torch.cat(yv)


def _final(y: torch.Tensor) -> torch.Tensor:
    return y[:, -1] if y.ndim == 2 else y


def _loss(pred: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """MSE over the available steps: with the full-horizon target, NaN steps are masked."""
    if y.ndim == 1:
        return torch.mean((_final(pred) - y) ** 2)
    m = torch.isfinite(y)
    return torch.sum(((pred - torch.nan_to_num(y)) ** 2) * m) / m.sum()


def train_one(x_fit, y_fit, x_val, y_val, seed: int, hidden: int = 64, layers: int = 1,
              lr: float = 1e-3, batch: int = 256, max_epochs: int = 100, patience: int = 10,
              device: str = "cpu", family: str = "gru", dropout: float = 0.0, init_state=None,
              **cfg) -> tuple[torch.nn.Module, dict]:
    """Fit on x_fit, stop on x_val's final-step RMSE (the scored quantity). `init_state`
    warm-starts fine-tuning from a pre-trained model (D-035)."""
    set_seed(seed)
    n_out = y_fit.shape[1] if y_fit.ndim == 2 else 1
    model = build_model(family, x_fit.shape[-1], n_out, hidden=hidden, layers=layers, dropout=dropout, **cfg).to(device)
    if init_state is not None:
        model.load_state_dict(init_state)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    best, best_state, bad, history = float("inf"), None, 0, []
    g = torch.Generator().manual_seed(seed)
    yv_final = _final(y_val).to(device)
    for epoch in range(max_epochs):
        model.train()
        perm = torch.randperm(len(x_fit), generator=g)
        for i in range(0, len(perm), batch):
            j = perm[i:i + batch]
            opt.zero_grad()
            loss = _loss(model(x_fit[j].to(device)), y_fit[j].to(device))
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            val = float(torch.sqrt(torch.mean((_final(model(x_val.to(device))) - yv_final) ** 2)))
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


def predict(model, ws: WindowSet, covariates: bool, device: str = "cpu") -> tuple[np.ndarray, np.ndarray]:
    """(target, prediction) on every window of ws that has finite covariates; final step only."""
    x, y, keep = tensors(ws, np.arange(len(ws)), covariates)
    model.eval()
    with torch.no_grad():
        p = _final(model(x.to(device))).cpu().numpy()
    return y.numpy(), p


def score_test(model: GlucoseGRU, patient: str, horizon: int, history: int, covariates: bool,
               device: str = "cpu", variant: str = "primary") -> dict:
    """The only place a test file is read, after fitting and selection are over."""
    ws = windows_for(patient, "test", horizon, history, variant=variant)
    yv, p = predict(model, ws, covariates, device)
    return dict(n_scored=len(yv), rmse=round(rmse(yv, p), 3), mae=round(mae(yv, p), 3),
                mape=round(mape(yv, p), 3), n_dropped_nan_covariate=len(ws) - len(yv),
                n_fallback=int(ws.fallback.sum()))


def score_both_variants(model, patient, horizon, history, covariates, device="cpu") -> dict:
    """Primary keys unprefixed (the frozen protocol); the conforming variant prefixed bglp_."""
    out = dict(score_test(model, patient, horizon, history, covariates, device, "primary"))
    out.update({f"bglp_{k}": v for k, v in
                score_test(model, patient, horizon, history, covariates, device, "bglp").items()})
    return out


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
    ap.add_argument("--fit-variant", default="primary", choices=VARIANTS,
                    help="which windows to fit on; the test is always scored under both variants")
    ap.add_argument("--out", default=str(REPO_ROOT / "results" / "sequence_runs.jsonl"))
    a = ap.parse_args(argv)
    if not a.patient:
        ap.error("--patient is required (the held-out patient when --population)")
    cov = not a.no_covariates
    t0 = time.time()
    xf, yf, xv, yv = fit_sets(a.patient, a.population, a.horizon, a.history, cov, a.fit_variant)
    model, fit_info = train_one(xf, yf, xv, yv, a.seed, a.hidden, a.layers, a.lr,
                                max_epochs=a.max_epochs, patience=a.patience)
    test = score_both_variants(model, a.patient, a.horizon, a.history, cov)
    row = dict(patient=a.patient, fit="population (LOPO)" if a.population else "per-patient",
               fit_variant=a.fit_variant,
               horizon_steps=a.horizon, history_steps=a.history, covariates=cov,
               hidden=a.hidden, layers=a.layers, lr=a.lr, val_fraction=VAL_FRACTION,
               n_fit=len(yf), n_val=len(yv), seconds=round(time.time() - t0, 1), **fit_info, **test)
    with open(a.out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps({k: v for k, v in row.items() if k != "val_curve"}))


if __name__ == "__main__":
    main()
