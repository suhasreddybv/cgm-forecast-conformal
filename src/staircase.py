"""The staircase: the sequence model built in steps, every step evaluated under both
variants, nothing selected on test (D-032 onward).

    python -m src.staircase step1                 # single configuration, 3 seeds
    python -m src.staircase step2                 # coordinate search on validation
    python -m src.staircase step3 | step4 | step5 | step6 | step7

Every fit is a job; a job's row (aggregate metrics only, no timestamps, no values) is
appended to results/sequence_runs.jsonl and its predictions go to the out-of-tree cache
(CACHE_DIR/predictions) for ensembling and clinical scoring. Re-running a step skips jobs
whose rows exist. Selections are written to results/sequence_selection.json.

Validation only: a configuration is chosen by its hold-out RMSE (final step); test figures
are computed for every job and never consulted for a choice.
"""
from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

from src.data.loader import CACHE_DIR, COHORT_2018, COHORT_2020, COHORT_OF, PATIENTS, REPO_ROOT
from src.data.windows import HORIZONS, windows_for
from src.eval.metrics import mae, mape, rmse
from src.train import fit_sets, predict, train_one

RUNS = REPO_ROOT / "results" / "sequence_runs.jsonl"
SELECTION = REPO_ROOT / "results" / "sequence_selection.json"
PRED_DIR = CACHE_DIR / "predictions"
WORKERS = int(os.environ.get("STAIRCASE_WORKERS", "3"))
THREADS = int(os.environ.get("STAIRCASE_THREADS", "3"))

STEP1 = dict(family="gru", hidden=64, layers=1, dropout=0.2, H=24, lr=1e-3)
PATIENCE, MAX_EPOCHS, BATCH = 20, 300, 256
SEEDS3 = (0, 1, 2)
SEEDS5 = (20261010, 20261011, 20261012, 20261013, 20261014)
ENSEMBLE_H = (12, 18, 24, 30, 36, 42)
# the brief's grid, searched one axis at a time in this order (D-033)
GRID = [("H", (12, 18, 24, 36, 42)), ("hidden", (32, 64, 128, 256)), ("family", ("gru", "lstm")),
        ("layers", (1, 2)), ("dropout", (0.0, 0.1, 0.2, 0.3)), ("lr", (3e-4, 1e-3))]
TCN_GRID = [("H", (12, 18, 24, 36, 42)), ("hidden", (32, 64, 128)), ("levels", (2, 3, 4)),
            ("kernel", (2, 3, 5)), ("dropout", (0.0, 0.1, 0.2, 0.3)), ("lr", (3e-4, 1e-3))]


# ------------------------------------------------------------------ jobs

def job_key(j: dict) -> str:
    cfg = "-".join(f"{k}{j['cfg'][k]}" for k in sorted(j["cfg"]))
    return (f"{j['step']}_{j['fit']}_{j['patient']}_h{j['horizon']}_{cfg}_fh{int(j['full_horizon'])}"
            f"_cov{int(j['covariates'])}_s{j['seed']}" + (f"_from-{j['init_from']}" if j.get("init_from") else ""))


def make_job(step, fit, patient, horizon, cfg, seed, full_horizon=False, covariates=True, init_from=None,
             patients=None, extra_2018_test=False) -> dict:
    j = dict(step=step, fit=fit, patient=patient, horizon=horizon, cfg=dict(cfg), seed=seed,
             full_horizon=full_horizon, covariates=covariates, init_from=init_from,
             patients=patients, extra_2018_test=extra_2018_test)
    j["key"] = job_key(j)
    return j


def existing_rows() -> dict[str, dict]:
    if not RUNS.exists():
        return {}
    return {r["key"]: r for r in (json.loads(line) for line in RUNS.read_text().splitlines() if line.strip())}


def run_job(j: dict) -> dict:
    torch.set_num_threads(THREADS)
    t0 = time.time()
    cfg = dict(j["cfg"])
    H, family = cfg.pop("H"), cfg.pop("family")
    lr = cfg.pop("lr")
    population = j["fit"] in ("population-LOPO", "pretrain-all12")
    patients = j.get("patients")
    if j["fit"] == "pretrain-all12":
        patients = list(PATIENTS)
    xf, yf, xv, yv = fit_sets(None if j["fit"] == "pretrain-all12" else j["patient"], population, j["horizon"], H,
                              j["covariates"], "primary", j["full_horizon"], patients, j.get("extra_2018_test", False))
    init_state = None
    if j.get("init_from"):
        init_state = torch.load(PRED_DIR / f"{j['init_from']}.pt", weights_only=True)
    model, info = train_one(xf, yf, xv, yv, j["seed"], lr=lr, batch=BATCH, max_epochs=MAX_EPOCHS,
                            patience=PATIENCE, family=family, init_state=init_state, **cfg)
    row = dict(key=j["key"], step=j["step"], fit=j["fit"], patient=j["patient"], cohort=COHORT_OF.get(j["patient"], ""),
               horizon_steps=j["horizon"], family=family, H=H, lr=lr, **cfg, full_horizon=j["full_horizon"],
               covariates=j["covariates"], init_from=j.get("init_from"), seed=j["seed"], n_fit=len(yf), n_val=len(yv),
               epochs_run=info["epochs_run"], best_epoch=info["best_epoch"], val_rmse=info["val_rmse"])
    preds = {}
    for variant in ("primary", "bglp"):
        pids = [j["patient"]] if j["patient"] else list(PATIENTS)
        for pid in pids:
            ws = windows_for(pid, "test", j["horizon"], H, variant=variant)
            y, p = predict(model, ws, j["covariates"])
            preds[f"{variant}_{pid}_ts"] = ws.target_ts
            preds[f"{variant}_{pid}_pred"] = p
            if j["patient"]:
                row.update({f"{variant}_n": len(y), f"{variant}_rmse": round(rmse(y, p), 3),
                            f"{variant}_mae": round(mae(y, p), 3), f"{variant}_mape": round(mape(y, p), 3),
                            f"{variant}_fallback": int(ws.fallback.sum())})
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(PRED_DIR / f"{j['key']}.npz", **preds)
    if j["fit"] == "pretrain-all12":
        torch.save(model.state_dict(), PRED_DIR / f"{j['key']}.pt")
    row["seconds"] = round(time.time() - t0, 1)
    row["val_curve_tail"] = info["val_curve"][-5:]
    return row


def run_jobs(jobs: list[dict], label: str) -> dict[str, dict]:
    rows = existing_rows()
    todo = [j for j in jobs if j["key"] not in rows]
    print(f"[{label}] {len(jobs)} jobs, {len(todo)} to run, {WORKERS} workers", flush=True)
    if todo:
        t0 = time.time()
        with mp.get_context("spawn").Pool(WORKERS) as pool, open(RUNS, "a") as f:
            for i, row in enumerate(pool.imap_unordered(run_job, todo), 1):
                f.write(json.dumps(row) + "\n"); f.flush()
                rows[row["key"]] = row
                flag = " **conforming easier than primary**" if row.get("bglp_rmse", 1e9) < row.get("primary_rmse", 0) - 0.5 else ""
                print(f"  {i}/{len(todo)} {row['patient'] or 'all'} h={row['horizon_steps']} s={row['seed']} "
                      f"val {row['val_rmse']:.2f} primary {row.get('primary_rmse', float('nan')):.2f} "
                      f"bglp {row.get('bglp_rmse', float('nan')):.2f} ({row['epochs_run']} ep, {row['seconds']:.0f}s)"
                      f"{flag}  [{(time.time() - t0) / 60:.1f} min]", flush=True)
    return {j["key"]: rows[j["key"]] for j in jobs}


def load_selection() -> dict:
    return json.loads(SELECTION.read_text()) if SELECTION.exists() else {}


def save_selection(sel: dict) -> None:
    SELECTION.write_text(json.dumps(sel, indent=1, sort_keys=True))


def summarize(rows: list[dict], label: str) -> None:
    for variant in ("bglp", "primary"):
        for h in HORIZONS:
            for name, pids in (("2020", COHORT_2020), ("all12", PATIENTS)):
                per_seed = {}
                for r in rows:
                    if r["horizon_steps"] == h and r["patient"] in pids:
                        per_seed.setdefault(r["seed"], []).append(r[f"{variant}_rmse"])
                means = [np.mean(v) for v in per_seed.values() if len(v) == len(pids)]
                if means:
                    print(f"  {label} {variant:<8} h={h:>2} {name:<6} RMSE mean of patients {np.mean(means):.2f} "
                          f"± {np.std(means, ddof=1) if len(means) > 1 else 0:.2f} over {len(means)} seeds")


# ------------------------------------------------------------------ steps

def step1(args) -> None:
    jobs = [make_job("step1", "per-patient", pid, h, STEP1, s) for pid in PATIENTS for h in HORIZONS for s in SEEDS3]
    rows = run_jobs(jobs, "step1")
    summarize(list(rows.values()), "step1")
    v = [np.mean([rows[make_job("step1", "per-patient", pid, 6, STEP1, s)["key"]]["bglp_rmse"] for pid in COHORT_2020])
         for s in SEEDS3]
    m = float(np.mean(v))
    print(f"\nPRE-REGISTERED CHECK (single configuration, conforming, 2020, 30 min): {m:.2f} mg/dL "
          f"({'inside' if 18.5 <= m <= 19.5 else 'OUTSIDE'} the 18.5-19.5 range registered in D-030)")


def coordinate_search(step: str, base: dict, grid, full_horizon: bool, seed: int = 0) -> dict:
    """One axis at a time, in grid order, selecting per patient-horizon on validation."""
    sel = load_selection()
    chosen = {f"{pid}_h{h}": dict(base) for pid in PATIENTS for h in HORIZONS}
    for axis, values in grid:
        jobs = []
        for pid in PATIENTS:
            for h in HORIZONS:
                for v in values:
                    cfg = dict(chosen[f"{pid}_h{h}"]); cfg[axis] = v
                    jobs.append(make_job(step, "per-patient", pid, h, cfg, seed, full_horizon))
        rows = run_jobs(jobs, f"{step} axis={axis}")
        for pid in PATIENTS:
            for h in HORIZONS:
                cands = [(rows[j["key"]]["val_rmse"], j["cfg"]) for j in jobs if j["patient"] == pid and j["horizon"] == h]
                cands.sort(key=lambda c: c[0])
                chosen[f"{pid}_h{h}"] = dict(cands[0][1])
                print(f"    {pid} h={h} {axis}: " + ", ".join(f"{c[1][axis]}={c[0]:.2f}" for c in cands) + f" -> {cands[0][1][axis]}")
    sel[step] = {k: dict(cfg=v, val_rmse=None) for k, v in chosen.items()}
    # record the chosen config's validation RMSE
    for k, v in sel[step].items():
        pid, h = k.split("_h")
        v["val_rmse"] = existing_rows()[make_job(step, "per-patient", pid, int(h), v["cfg"], seed, full_horizon)["key"]]["val_rmse"]
    save_selection(sel)
    return chosen


def step2(args) -> None:
    chosen = coordinate_search("step2", STEP1, GRID, False)
    rows = existing_rows()
    print("\nchosen per patient-horizon (validation), with its test figures for the record:")
    for k, cfg in chosen.items():
        pid, h = k.split("_h")
        r = rows[make_job("step2", "per-patient", pid, int(h), cfg, 0)["key"]]
        print(f"  {k}: {cfg}  val {r['val_rmse']:.2f} primary {r['primary_rmse']:.2f} bglp {r['bglp_rmse']:.2f}")


def _chosen(step: str) -> dict:
    sel = load_selection()
    if step not in sel:
        sys.exit(f"{step} has not been run; nothing to build on")
    return {k: v["cfg"] for k, v in sel[step].items()}


def step3(args) -> None:
    """Full-horizon loss on the Step 2 configuration; keep whichever wins on validation."""
    base = _chosen("step2")
    jobs = [make_job("step3", "per-patient", pid, h, base[f"{pid}_h{h}"], 0, True) for pid in PATIENTS for h in HORIZONS]
    rows = run_jobs(jobs, "step3 full-horizon")
    sel = load_selection(); prev = sel["step2"]; sel["step3"] = {}
    wins = 0
    for j in jobs:
        k = f"{j['patient']}_h{j['horizon']}"
        v_full, v_single = rows[j["key"]]["val_rmse"], prev[k]["val_rmse"]
        use_full = v_full < v_single
        wins += use_full
        sel["step3"][k] = dict(cfg=base[k], full_horizon=bool(use_full), val_rmse=min(v_full, v_single),
                               val_full=v_full, val_single=v_single)
    save_selection(sel)
    print(f"\nfull-horizon loss wins on validation for {wins}/{len(jobs)} patient-horizons")


def step4(args) -> None:
    """5 seeds of the chosen configuration, and the history-length ensemble (1 seed per H)."""
    sel = load_selection()["step3"]
    jobs = []
    for k, v in sel.items():
        pid, h = k.split("_h")
        for s in SEEDS5:
            jobs.append(make_job("step4", "per-patient", pid, int(h), v["cfg"], s, v["full_horizon"]))
        for H in ENSEMBLE_H:
            cfg = dict(v["cfg"]); cfg["H"] = H
            jobs.append(make_job("step4ens", "per-patient", pid, int(h), cfg, SEEDS5[0], v["full_horizon"]))
    run_jobs(jobs, "step4 seeds + H-ensemble")


def step5(args) -> None:
    """Pre-train on all twelve training files (+ 2018 test files), fine-tune per patient;
    and the LOPO population model alone."""
    sel = load_selection()["step3"]
    pre = {}
    for h in HORIZONS:
        # one configuration for the population: the most-chosen per-patient configuration at this horizon
        cfgs = [json.dumps(v["cfg"], sort_keys=True) for k, v in sel.items() if k.endswith(f"_h{h}")]
        cfg = json.loads(max(set(cfgs), key=cfgs.count))
        fh = sum(v["full_horizon"] for k, v in sel.items() if k.endswith(f"_h{h}")) > 6
        for s in SEEDS5:
            pre[(h, s)] = make_job("step5pre", "pretrain-all12", "", h, cfg, s, fh, extra_2018_test=True)
    rows = run_jobs(list(pre.values()), "step5 pre-train all twelve (+2018 test)")
    jobs = []
    for h in HORIZONS:
        for s in SEEDS5:
            p = pre[(h, s)]
            for pid in PATIENTS:
                jobs.append(make_job("step5ft", "fine-tuned", pid, h, p["cfg"], s, p["full_horizon"], init_from=p["key"]))
    run_jobs(jobs, "step5 fine-tune per patient")
    jobs = []
    for h in HORIZONS:
        p = pre[(h, SEEDS5[0])]
        for pid in PATIENTS:
            jobs.append(make_job("step5pop", "population-LOPO", pid, h, p["cfg"], SEEDS5[0], p["full_horizon"]))
    run_jobs(jobs, "step5 population LOPO (1 seed)")


def step6(args) -> None:
    """The second family: a dilated causal CNN, searched the same way, then steps 3-5 for it."""
    base = dict(family="tcn", hidden=64, levels=3, kernel=3, dropout=0.2, H=24, lr=1e-3)
    chosen = coordinate_search("step6", base, TCN_GRID, False)
    jobs = [make_job("step6fh", "per-patient", pid, h, chosen[f"{pid}_h{h}"], 0, True) for pid in PATIENTS for h in HORIZONS]
    rows = run_jobs(jobs, "step6 full-horizon")
    sel = load_selection(); sel["step6fh"] = {}
    for j in jobs:
        k = f"{j['patient']}_h{j['horizon']}"
        v_full, v_single = rows[j["key"]]["val_rmse"], sel["step6"][k]["val_rmse"]
        sel["step6fh"][k] = dict(cfg=chosen[k], full_horizon=bool(v_full < v_single), val_rmse=min(v_full, v_single))
    save_selection(sel)
    jobs = []
    for k, v in sel["step6fh"].items():
        pid, h = k.split("_h")
        for s in SEEDS5:
            jobs.append(make_job("step6seeds", "per-patient", pid, int(h), v["cfg"], s, v["full_horizon"]))
        for H in ENSEMBLE_H:
            cfg = dict(v["cfg"]); cfg["H"] = H
            jobs.append(make_job("step6ens", "per-patient", pid, int(h), cfg, SEEDS5[0], v["full_horizon"]))
    run_jobs(jobs, "step6 seeds + H-ensemble")


STEPS = dict(step1=step1, step2=step2, step3=step3, step4=step4, step5=step5, step6=step6)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=list(STEPS))
    a = ap.parse_args()
    RUNS.parent.mkdir(exist_ok=True)
    STEPS[a.step](a)


if __name__ == "__main__":
    main()
