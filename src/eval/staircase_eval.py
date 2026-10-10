"""Part D: evaluate every step of the staircase under both variants, from the aggregate rows
in results/sequence_runs.jsonl and the cached predictions.

Writes
  results/sequence_steps.csv        the staircase: step x variant x cohort x horizon, RMSE/MAE/MAPE
                                    as mean of per-patient values over seeds (mean +- SD), n
  results/sequence_per_patient.csv  per patient x step x variant x horizon, seeds mean +- SD,
                                    n scored and fallback count
  results/sequence_paired.csv       paired gains of a step over persistence and the linear AR,
                                    per patient and bootstrap over patients (n stated)
  results/sequence_ranking.csv      the final model's conforming four-score sum placed in the
                                    official BGLP 2020 ranking (docs/published_comparison.md)
Aggregate numbers only.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict

import numpy as np

from src.data.loader import COHORT_2020, COHORT_OF, PATIENTS, REPO_ROOT
from src.data.windows import HORIZONS, windows_for
from src.eval.baselines import design, fit_ols, persistence
from src.eval.metrics import mae, mape, rmse
from src.staircase import ENSEMBLE_H, PRED_DIR, RUNS, SEEDS5, load_selection

OFFICIAL = [("Rubin-Falcone", 86.31), ("Hameed", 87.15), ("Zhu", 88.12), ("Yang", 88.41), ("Bevan", 89.45),
            ("Sun", 90.36), ("Joedicke", 93.96), ("Ma", 95.85)]          # docs/published_comparison.md, Table 1


def rows_of(step: str) -> list[dict]:
    return [r for r in (json.loads(l) for l in RUNS.read_text().splitlines() if l.strip()) if r["step"] == step]


def preds_of(key: str, variant: str, pid: str):
    z = np.load(PRED_DIR / f"{key}.npz")
    return z[f"{variant}_{pid}_ts"], z[f"{variant}_{pid}_pred"]


def baseline_preds(pid: str, h: int, variant: str, H: int = 24):
    """persistence and the linear AR (per-patient, H=24, +cov) on the same targets."""
    ws = windows_for(pid, "test", h, H, variant=variant)
    tr = windows_for(pid, "train", h, H, variant=variant)
    return ws.target_ts, ws.target, persistence(ws), design(ws, True) @ fit_ols([tr], True)


# ------------------------------------------------------------------ per patient per seed

def per_patient_seed(step: str, variant: str, ensemble: bool = False) -> dict:
    """{(pid, h, seed): (ts, y, pred)} for a step. For the H-ensemble, the median over the H
    models' predictions at every target both variants share (bglp: all targets; primary: the
    intersection of each H's kept targets)."""
    out = {}
    rows = rows_of(step)
    if not ensemble:
        for r in rows:
            ts, p = preds_of(r["key"], variant, r["patient"])
            ws = windows_for(r["patient"], "test", r["horizon_steps"], r["H"], variant=variant)
            out[(r["patient"], r["horizon_steps"], r["seed"])] = (ts, ws.target, p)
        return out
    by = defaultdict(list)
    for r in rows:
        by[(r["patient"], r["horizon_steps"], r["seed"])].append(r)
    for (pid, h, seed), rs in by.items():
        series = [preds_of(r["key"], variant, pid) for r in rs]
        common = series[0][0]
        for ts, _ in series[1:]:
            common = np.intersect1d(common, ts)
        stack = np.stack([p[np.isin(ts, common)] for ts, p in series])
        ws = windows_for(pid, "test", h, min(r["H"] for r in rs), variant=variant)
        y = ws.target[np.isin(ws.target_ts, common)]
        out[(pid, h, seed)] = (common, y, np.median(stack, axis=0))
    return out


def summarize(step_label: str, pp: dict, variant: str) -> tuple[list[dict], list[dict]]:
    """Per-patient rows (seed mean +- SD) and cohort rows (mean of patients, then over seeds)."""
    per, coh = [], []
    for h in HORIZONS:
        pats = sorted({k[0] for k in pp if k[1] == h})
        seeds = sorted({k[2] for k in pp if k[1] == h})
        metric = {}
        for pid in pats:
            for s in seeds:
                if (pid, h, s) in pp:
                    ts, y, p = pp[(pid, h, s)]
                    ws = windows_for(pid, "test", h, 24, variant=variant)
                    fb = int(ws.fallback[np.isin(ws.target_ts, ts)].sum())
                    metric[(pid, s)] = dict(rmse=rmse(y, p), mae=mae(y, p), mape=mape(y, p), n=len(y), fb=fb)
            vals = [metric[(pid, s)] for s in seeds if (pid, s) in metric]
            per.append(dict(step=step_label, variant=variant, horizon_min=h * 5, patient=pid, cohort=COHORT_OF[pid],
                            n_seeds=len(vals), n_scored=vals[0]["n"], n_fallback=vals[0]["fb"],
                            rmse_mean=round(np.mean([v["rmse"] for v in vals]), 3),
                            rmse_sd=round(np.std([v["rmse"] for v in vals], ddof=1), 3) if len(vals) > 1 else 0.0,
                            mae_mean=round(np.mean([v["mae"] for v in vals]), 3),
                            mape_mean=round(np.mean([v["mape"] for v in vals]), 3)))
        for cohort, pids in (("2020", [p for p in pats if p in COHORT_2020]), ("all12", pats)):
            if len(pids) < (6 if cohort == "2020" else 12):
                continue
            seed_means = {m: [np.mean([metric[(pid, s)][m] for pid in pids]) for s in seeds
                              if all((pid, s) in metric for pid in pids)] for m in ("rmse", "mae", "mape")}
            if not seed_means["rmse"]:
                continue
            coh.append(dict(step=step_label, variant=variant, horizon_min=h * 5, cohort=cohort, n_patients=len(pids),
                            n_seeds=len(seed_means["rmse"]), n_scored=sum(metric[(pid, seeds[0])]["n"] for pid in pids),
                            n_fallback=sum(metric[(pid, seeds[0])]["fb"] for pid in pids),
                            rmse_mean=round(np.mean(seed_means["rmse"]), 3),
                            rmse_sd=round(np.std(seed_means["rmse"], ddof=1), 3) if len(seed_means["rmse"]) > 1 else 0.0,
                            mae_mean=round(np.mean(seed_means["mae"]), 3), mape_mean=round(np.mean(seed_means["mape"]), 3)))
    return per, coh


def paired(step_label: str, pp: dict, variant: str, n_boot: int = 10000, seed: int = 0) -> list[dict]:
    """Per-patient gain of the step (seed-mean prediction) over persistence and the linear AR
    on identical targets, and a bootstrap over patients of the mean gain."""
    rng = np.random.default_rng(seed)
    out = []
    for h in HORIZONS:
        pats = sorted({k[0] for k in pp if k[1] == h})
        gains = {"persistence": [], "linear AR": []}
        per_pat = []
        for pid in pats:
            seeds = [k[2] for k in pp if k[0] == pid and k[1] == h]
            ts0, y0, _ = pp[(pid, h, seeds[0])]
            p_model = np.mean([pp[(pid, h, s)][2] for s in seeds], axis=0)
            bts, by, bp0, bl2 = baseline_preds(pid, h, variant)
            m = np.isin(bts, ts0)
            assert np.array_equal(bts[m], ts0) and np.allclose(by[m], y0)
            r_model = rmse(y0, p_model)
            g0, g2 = rmse(y0, bp0[m]) - r_model, rmse(y0, bl2[m]) - r_model
            gains["persistence"].append(g0); gains["linear AR"].append(g2)
            per_pat.append(dict(step=step_label, variant=variant, horizon_min=h * 5, patient=pid, cohort=COHORT_OF[pid],
                                n_scored=len(y0), rmse_model=round(r_model, 3), rmse_persistence=round(rmse(y0, bp0[m]), 3),
                                rmse_linear_ar=round(rmse(y0, bl2[m]), 3), gain_vs_persistence=round(g0, 3),
                                gain_vs_linear_ar=round(g2, 3)))
        out += per_pat
        for cohort, idx in (("2020", [i for i, p in enumerate(pats) if p in COHORT_2020]), ("all12", list(range(len(pats))))):
            for ref, g in gains.items():
                g = np.array(g)[idx]
                boots = [np.mean(g[rng.integers(0, len(g), len(g))]) for _ in range(n_boot)]
                out.append(dict(step=step_label, variant=variant, horizon_min=h * 5, patient=f"bootstrap over {len(g)} patients",
                                cohort=cohort, n_scored=len(g), rmse_model="", rmse_persistence="", rmse_linear_ar="",
                                gain_vs_persistence=round(float(np.mean(g)), 3) if ref == "persistence" else "",
                                gain_vs_linear_ar=round(float(np.mean(g)), 3) if ref == "linear AR" else "",
                                reference=ref, ci_low=round(float(np.percentile(boots, 2.5)), 3),
                                ci_high=round(float(np.percentile(boots, 97.5)), 3),
                                patients_improved=int((g > 0).sum())))
    return out


def four_score(coh_rows: list[dict], step_label: str) -> dict | None:
    r = {c["horizon_min"]: c for c in coh_rows if c["step"] == step_label and c["variant"] == "bglp" and c["cohort"] == "2020"}
    if set(r) != {30, 60}:
        return None
    total = r[30]["rmse_mean"] + r[30]["mae_mean"] + r[60]["rmse_mean"] + r[60]["mae_mean"]
    rank = 1 + sum(1 for _, s in OFFICIAL if s < total)
    return dict(step=step_label, rmse30=r[30]["rmse_mean"], mae30=r[30]["mae_mean"], rmse60=r[60]["rmse_mean"],
                mae60=r[60]["mae_mean"], four_score_sum=round(total, 3), position_in_official_ranking=rank,
                of=len(OFFICIAL) + 1)


STEP_LABELS = [("step1", "step1", False, "1 single configuration (GRU 64, H=24, 3 seeds)"),
               ("step2", "step2", False, "2 validation-selected configuration (1 seed)"),
               ("step3", "step3", False, "3 full-horizon loss (1 seed; kept where it won on validation)"),
               ("step4", "step4", False, "4 chosen configuration, 5 seeds"),
               ("step4ens", "step4ens", True, "4e history-length ensemble (median over H)"),
               ("step5ft", "step5ft", False, "5 pre-trained on all twelve, fine-tuned per patient (5 seeds)"),
               ("step5pop", "step5pop", False, "5p population model alone, LOPO (1 seed)"),
               ("step6seeds", "step6seeds", False, "6 TCN, chosen configuration, 5 seeds"),
               ("step6ens", "step6ens", True, "6e TCN history-length ensemble"),
               ("step7", "step7", False, "7 final model, 5 seeds"),
               ("step7ens", "step7ens", True, "7e final model, H-ensemble"),
               ("step7nocov", "step7nocov", False, "7a final configuration without covariates (5 seeds)")]


def main() -> None:
    per_all, coh_all, paired_all, ranks = [], [], [], []
    for step, label_key, ens, label in STEP_LABELS:
        if not rows_of(step):
            continue
        for variant in ("primary", "bglp"):
            pp = per_patient_seed(step, variant, ens)
            per, coh = summarize(label, pp, variant)
            per_all += per; coh_all += coh
            paired_all += paired(label, pp, variant)
        fs = four_score(coh_all, label)
        if fs:
            ranks.append(fs)
        print(f"{label}: " + "; ".join(f"{c['variant']} {c['cohort']} h={c['horizon_min']} {c['rmse_mean']:.2f}±{c['rmse_sd']:.2f}"
                                       for c in coh_all if c["step"] == label and c["cohort"] == "2020"))
    res = REPO_ROOT / "results"
    for data, fn in ((coh_all, "sequence_steps.csv"), (per_all, "sequence_per_patient.csv"),
                     (paired_all, "sequence_paired.csv"), (ranks, "sequence_ranking.csv")):
        if data:
            keys = list({k: None for r in data for k in r})
            with open(res / fn, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(data)
            print(f"wrote results/{fn}")
    # stop conditions (D-032): 60-min error below 30-min; conforming easier than primary
    bad = []
    for r in per_all:
        if r["horizon_min"] == 30:
            r60 = next((q for q in per_all if (q["step"], q["variant"], q["patient"]) == (r["step"], r["variant"], r["patient"])
                        and q["horizon_min"] == 60), None)
            if r60 and r60["rmse_mean"] <= r["rmse_mean"]:
                bad.append(("60 <= 30", r["step"], r["variant"], r["patient"]))
    for c in coh_all:
        if c["variant"] == "bglp":
            p = next((q for q in coh_all if (q["step"], q["cohort"], q["horizon_min"]) == (c["step"], c["cohort"], c["horizon_min"])
                      and q["variant"] == "primary"), None)
            if p and c["rmse_mean"] < p["rmse_mean"] - 0.5:
                bad.append(("conforming easier than primary", c["step"], c["cohort"], c["horizon_min"]))
    print("stop-condition checks:", "none fired" if not bad else bad)
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
