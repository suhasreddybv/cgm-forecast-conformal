"""The first fits of the repository: three baselines on the frozen protocol.

  p0  persistence - the value at the last history slot, T - h. Every later method sees
      exactly this reading, so this is the floor they are measured from. The slot value
      is the raw reading wherever a reading sits on the slot, and within 1/301 of a step
      of it otherwise (D-013).
  l1  linear extrapolation - a least-squares line through the last k history slots,
      k in {2, 3, 6}, extended h steps. Shows how much of the task is local trend. All k
      reported, none selected.
  l2  linear autoregression - ordinary least squares on the H history slots, with and
      without the covariates the primary models may use (D-010, D-015), fitted
        per patient   on that patient's training windows;
        population    leave-one-patient-out on the other eleven patients' training
                      windows (protocol, "Fitting"), the evaluated patient held out.
      Coefficients come from training windows only; nothing is scaled, because OLS is
      invariant to feature scaling and a scaler would be one more fitted object.

Evaluated on the test split only, at 30 and 60 minutes, on the Day 2 windows unchanged.
Every figure carries the number of real targets it was computed on (D-009). Windows with
a NaN covariate are dropped from the with-covariate rows and kept in the without rows,
so those two are compared on stated, not identical, n.

Writes results/baselines_per_patient.csv and results/baselines.csv, then runs the three
sanity checks of D-021 and exits non-zero if any fails. Aggregate numbers only.
"""
from __future__ import annotations

import csv
import sys

import numpy as np

from src.data.loader import COHORT_OF, PATIENTS, REPO_ROOT, load_patient
from src.data.windows import HISTORY_LENGTHS, HORIZONS, WindowSet, windows_for
from src.eval.metrics import mae, mape, rmse

K_VALUES = (2, 3, 6)
COVARIATES = ("basal", "bolus", "carbs")


# ------------------------------------------------------------------ predictors

def persistence(ws: WindowSet) -> np.ndarray:
    return ws.history[:, -1]


def linear_extrapolation(ws: WindowSet, k: int) -> np.ndarray:
    """Least-squares line through the last k slots (x = -(k-1) ... 0 steps), at x = h."""
    y = ws.history[:, -k:]
    x = np.arange(-(k - 1), 1, dtype=float)
    xc = x - x.mean()
    slope = (y - y.mean(axis=1, keepdims=True)) @ xc / (xc @ xc)     # per step
    intercept = y.mean(axis=1) - slope * x.mean()                     # value at x = 0
    return intercept + slope * ws.horizon


def design(ws: WindowSet, covariates: bool) -> np.ndarray:
    cols = [ws.history] + ([getattr(ws, c) for c in COVARIATES] if covariates else [])
    X = np.concatenate(cols, axis=1)
    return np.concatenate([np.ones((len(X), 1)), X], axis=1)


def covariate_ok(ws: WindowSet) -> np.ndarray:
    return np.all([np.isfinite(getattr(ws, c)).all(axis=1) for c in COVARIATES], axis=0)


def fit_ols(train: list[WindowSet], covariates: bool) -> np.ndarray:
    X = np.concatenate([design(w, covariates)[covariate_ok(w) if covariates else slice(None)] for w in train])
    y = np.concatenate([w.target[covariate_ok(w) if covariates else slice(None)] for w in train])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


# ------------------------------------------------------------------ evaluation

def _score(row: dict, y: np.ndarray, p: np.ndarray) -> dict:
    return dict(row, n_scored=len(y), rmse=round(rmse(y, p), 3), mae=round(mae(y, p), 3),
                mape=round(mape(y, p), 3))


def per_patient_rows(eval_variant: str = "primary") -> list[dict]:
    rows = []
    test = {(pid, h, H): windows_for(pid, "test", h, H, variant=eval_variant) for pid in PATIENTS for h in HORIZONS for H in HISTORY_LENGTHS}
    train = {(pid, h, H): windows_for(pid, "train", h, H, variant=eval_variant) for pid in PATIENTS for h in HORIZONS for H in HISTORY_LENGTHS}
    if eval_variant == "bglp":
        # every challenge test point is a target: the counts must equal Table 2 exactly
        for (pid, h, H), ws in test.items():
            assert len(ws) == ws.n_candidates == ws.n_readings - load_patient(pid, "test").eval_start_index, (pid, h, H, len(ws))

    # targets are real readings - guaranteed by the windows, asserted anyway
    for (pid, h, H), ws in test.items():
        raw = load_patient(pid, "test").cgm.ts
        assert np.isin(ws.target_ts, raw).all(), f"{pid} h={h} H={H}: a target is not a raw reading"

    for h in HORIZONS:
        for H in HISTORY_LENGTHS:
            # population OLS, leave-one-patient-out, once per held-out patient
            pop = {}
            for cov in (False, True):
                for held in PATIENTS:
                    others = [train[(p, h, H)] for p in PATIENTS if p != held]
                    pop[(held, cov)] = fit_ols(others, cov)
            for pid in PATIENTS:
                ws, tr = test[(pid, h, H)], train[(pid, h, H)]
                base = dict(patient=pid, cohort=COHORT_OF[pid], horizon_min=h * 5, horizon_steps=h,
                            window_set_H=H, n_dropped_nan_covariate=0, n_train_windows="", n_train_dropped_nan="",
                            eval_variant=eval_variant, n_fallback=int(ws.fallback.sum()))
                y = ws.target
                rows.append(_score(dict(base, method="p0", variant="persistence", k="", fit="none", covariates=False),
                                   y, persistence(ws)))
                if H == min(HISTORY_LENGTHS):
                    for k in K_VALUES:
                        rows.append(_score(dict(base, method="l1", variant=f"linear extrapolation k={k}", k=k,
                                                fit="none", covariates=False), y, linear_extrapolation(ws, k)))
                for cov in (False, True):
                    keep = covariate_ok(ws) if cov else np.ones(len(ws), bool)
                    X = design(ws, cov)[keep]
                    n_tr_ok = int(covariate_ok(tr).sum()) if cov else len(tr)
                    info = dict(base, k="", covariates=cov, n_dropped_nan_covariate=int((~keep).sum()),
                                n_train_dropped_nan=len(tr) - n_tr_ok)
                    beta = fit_ols([tr], cov)
                    rows.append(_score(dict(info, method="l2", variant=f"linear AR H={H}", fit="per-patient",
                                            n_train_windows=n_tr_ok), y[keep], X @ beta))
                    n_pop = sum(int(covariate_ok(train[(p, h, H)]).sum()) if cov else len(train[(p, h, H)])
                                for p in PATIENTS if p != pid)
                    rows.append(_score(dict(info, method="l2", variant=f"linear AR H={H}", fit="population (LOPO)",
                                            n_train_windows=n_pop), y[keep], X @ pop[(pid, cov)]))
            print(f"  h={h} H={H} done")
    return rows


def challenge_scores(rows: list[dict]) -> list[dict]:
    """The BGLP 2020 scoring (D-027): per-patient RMSE and MAE over all that patient's points,
    the mean over the six contributors per cohort, and the four-score sum
    RMSE30 + MAE30 + RMSE60 + MAE60 of those means. Computed for both evaluation variants."""
    keys = sorted({(r["method"], r["variant"], r["k"], r["fit"], r["covariates"], r["window_set_H"], r["eval_variant"])
                   for r in rows}, key=lambda t: (t[0], str(t[2]), t[5], t[3], t[4], t[6]))
    out = []
    for method, variant, k, fit, cov, H, ev in keys:
        for cohort in ("2018", "2020"):
            g = [r for r in rows if (r["method"], r["variant"], r["k"], r["fit"], r["covariates"], r["window_set_H"],
                                     r["eval_variant"], r["cohort"]) == (method, variant, k, fit, cov, H, ev, cohort)]
            m = {h: dict(rmse=np.mean([r["rmse"] for r in g if r["horizon_steps"] == h]),
                         mae=np.mean([r["mae"] for r in g if r["horizon_steps"] == h]),
                         n=sum(r["n_scored"] for r in g if r["horizon_steps"] == h),
                         fb=sum(r["n_fallback"] for r in g if r["horizon_steps"] == h)) for h in HORIZONS}
            out.append(dict(method=method, variant=variant, k=k, fit=fit, covariates=cov, window_set_H=H,
                            eval_variant=ev, cohort=cohort, n_patients=len({r["patient"] for r in g}),
                            rmse30_mean_of_six=round(m[6]["rmse"], 3), mae30_mean_of_six=round(m[6]["mae"], 3),
                            rmse60_mean_of_six=round(m[12]["rmse"], 3), mae60_mean_of_six=round(m[12]["mae"], 3),
                            four_score_sum=round(m[6]["rmse"] + m[6]["mae"] + m[12]["rmse"] + m[12]["mae"], 3),
                            n_scored_30=m[6]["n"], n_scored_60=m[12]["n"],
                            n_fallback_30=m[6]["fb"], n_fallback_60=m[12]["fb"]))
    return out


def aggregate(rows: list[dict]) -> list[dict]:
    """Per cohort and pooled: mean of per-patient metrics (the stated aggregation), the
    pooled figure beside it, SD and worst patient, and total scored n."""
    keys = sorted({(r["method"], r["variant"], r["k"], r["fit"], r["covariates"], r["horizon_steps"], r["window_set_H"])
                   for r in rows}, key=lambda t: (t[0], str(t[2]), t[6], t[3], t[4], t[5]))
    out = []
    for method, variant, k, fit, cov, h, H in keys:
        grp = [r for r in rows if (r["method"], r["variant"], r["k"], r["fit"], r["covariates"],
                                   r["horizon_steps"], r["window_set_H"]) == (method, variant, k, fit, cov, h, H)]
        for cohort in ("2018", "2020", "pooled"):
            g = [r for r in grp if cohort == "pooled" or r["cohort"] == cohort]
            n = sum(r["n_scored"] for r in g)
            worst = max(g, key=lambda r: r["rmse"])
            out.append(dict(
                method=method, variant=variant, k=k, fit=fit, covariates=cov, horizon_min=h * 5,
                window_set_H=H, eval_variant=grp[0]["eval_variant"], cohort=cohort, n_patients=len(g), n_scored=n,
                rmse_mean_of_patients=round(np.mean([r["rmse"] for r in g]), 3),
                rmse_pooled=round(float(np.sqrt(sum(r["rmse"] ** 2 * r["n_scored"] for r in g) / n)), 3),
                mae_mean_of_patients=round(np.mean([r["mae"] for r in g]), 3),
                mae_pooled=round(sum(r["mae"] * r["n_scored"] for r in g) / n, 3),
                mape_mean_of_patients=round(np.mean([r["mape"] for r in g]), 3),
                mape_pooled=round(sum(r["mape"] * r["n_scored"] for r in g) / n, 3),
                rmse_sd_across_patients=round(float(np.std([r["rmse"] for r in g], ddof=1)), 3),
                worst_patient=worst["patient"], worst_patient_rmse=worst["rmse"],
                n_dropped_nan_covariate=sum(r["n_dropped_nan_covariate"] for r in g),
                n_fallback=sum(r["n_fallback"] for r in g),
                published_rmse="", published_mae="", published_source="", published_inputs="",
                like_for_like=""))
    return out


# --------------------------------------------------------------- sanity checks

def sanity_checks(rows: list[dict]) -> list[tuple[str, bool, str]]:
    out = []
    p0 = {(r["patient"], r["horizon_steps"], r["window_set_H"]): r["rmse"] for r in rows if r["method"] == "p0"}
    # 1. persistence at 30 min lands where the literature puts it
    v = [p0[(p, 6, 6)] for p in PATIENTS]
    out.append(("p0 RMSE at 30 min within [10, 40] mg/dL for every patient",
                all(10 <= x <= 40 for x in v), f"range {min(v):.2f}-{max(v):.2f}"))
    # 2. 60-minute error exceeds 30-minute error, every baseline, every patient
    key = lambda r: (r["patient"], r["method"], r["variant"], r["fit"], r["covariates"], r["window_set_H"])  # noqa: E731
    by = {}
    for r in rows:
        by.setdefault(key(r), {})[r["horizon_steps"]] = r["rmse"]
    bad = [k for k, d in by.items() if not d[12] > d[6]]
    out.append(("60-min RMSE > 30-min RMSE for every baseline and patient", not bad,
                f"{len(by)} comparisons, {len(bad)} violations" + (f": {bad[:3]}" if bad else "")))
    # 3. two-point extrapolation an hour out is noisier than persistence for most patients
    l1 = {r["patient"]: r["rmse"] for r in rows if r["method"] == "l1" and r["k"] == 2 and r["horizon_steps"] == 12}
    worse = sum(l1[p] > p0[(p, 12, 6)] for p in PATIENTS)
    out.append(("l1 k=2 noisier than p0 at 60 min for most patients", worse > len(PATIENTS) / 2,
                f"{worse}/{len(PATIENTS)} patients"))
    return out


def main() -> None:
    rows = per_patient_rows("primary")
    agg = aggregate(rows)
    rows_b = per_patient_rows("bglp")
    agg_b = aggregate(rows_b)
    chall = challenge_scores(rows + rows_b)
    res = REPO_ROOT / "results"
    for data, fn in ((rows, "baselines_per_patient.csv"), (agg, "baselines.csv"),
                     (rows_b, "baselines_bglp_per_patient.csv"), (agg_b, "baselines_bglp.csv"),
                     (chall, "baselines_challenge_scores.csv")):
        with open(res / fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
        print(f"wrote results/{fn}")

    print(f"\n{'method':<34}{'h':>4}{'H':>4}{'cohort':<8}{'n':>7}{'RMSE mean':>10}{'pooled':>8}"
          f"{'MAE':>8}{'MAPE':>7}{'SD':>7}{'worst':>12}")
    for r in agg:
        label = f"{r['method']} {r['variant']}" + (f" {r['fit']}" if r["method"] == "l2" else "") \
                + (" +cov" if r["covariates"] else "")
        print(f"{label:<34}{r['horizon_min']:>4}{r['window_set_H']:>4}{r['cohort']:<8}{r['n_scored']:>7}"
              f"{r['rmse_mean_of_patients']:>10.2f}{r['rmse_pooled']:>8.2f}{r['mae_mean_of_patients']:>8.2f}"
              f"{r['mape_mean_of_patients']:>7.1f}{r['rmse_sd_across_patients']:>7.2f}"
              f"{r['worst_patient'] + ' ' + str(r['worst_patient_rmse']):>12}")

    print("\nprimary vs conforming (bglp): mean of six per cohort; sum = RMSE30+MAE30+RMSE60+MAE60")
    print(f"{'method':<42}{'cohort':<7}{'variant':<9}{'RMSE30':>8}{'MAE30':>8}{'RMSE60':>8}{'MAE60':>8}{'sum':>8}{'n30':>7}{'fb30':>6}{'fb60':>6}")
    for r in chall:
        label = f"{r['method']} {r['fit']}" + (f" k={r['k']}" if r["k"] != "" else f" H={r['window_set_H']}") + (" +cov" if r["covariates"] else "")
        print(f"{label:<42}{r['cohort']:<7}{r['eval_variant']:<9}{r['rmse30_mean_of_six']:>8.2f}{r['mae30_mean_of_six']:>8.2f}"
              f"{r['rmse60_mean_of_six']:>8.2f}{r['mae60_mean_of_six']:>8.2f}{r['four_score_sum']:>8.2f}{r['n_scored_30']:>7}"
              f"{r['n_fallback_30']:>6}{r['n_fallback_60']:>6}")
    print("\nsanity checks (primary, then bglp)")
    checks = sanity_checks(rows) + sanity_checks(rows_b)
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name} - {detail}")
    if not all(ok for _, ok, _ in checks):
        print("a sanity check failed: the numbers above are not to be trusted or committed")
        sys.exit(1)


if __name__ == "__main__":
    main()
