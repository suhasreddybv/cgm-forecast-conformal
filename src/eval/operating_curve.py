"""The operating curve (D-050): hypoglycaemia onset sensitivity against false alarms per
day, for the conformal lower bounds across nominal levels and for the control that any
constant-width lower bound reduces to - the point forecast minus an offset.

c0's lower bound at any level is the point forecast minus a per-patient constant, so an
alarm on it is an alarm on the point forecast with an offset. The curve's baseline is
therefore the sweep "alarm when point forecast < threshold + delta", delta from 0 to 60
mg/dL in 2-mg/dL steps, for the final model and for the linear AR. The conformal curves
are c0-c4's lower bounds at nominal levels 5-95% in 5-point steps (c4 one-sided, D-062). Persistence is
excluded: its onset sensitivity is zero by construction (D-047).

The claim is read off the figure at 1 and 2 false alarms per day: c1/c3 above the sweep
means adaptive width detects onsets a fixed offset cannot; coinciding means conformal
supplied a principled offset and nothing more. Either is stated as found.

    python -m src.eval.operating_curve [--out-dir DIR]

Writes results/operating_curve.csv and figures/operating_curve.png (70 mg/dL, both
horizons, conforming, 2020 cohort pooled; 54 mg/dL as a second figure with the censoring
flag). Aggregate only. Scheduled for Saturday 17 Oct; defined on 12 Oct.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.conformal import c0, c1, c2, c3, c4, c5, choose_gamma, onset_detection  # noqa: E402
from src.data.loader import COHORT_2020, REPO_ROOT  # noqa: E402
from src.data.windows import HORIZONS  # noqa: E402
from src.eval.coverage import aligned_windows, frozen_predictions, member_keys  # noqa: E402
from src.data.windows import windows_for  # noqa: E402
from src.eval.baselines import design, fit_ols  # noqa: E402

DELTAS = np.arange(0, 61, 2, dtype=float)
# 5-95% in 5-point steps: a dry run on 12 Oct showed that at 70 mg/dL no lower bound at a level of
# 50% or above reaches 1-2 false alarms per day, so the grid extends down to where the control lives
LEVELS_FINE = np.round(np.arange(0.05, 0.951, 0.05), 2)
THRESHOLDS = (70.0, 54.0)
FA_MARKS = (1.0, 2.0)


N_BOOT = 2000


def pooled_onset(parts: list[dict]) -> dict:
    tp, fn, fp, tn = (sum(p[k] for p in parts) for k in ("tp", "fn", "fp", "tn"))
    days = sum(p["n_at_risk"] for p in parts) * 5 / 1440
    floor = sum(p["n_onsets_at_sensor_floor"] for p in parts)
    return dict(n_onsets=tp + fn, tp=tp, fp=fp, sensitivity=tp / (tp + fn) if tp + fn else float("nan"),
                false_alarms_per_day=fp / days if days else float("nan"), n_onsets_at_sensor_floor=floor, censored=floor > 0)


def curves(variant: str = "bglp", patients=COHORT_2020) -> list[dict]:
    """One row per (threshold, horizon, curve, parameter): pooled onset sensitivity and false
    alarms per day over the given patients."""
    keys = member_keys()
    per = {}      # (pid, h) -> dict of arrays
    for h in HORIZONS:
        for pid in patients:
            (ts_tr, y_tr, p_tr, is_val), (ts_te, y_te, p_te) = frozen_predictions(pid, h, variant, keys[(pid, h)])
            hist_tr, _, keep_tr = aligned_windows(pid, "train", h, variant, ts_tr)
            hist_te, origin_te, keep_te = aligned_windows(pid, "test", h, variant, ts_te)
            y_tr, p_tr, is_val = y_tr[keep_tr], p_tr[keep_tr], is_val[keep_tr]
            ts_te, y_te, p_te = ts_te[keep_te], y_te[keep_te], p_te[keep_te]
            ws24, tr24 = windows_for(pid, "test", h, 24, variant=variant), windows_for(pid, "train", h, 24, variant=variant)
            m24 = np.isin(ws24.target_ts, ts_te)
            p_lin = (design(ws24, True) @ fit_ols([tr24], True))[m24]
            per[(pid, h)] = dict(y=y_te, p=p_te, p_lin=p_lin, origin=origin_te, hist=hist_te,
                                 cal_y=y_tr[is_val], cal_p=p_tr[is_val], cal_hist=hist_tr[is_val],
                                 fit_y=y_tr[~is_val], fit_p=p_tr[~is_val])
    rows = []
    for thr in THRESHOLDS:
        for h in HORIZONS:
            sets = [per[(pid, h)] for pid in patients]
            for name, key in (("point forecast − δ", "p"), ("linear AR − δ", "p_lin")):
                for d in DELTAS:
                    parts = [onset_detection(s["y"], s["origin"], s[key] - d, thr) for s in sets]
                    rows.append(dict(threshold=thr, horizon_min=h * 5, curve=name, parameter=float(d), **pooled_onset(parts), _parts=parts))
            for method in ("c0", "c1", "c2", "c3", "c4", "c5"):
                for level in LEVELS_FINE:
                    alpha = 1 - level
                    parts = []
                    for s in sets:
                        if method == "c0":
                            iv = c0(s["cal_y"], s["cal_p"], s["p"], alpha)
                        elif method == "c4":
                            iv = c4(s["cal_y"], s["cal_p"], s["p"], alpha)
                        elif method == "c5":
                            iv = c5(s["cal_y"], s["cal_p"], s["cal_hist"], s["p"], s["hist"], alpha)
                        elif method == "c1":
                            iv = c1(s["cal_y"], s["cal_p"], s["cal_hist"], s["p"], s["hist"], alpha)
                        elif method == "c2":
                            iv = c2(s["cal_y"], s["cal_p"], s["p"], alpha)
                        else:
                            g = choose_gamma(s["fit_y"], s["fit_p"], s["cal_y"], s["cal_p"], alpha)
                            iv = c3(s["cal_y"], s["cal_p"], s["y"], s["p"], alpha, g)
                        parts.append(onset_detection(s["y"], s["origin"], iv.lower, thr))
                    rows.append(dict(threshold=thr, horizon_min=h * 5, curve=f"{method} lower bound" + (" (one-sided)" if method == "c4" else ""),
                                     parameter=float(level), **pooled_onset(parts), _parts=parts))
    return rows


def sensitivity_at(rows: list[dict], thr: float, h: int, curve: str, fa: float) -> float:
    """Sensitivity on a curve at a false-alarm rate, linearly interpolated along the curve;
    NaN if the curve never reaches that rate."""
    pts = sorted((r["false_alarms_per_day"], r["sensitivity"]) for r in rows
                 if r["threshold"] == thr and r["horizon_min"] == h and r["curve"] == curve and np.isfinite(r["false_alarms_per_day"]))
    xs, ys = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    if len(xs) < 2 or fa < xs.min() or fa > xs.max():
        return float("nan")
    return float(np.interp(fa, xs, ys))


def _count_array(rows: list[dict], thr: float, h: int, curve: str) -> np.ndarray:
    """[n_settings, n_patients, 4] of (tp, fn, fp, n_at_risk) for one curve, settings in
    row order."""
    sel = [r for r in rows if r["threshold"] == thr and r["horizon_min"] == h and r["curve"] == curve]
    return np.array([[[q["tp"], q["fn"], q["fp"], q["n_at_risk"]] for q in r["_parts"]] for r in sel], float)


def _curve_from_counts(a: np.ndarray, idx) -> tuple[np.ndarray, np.ndarray]:
    """(false alarms/day, sensitivity) re-pooled over the patients in idx, sorted by rate."""
    c = a[:, idx, :].sum(axis=1)                       # [n_settings, 4]
    tp, fn, fp, risk = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        sens = tp / (tp + fn)
        fa = fp / (risk * 5 / 1440)
    ok = np.isfinite(sens) & np.isfinite(fa)
    order = np.argsort(fa[ok])
    return fa[ok][order], sens[ok][order]


def sens_at(xs, ys, fa):
    if len(xs) < 2 or fa < xs.min() or fa > xs.max():
        return float("nan")
    return float(np.interp(fa, xs, ys))


def claim_rule(rows: list[dict], n_patients: int, seed: int = 0) -> list[dict]:
    """D-064: each method's sensitivity minus the point-forecast sweep's at 1 and 2 false
    alarms per day, bootstrapped over patients (numpy sums over a precomputed count array)."""
    rng = np.random.default_rng(seed)
    out = []
    methods = sorted({r["curve"] for r in rows if "lower bound" in r["curve"]} | {"linear AR − δ"})
    for thr in THRESHOLDS:
        for h in (30, 60):
            ref = _count_array(rows, thr, h, "point forecast − δ")
            full = np.arange(n_patients)
            ref_x, ref_y = _curve_from_counts(ref, full)
            boots = [rng.integers(0, n_patients, n_patients) for _ in range(N_BOOT)]
            ref_b = {fa: np.array([sens_at(*_curve_from_counts(ref, b), fa) for b in boots]) for fa in FA_MARKS}
            for m in methods:
                a = _count_array(rows, thr, h, m)
                mx, my = _curve_from_counts(a, full)
                for fa in FA_MARKS:
                    point = sens_at(mx, my, fa) - sens_at(ref_x, ref_y, fa)
                    mb = np.array([sens_at(*_curve_from_counts(a, b), fa) for b in boots])
                    diffs = mb - ref_b[fa]; fin = diffs[np.isfinite(diffs)]
                    lo, hi = (np.percentile(fin, [2.5, 97.5]) if len(fin) > 100 else (float("nan"), float("nan")))
                    verdict = ("not reached" if not np.isfinite(point) else "above the sweep" if lo > 0 else "below the sweep" if hi < 0 else "within variability")
                    out.append(dict(threshold=thr, horizon_min=h, method=m, false_alarms_per_day=fa,
                                    sensitivity=round(sens_at(mx, my, fa), 4), sweep_sensitivity=round(sens_at(ref_x, ref_y, fa), 4),
                                    difference=round(point, 4) if np.isfinite(point) else "", ci_low=round(float(lo), 4) if np.isfinite(lo) else "",
                                    ci_high=round(float(hi), 4) if np.isfinite(hi) else "", n_bootstrap_resolved=int(len(fin)), verdict=verdict))
    return out


def operating_points(rows: list[dict], n_patients: int) -> list[dict]:
    """C5: at 1 and 2 false alarms per day, every method's sensitivity, specificity, PPV and
    the setting that achieves it - the last setting on the curve at or below the rate -
    pooled and per patient."""
    out = []
    curves = sorted({r["curve"] for r in rows})
    for thr in THRESHOLDS:
        for h in (30, 60):
            for m in curves:
                sel = [r for r in rows if r["threshold"] == thr and r["horizon_min"] == h and r["curve"] == m]
                a = np.array([[[q["tp"], q["fn"], q["fp"], q["tn"], q["n_at_risk"]] for q in r["_parts"]] for r in sel], float)
                params = [r["parameter"] for r in sel]
                for scope, idx in [("pooled", np.arange(n_patients))] + [(f"patient {i}", np.array([i])) for i in range(n_patients)]:
                    c = a[:, idx, :].sum(axis=1)
                    tp, fn, fp, tn, risk = (c[:, k] for k in range(5))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        fa = fp / (risk * 5 / 1440); sens = tp / (tp + fn); spec = tn / (tn + fp); ppv = tp / (tp + fp)
                    for mark in FA_MARKS:
                        ok = np.flatnonzero(np.isfinite(fa) & (fa <= mark))
                        if len(ok) == 0:
                            out.append(dict(threshold=thr, horizon_min=h, method=m, scope=scope, false_alarms_per_day=mark, setting="not reached",
                                            sensitivity="", specificity="", ppv="", achieved_false_alarms_per_day="", n_onsets=""))
                            continue
                        k = ok[np.argmax(fa[ok])]
                        out.append(dict(threshold=thr, horizon_min=h, method=m, scope=scope, false_alarms_per_day=mark,
                                        setting=(f"δ = {params[k]:.0f} mg/dL" if "δ" in m else f"level {params[k]:.2f}"),
                                        sensitivity=round(float(sens[k]), 4) if np.isfinite(sens[k]) else "", specificity=round(float(spec[k]), 4) if np.isfinite(spec[k]) else "",
                                        ppv=round(float(ppv[k]), 4) if np.isfinite(ppv[k]) else "", achieved_false_alarms_per_day=round(float(fa[k]), 3),
                                        n_onsets=int(tp[k] + fn[k])))
    return out


def figure(rows: list[dict], thr: float, out: Path, censored: bool) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    styles = {"point forecast − δ": ("#222222", "-"), "linear AR − δ": ("#999999", "--"),
              "c0 lower bound": ("#1f77b4", "-"), "c1 lower bound": ("#ff7f0e", "-"), "c2 lower bound": ("#2ca02c", "-"), "c3 lower bound": ("#9467bd", "-"),
              "c4 lower bound (one-sided)": ("#d62728", "-"), "c5 lower bound": ("#8c564b", "-")}
    for ax, h in zip(axes, (30, 60)):
        for curve, (col, ls) in styles.items():
            pts = sorted((r["false_alarms_per_day"], r["sensitivity"]) for r in rows if r["threshold"] == thr and r["horizon_min"] == h and r["curve"] == curve)
            xs, ys = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
            s2 = sens_at(xs, ys, 2.0)
            ax.plot(xs, ys, color=col, ls=ls, lw=3.0 if curve == "point forecast − δ" else 1.5, marker="o", ms=2.5,
                    label=f"{curve} — {s2:.2f} at 2/day" if np.isfinite(s2) else f"{curve} — not reached at 2/day")
        for fa in FA_MARKS:
            ax.axvline(fa, color="black", lw=0.7, ls=":")
        ax.set_xscale("log"); ax.set_xlim(0.05, 200); ax.set_ylim(0, 1)
        ax.set_xlabel("false alarms per day (log scale)"); ax.set_title(f"{h} min", fontsize=10)
        ax.grid(True, which="both", lw=0.3, alpha=0.5)
    axes[0].set_ylabel(f"onset sensitivity at {thr:.0f} mg/dL"); axes[0].legend(fontsize=8, frameon=False, loc="lower right")
    n_on = next(r["n_onsets"] for r in rows if r["threshold"] == thr and r["horizon_min"] == 30 and r["curve"] == "c0 lower bound")
    fig.suptitle(f"Hypoglycaemia onset detection at {thr:.0f} mg/dL — conforming evaluation, 2020 cohort pooled ({n_on} onsets at 30 min)\n"
                 f"Point-forecast sweep drawn heavier as the reference. Persistence excluded: its onset sensitivity is zero by construction."
                 + ("\nReference events at the sensor floor are censored (D-023)." if censored else ""), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(out, dpi=130)


def redraw_from_csv(root: Path = REPO_ROOT) -> None:
    """Redraw both figures from results/operating_curve.csv without recomputing."""
    import csv as _csv
    rows = []
    for r in _csv.DictReader(open(root / "results" / "operating_curve.csv")):
        rows.append(dict(threshold=float(r["threshold"]), horizon_min=int(r["horizon_min"]), curve=r["curve"], parameter=float(r["parameter"]),
                         n_onsets=int(r["n_onsets"]), false_alarms_per_day=float(r["false_alarms_per_day"]), sensitivity=float(r["sensitivity"]),
                         censored=r["censored"] == "True"))
    for thr, fn in ((70.0, "operating_curve.png"), (54.0, "operating_curve_54.png")):
        figure(rows, thr, root / "figures" / fn, any(r["censored"] for r in rows if r["threshold"] == thr))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(REPO_ROOT), help="root for results/ and figures/ (default: the repository)")
    a = ap.parse_args()
    root = Path(a.out_dir)
    (root / "results").mkdir(parents=True, exist_ok=True); (root / "figures").mkdir(parents=True, exist_ok=True)
    rows = curves()
    n_pat = len(rows[0]["_parts"])
    public = [{k: v for k, v in r.items() if k != "_parts"} for r in rows]
    with open(root / "results" / "operating_curve.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(public[0])); w.writeheader(); w.writerows(public)
    claims = claim_rule(rows, n_pat)
    with open(root / "results" / "operating_curve_claim.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(claims[0])); w.writeheader(); w.writerows(claims)
    points = operating_points(rows, n_pat)
    with open(root / "results" / "operating_points.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(points[0])); w.writeheader(); w.writerows(points)
    print("\nclaim rule (D-064): sensitivity minus the point-forecast sweep at the same false-alarm rate, bootstrap over 6 patients")
    for c in claims:
        if c["threshold"] == 70.0 and c["false_alarms_per_day"] == 2.0:
            print(f"  {c['horizon_min']:>3} min {c['method']:<28} sens {c['sensitivity']}  sweep {c['sweep_sensitivity']}  diff {c['difference']} [{c['ci_low']}, {c['ci_high']}]  -> {c['verdict']}")
    for thr, fn in ((70.0, "operating_curve.png"), (54.0, "operating_curve_54.png")):
        cens = any(r["censored"] for r in rows if r["threshold"] == thr)
        figure(rows, thr, root / "figures" / fn, cens)
    print(f"wrote {root / 'results' / 'operating_curve.csv'} and the two figures")
    print(f"\n{'curve':<22}" + "".join(f"{f'{h} min @ {fa:.0f}/day':>16}" for h in (30, 60) for fa in FA_MARKS))
    for thr in THRESHOLDS:
        print(f"-- {thr:.0f} mg/dL")
        for curve in ("point forecast − δ", "linear AR − δ", "c0 lower bound", "c1 lower bound", "c2 lower bound", "c3 lower bound", "c4 lower bound (one-sided)", "c5 lower bound"):
            print(f"{curve:<22}" + "".join(f"{sensitivity_at(rows, thr, h, curve, fa):>16.3f}" for h in (30, 60) for fa in FA_MARKS))


if __name__ == "__main__":
    main()
