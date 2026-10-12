"""Week 5, Monday: conformal intervals around the frozen final model (D-045 to D-048).

Calibration set: the frozen validation tail of each patient's training file (consumed by
early stopping; disclosed, and its consequence measured). In-sample bracket: calibration
on the fit portion's residuals, known optimistic. Test: both evaluation variants, both
horizons, every patient, levels 80 / 90 / 95%, methods c0-c3.

The final model's predictions on the training file come from the regenerated ensemble
members (src/staircase.py regen), each proven identical to the frozen model on its test
predictions before use; the ensemble is the median over members on the targets all
members share. No test target is used for calibration or for choosing gamma.

Writes results/coverage.csv, results/onset_first_pass.csv, figures/coverage_calibration.png
and figures/interval_widths.png. Aggregate only.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.conformal import LEVELS, c0, c1, c2, c3, choose_gamma, onset_detection, summarize  # noqa: E402
from src.data.loader import COHORT_2020, COHORT_OF, PATIENTS, REPO_ROOT  # noqa: E402
from src.data.windows import HORIZONS, windows_for  # noqa: E402
from src.eval.baselines import design, fit_ols, persistence  # noqa: E402
from src.eval.metrics import rmse  # noqa: E402
from src.staircase import PRED_DIR, SEEDS5, existing_rows  # noqa: E402

REF_H = 24          # history used for c1's sigma and for the origin reading


def member_keys() -> dict[tuple[str, int], list[str]]:
    rows = existing_rows()
    out = defaultdict(list)
    for k, r in rows.items():
        if r["step"] == "step7ens" or (r["step"] == "step7" and r["seed"] == SEEDS5[0]):
            out[(r["patient"], r["horizon_steps"])].append(k)
    return out


def _median_on_common(series: list[tuple[np.ndarray, np.ndarray, np.ndarray]]):
    """series: [(ts, y, pred)] per member -> (ts, y, median pred) on the shared targets."""
    common = series[0][0]
    for ts, _, _ in series[1:]:
        common = np.intersect1d(common, ts)
    preds = np.stack([p[np.isin(ts, common)] for ts, _, p in series])
    y = series[0][1][np.isin(series[0][0], common)]
    return common, y, np.median(preds, axis=0)


def frozen_predictions(pid: str, h: int, variant: str, keys: list[str]):
    """(train-file ts, y, pred, is_val), (test ts, y, pred) for the ensemble."""
    tr, te = [], []
    for k in keys:
        v = np.load(PRED_DIR / f"{k}-val.npz")
        tr.append((v[f"{variant}_ts"], v[f"{variant}_y"], v[f"{variant}_pred"], v[f"{variant}_is_val"]))
        t = np.load(PRED_DIR / f"{k}.npz")
        ws = windows_for(pid, "test", h, int(existing_rows()[k]["H"]), variant=variant)
        te.append((t[f"{variant}_{pid}_ts"], ws.target, t[f"{variant}_{pid}_pred"]))
    ts_tr, y_tr, p_tr = _median_on_common([(a, b, c) for a, b, c, _ in tr])
    is_val = tr[0][3][np.isin(tr[0][0], ts_tr)]
    ts_te, y_te, p_te = _median_on_common(te)
    return (ts_tr, y_tr, p_tr, is_val), (ts_te, y_te, p_te)


def aligned_windows(pid: str, split: str, h: int, variant: str, ts: np.ndarray):
    """History (for sigma) and origin reading at the given targets, from the H=24 windows."""
    ws = windows_for(pid, split, h, REF_H, variant=variant)
    m = np.isin(ws.target_ts, ts)
    assert np.array_equal(ws.target_ts[m], ts[np.isin(ts, ws.target_ts)])
    keep = np.isin(ts, ws.target_ts)
    return ws.history[m], ws.history[m, -1], keep


def main() -> None:
    keys = member_keys()
    rows, onset_rows, calib_plot, width_plot = [], [], defaultdict(list), defaultdict(list)
    fired = []
    for variant in ("bglp", "primary"):
        for h in HORIZONS:
            for pid in PATIENTS:
                (ts_tr, y_tr, p_tr, is_val), (ts_te, y_te, p_te) = frozen_predictions(pid, h, variant, keys[(pid, h)])
                hist_tr, origin_tr, keep_tr = aligned_windows(pid, "train", h, variant, ts_tr)
                hist_te, origin_te, keep_te = aligned_windows(pid, "test", h, variant, ts_te)
                ts_tr, y_tr, p_tr, is_val = ts_tr[keep_tr], y_tr[keep_tr], p_tr[keep_tr], is_val[keep_tr]
                ts_te, y_te, p_te = ts_te[keep_te], y_te[keep_te], p_te[keep_te]
                fit, val = ~is_val, is_val
                assert val.sum() > 100 and fit.sum() > 100
                # baselines on the same test targets, for the onset comparison
                ws24 = windows_for(pid, "test", h, REF_H, variant=variant); tr24 = windows_for(pid, "train", h, REF_H, variant=variant)
                m24 = np.isin(ws24.target_ts, ts_te)
                p_persist, p_lin = persistence(ws24)[m24], (design(ws24, True) @ fit_ols([tr24], True))[m24]
                for alpha in (1 - l for l in LEVELS):
                    level = round(1 - alpha, 2)
                    gamma = choose_gamma(y_tr[fit], p_tr[fit], y_tr[val], p_tr[val], alpha)
                    for cal_name, cm in (("validation", val), ("train (in-sample bracket)", fit)):
                        ivs = {"c0": c0(y_tr[cm], p_tr[cm], p_te, alpha),
                               "c1": c1(y_tr[cm], p_tr[cm], hist_tr[cm], p_te, hist_te, alpha),
                               "c2": c2(y_tr[cm], p_tr[cm], p_te, alpha),
                               "c3": c3(y_tr[cm], p_tr[cm], y_te, p_te, alpha, gamma)}
                        for method, iv in ivs.items():
                            s = summarize(y_te, iv)
                            rows.append(dict(method=method, level=level, calibration=cal_name, variant=variant, horizon_min=h * 5,
                                             patient=pid, cohort=COHORT_OF[pid], n_cal=int(cm.sum()), gamma=gamma if method == "c3" else "",
                                             point_rmse=round(rmse(y_te, p_te), 3), **s))
                            if cal_name == "validation":
                                calib_plot[(method, variant, h)].append((level, s["coverage"]))
                                if level == 0.90 and variant == "bglp":
                                    width_plot[(method, h)].append(iv.width[np.isfinite(iv.width)])
                            if method == "c0" and cal_name == "validation" and level == 0.90 and s["coverage"] < 0.85:
                                fired.append((variant, h, pid, s["coverage"]))
                        # onset detection, c0 at 90% (today's first pass) beside the point forecasts
                        if cal_name == "validation" and level == 0.90:
                            for thr in (70.0, 54.0):
                                for name, lower in (("final model, point", p_te), ("c0 lower bound, 90%", ivs["c0"].lower),
                                                    ("linear AR, point", p_lin), ("persistence", p_persist)):
                                    onset_rows.append(dict(detector=name, variant=variant, horizon_min=h * 5, patient=pid,
                                                           cohort=COHORT_OF[pid], **onset_detection(y_te, origin_te, lower, thr)))
                print(f"  {variant} h={h} {pid}: n_cal {int(val.sum())}, c0@90 coverage {next(r['coverage'] for r in rows if r['patient']==pid and r['variant']==variant and r['horizon_min']==h*5 and r['method']=='c0' and r['level']==0.9 and r['calibration']=='validation'):.3f}", flush=True)

    # cohort and pooled rows: coverage pooled over points (weighted by n), widths as mean of patient medians
    agg = []
    for key in sorted({(r["method"], r["level"], r["calibration"], r["variant"], r["horizon_min"]) for r in rows}):
        g = [r for r in rows if (r["method"], r["level"], r["calibration"], r["variant"], r["horizon_min"]) == key]
        for cohort, pids in (("2020", COHORT_2020), ("all12", PATIENTS)):
            gg = [r for r in g if r["patient"] in pids]
            n = sum(r["n"] for r in gg)
            agg.append(dict(method=key[0], level=key[1], calibration=key[2], variant=key[3], horizon_min=key[4], patient="", cohort=cohort,
                            n_cal="", gamma="", point_rmse="", n=n,
                            coverage=round(sum(r["coverage"] * r["n"] for r in gg) / n, 4),
                            mean_width=round(float(np.mean([r["mean_width"] for r in gg if np.isfinite(r["mean_width"])] or [np.inf])), 2),
                            median_width=round(float(np.median([r["median_width"] for r in gg if np.isfinite(r["median_width"])] or [np.inf])), 2),
                            useless_fraction=round(sum(r["useless_fraction"] * r["n"] for r in gg) / n, 4),
                            note="; ".join(sorted({r["note"] for r in gg if r["note"]}))[:200]))
    res = REPO_ROOT / "results"
    for data, fn in ((rows + agg, "coverage.csv"), (onset_rows, "onset_first_pass.csv")):
        with open(res / fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
        print(f"wrote results/{fn}")

    # ---- figures
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    colours = {"c0": "#1f77b4", "c1": "#ff7f0e", "c2": "#2ca02c", "c3": "#9467bd"}
    for ax, h in zip(axes, HORIZONS):
        ax.plot([0.78, 0.97], [0.78, 0.97], color="black", lw=0.8, ls="--", label="nominal")
        for method in ("c0", "c1", "c2", "c3"):
            pts = [(r["level"], r["coverage"]) for r in agg if (r["method"], r["variant"], r["horizon_min"], r["cohort"], r["calibration"]) == (method, "bglp", h * 5, "2020", "validation")]
            pts.sort()
            ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", color=colours[method], label=method, lw=1.5)
            br = [(r["level"], r["coverage"]) for r in agg if (r["method"], r["variant"], r["horizon_min"], r["cohort"], r["calibration"]) == (method, "bglp", h * 5, "2020", "train (in-sample bracket)")]
            br.sort()
            ax.plot([p[0] for p in br], [p[1] for p in br], marker="x", color=colours[method], lw=0.8, ls=":", alpha=0.7)
        ax.set_title(f"{h * 5} min — conforming, 2020 cohort pooled (dotted: in-sample bracket)", fontsize=9)
        ax.set_xlabel("nominal coverage"); ax.set_xticks([0.8, 0.9, 0.95])
    axes[0].set_ylabel("empirical test coverage"); axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("Empirical coverage against nominal, calibrated on the validation tail", fontsize=10)
    fig.tight_layout(); fig.savefig(REPO_ROOT / "figures" / "coverage_calibration.png", dpi=130)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, h in zip(axes, HORIZONS):
        data = [np.concatenate(width_plot[(m, h)]) for m in ("c0", "c1", "c2", "c3")]
        ax.boxplot(data, tick_labels=["c0", "c1", "c2", "c3"], showfliers=False, widths=0.5)
        ax.axhline(100, color="#d62728", lw=0.8, ls="--"); ax.text(4.35, 102, "useless", fontsize=7, color="#d62728")
        ax.set_title(f"{h * 5} min, 90% level, conforming, 2020 — finite widths", fontsize=9)
    axes[0].set_ylabel("interval width (mg/dL)")
    fig.suptitle("Interval width by method (infinite ACI intervals excluded; see useless_fraction)", fontsize=10)
    fig.tight_layout(); fig.savefig(REPO_ROOT / "figures" / "interval_widths.png", dpi=130)
    print("wrote figures/coverage_calibration.png, figures/interval_widths.png")

    print("\ncoverage, conforming, 2020 pooled, validation-calibrated (in-sample bracket in brackets):")
    print(f"{'method':<7}{'h':>4}" + "".join(f"{f'{int(l*100)}% cov (bracket)':>23}{'width':>8}{'useless':>9}" for l in LEVELS))
    for h in HORIZONS:
        for method in ("c0", "c1", "c2", "c3"):
            line = f"{method:<7}{h * 5:>4}"
            for l in LEVELS:
                r = next(r for r in agg if (r["method"], r["level"], r["calibration"], r["variant"], r["horizon_min"], r["cohort"]) == (method, round(l, 2), "validation", "bglp", h * 5, "2020"))
                b = next(r for r in agg if (r["method"], r["level"], r["calibration"], r["variant"], r["horizon_min"], r["cohort"]) == (method, round(l, 2), "train (in-sample bracket)", "bglp", h * 5, "2020"))
                line += f"  {r['coverage']:.3f} ({b['coverage']:.3f}){r['median_width']:>8.1f}{r['useless_fraction']:>9.3f}"
            print(line)
    print("\nonset detection at 70 mg/dL, conforming, 2020 pooled, 90%:")
    for h in HORIZONS:
        for det in ("persistence", "linear AR, point", "final model, point", "c0 lower bound, 90%"):
            g = [r for r in onset_rows if r["detector"] == det and r["threshold"] == 70.0 and r["variant"] == "bglp" and r["horizon_min"] == h * 5 and r["patient"] in COHORT_2020]
            tp, fn, fp, tn = (sum(r[k] for r in g) for k in ("tp", "fn", "fp", "tn"))
            days = sum(r["n_at_risk"] for r in g) * 5 / 1440
            print(f"  h={h * 5:>3} {det:<22} onsets {tp + fn:>4}  sens {tp / (tp + fn) if tp + fn else float('nan'):.3f}  spec {tn / (tn + fp) if tn + fp else float('nan'):.3f}  "
                  f"ppv {tp / (tp + fp) if tp + fp else float('nan'):.3f}  false alarms/day {fp / days:.2f}")
    if fired:
        print(f"\nSTOP CONDITION: c0 at 90% more than 5 points below nominal on {len(fired)} patient sets: {fired[:6]}")
        sys.exit(1)
    print("\nstop condition (c0@90 within 5 points of nominal on every patient set): clear")


if __name__ == "__main__":
    main()
