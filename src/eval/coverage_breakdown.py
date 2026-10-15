"""Wednesday's coverage breakdowns, computed Tuesday and interpreted Wednesday (D-051).

For c0-c3 at 80/90/95%, both horizons, both variants: coverage, median width and useless
fraction broken down by true-glucose bin, by rate of change at the origin, by time since
the last gap over 30 minutes, and by patient. For c3: coverage and width as a time series
across 552's 118-hour gap and the ten largest test-split gaps in the cohort. Width against
error: per method and patient, the Spearman correlation between interval width and
absolute error. No interpretation here; the tables are read on Wednesday.

Bins (pure functions below, tested at their edges):
  true glucose      < 54, 54-70, 70-180, 180-250, > 250 mg/dL   (lower edge inclusive)
  rate of change    |x_origin - x_previous| / 5 min:  <= 1, 1-2, 2-3, > 3 mg/dL/min
  time since gap    minutes from the end of the most recent interval > 30 min in the
                    patient's series (training tail included) to the forecast origin:
                    < 30, 30-120, > 120 (no gap inside a 2-hour window), or no prior gap
Writes results/coverage_breakdown.csv, results/width_vs_error.csv, results/c3_after_gaps.csv,
figures/c3_after_gaps.png. Aggregate only.
"""
from __future__ import annotations

import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from src.conformal import LEVELS, c0, c1, c2, c3, choose_gamma  # noqa: E402
from src.data.loader import COHORT_2020, COHORT_OF, PATIENTS, REPO_ROOT, load_patient  # noqa: E402
from src.data.windows import HORIZONS, STEP_S, _series_for  # noqa: E402
from src.eval.coverage import aligned_windows, frozen_predictions, member_keys  # noqa: E402

GLUCOSE_EDGES = (54.0, 70.0, 180.0, 250.0)
GLUCOSE_NAMES = ("<54", "54-70", "70-180", "180-250", ">250")
ROC_EDGES = (1.0, 2.0, 3.0)
ROC_NAMES = ("<=1", "1-2", "2-3", ">3")
GAP_EDGES_MIN = (30.0, 120.0)
GAP_NAMES = ("<30 min", "30-120 min", ">120 min (no gap in 2-h window)", "no prior gap")
LONG_GAP_S = 1800
ROLL = 12          # one hour of targets for the rolling coverage


def glucose_bin(y) -> np.ndarray:
    """Index into GLUCOSE_NAMES; lower edges inclusive: 54 -> '54-70', 70 -> '70-180'."""
    return np.searchsorted(np.array(GLUCOSE_EDGES), np.asarray(y, float), side="right")


def roc_bin(origin, previous) -> np.ndarray:
    """|origin - previous| per 5-minute step, in mg/dL per minute; edges inclusive on the
    lower side of the bin above: exactly 1 -> '<=1', exactly 2 -> '1-2'."""
    r = np.abs(np.asarray(origin, float) - np.asarray(previous, float)) / (STEP_S / 60.0)
    return np.searchsorted(np.array(ROC_EDGES), r, side="left")


def gap_bin(minutes_since_gap) -> np.ndarray:
    """minutes_since_gap: NaN where no gap precedes the origin. 30 -> '30-120', 120 -> '>120'."""
    m = np.asarray(minutes_since_gap, float)
    out = np.searchsorted(np.array(GAP_EDGES_MIN), np.nan_to_num(m, nan=-1.0), side="right")
    return np.where(np.isnan(m), 3, out)


def minutes_since_long_gap(series_ts, origin_ts) -> np.ndarray:
    """For each origin time, minutes since the END of the most recent interval > 30 min in
    the series; NaN if none precedes it."""
    t = np.asarray(series_ts).astype("datetime64[s]").astype(np.int64)
    o = np.asarray(origin_ts).astype("datetime64[s]").astype(np.int64)
    d = np.diff(t)
    gap_end = t[1:][d > LONG_GAP_S]
    i = np.searchsorted(gap_end, o, side="right") - 1
    out = np.where(i >= 0, (o - gap_end[np.maximum(i, 0)]) / 60.0, np.nan)
    return out


def _intervals(s: dict, method: str, alpha: float):
    if method == "c0":
        return c0(s["cal_y"], s["cal_p"], s["p"], alpha)
    if method == "c1":
        return c1(s["cal_y"], s["cal_p"], s["cal_hist"], s["p"], s["hist"], alpha)
    if method == "c2":
        return c2(s["cal_y"], s["cal_p"], s["p"], alpha)
    g = choose_gamma(s["fit_y"], s["fit_p"], s["cal_y"], s["cal_p"], alpha)
    return c3(s["cal_y"], s["cal_p"], s["y"], s["p"], alpha, g)


def patient_sets(variant: str) -> dict:
    keys = member_keys()
    out = {}
    for h in HORIZONS:
        for pid in PATIENTS:
            (ts_tr, y_tr, p_tr, is_val), (ts_te, y_te, p_te) = frozen_predictions(pid, h, variant, keys[(pid, h)])
            hist_tr, _, keep_tr = aligned_windows(pid, "train", h, variant, ts_tr)
            hist_te, origin_te, keep_te = aligned_windows(pid, "test", h, variant, ts_te)
            y_tr, p_tr, is_val = y_tr[keep_tr], p_tr[keep_tr], is_val[keep_tr]
            ts_te, y_te, p_te = ts_te[keep_te], y_te[keep_te], p_te[keep_te]
            series, _, _ = _series_for(load_patient(pid, "test"))
            origin_ts = ts_te - np.timedelta64(h * STEP_S, "s")
            out[(pid, h)] = dict(ts=ts_te, y=y_te, p=p_te, hist=hist_te, origin=origin_te, previous=hist_te[:, -2],
                                 since_gap=minutes_since_long_gap(series.ts, origin_ts),
                                 cal_y=y_tr[is_val], cal_p=p_tr[is_val], cal_hist=hist_tr[is_val],
                                 fit_y=y_tr[~is_val], fit_p=p_tr[~is_val])
    return out


MIN_N = 100       # D-053 rule 1: a cell is reportable only with at least this many targets


def _summ(y, iv, m):
    """Coverage with its miscoverage split by direction (D-053 rule 4): miss_below is a
    target under the lower bound - the lower bound too high, the dangerous direction for
    hypoglycaemia - miss_above a target over the upper bound. covered + miss_below +
    miss_above = n for every cell."""
    if not m.any():
        return dict(n=0, coverage=float("nan"), covered=0, miss_below=0, miss_above=0, miss_below_share=float("nan"),
                    median_width=float("nan"), useless_fraction=float("nan"), reportable=False)
    w = iv.width[m]; fin = np.isfinite(w); yy = y[m]
    below, above = yy < iv.lower[m], yy > iv.upper[m]
    cov = ~below & ~above
    return dict(n=int(m.sum()), coverage=round(float(cov.mean()), 4), covered=int(cov.sum()),
                miss_below=int(below.sum()), miss_above=int(above.sum()),
                miss_below_share=round(float(below.sum() / max(1, (~cov).sum())), 4),
                median_width=round(float(np.median(w[fin])), 2) if fin.any() else float("inf"),
                useless_fraction=round(float(np.mean(~fin | (w > 100))), 4), reportable=bool(m.sum() >= MIN_N))


def main() -> None:
    rows, wve, ts_rows, rec_rows = [], [], [], []
    scatter = defaultdict(list)
    for variant in ("bglp", "primary"):
        sets = patient_sets(variant)
        for h in HORIZONS:
            for level in LEVELS:
                alpha = 1 - level
                for method in ("c0", "c1", "c2", "c3"):
                    ivs = {pid: _intervals(sets[(pid, h)], method, alpha) for pid in PATIENTS}
                    for cohort, pids in (("2020", COHORT_2020), ("all12", PATIENTS)):
                        dims = {"true_glucose": (GLUCOSE_NAMES, lambda s: glucose_bin(s["y"])),
                                "forecast_glucose": (GLUCOSE_NAMES, lambda s: glucose_bin(s["p"])),   # what Mondrian c2 conditions on
                                "rate_of_change": (ROC_NAMES, lambda s: roc_bin(s["origin"], s["previous"])),
                                "time_since_gap": (GAP_NAMES, lambda s: gap_bin(s["since_gap"]))}
                        for dim, (names, fn) in dims.items():
                            for b, name in enumerate(names):
                                y = np.concatenate([sets[(pid, h)]["y"] for pid in pids])
                                lo = np.concatenate([ivs[pid].lower for pid in pids]); up = np.concatenate([ivs[pid].upper for pid in pids])
                                m = np.concatenate([fn(sets[(pid, h)]) == b for pid in pids])
                                from src.conformal import Intervals
                                rows.append(dict(method=method, level=round(level, 2), variant=variant, horizon_min=h * 5, cohort=cohort,
                                                 dimension=dim, bin=name, **_summ(y, Intervals(lo, up), m)))
                    for pid in PATIENTS:
                        s = sets[(pid, h)]; iv = ivs[pid]
                        rows.append(dict(method=method, level=round(level, 2), variant=variant, horizon_min=h * 5, cohort=COHORT_OF[pid],
                                         dimension="patient", bin=pid, **_summ(s["y"], iv, np.ones(len(s["y"]), bool))))
                        if level == 0.90 and variant == "bglp" and pid in COHORT_2020:
                            scatter[(method, h)].append((pid, iv.width, np.abs(s["y"] - s["p"])))
                        if level == 0.90:
                            w = iv.width; fin = np.isfinite(w); err = np.abs(s["y"] - s["p"])
                            # a constant-width method varies by an ulp across windows; that is not a width
                            rho = spearmanr(w[fin], err[fin]).statistic if fin.sum() > 10 and np.ptp(w[fin]) > 1e-6 else float("nan")
                            wve.append(dict(method=method, variant=variant, horizon_min=h * 5, patient=pid, cohort=COHORT_OF[pid],
                                            n=int(fin.sum()), spearman_width_vs_abs_error=round(float(rho), 4) if rho == rho else "",
                                            width_sd=round(float(np.std(w[fin])), 2) if fin.any() else ""))
                print(f"  {variant} h={h} level {level:.2f} done", flush=True)
        # c3 time series across gaps (90%, h=6 and 12), conforming only
        if variant == "bglp":
            gaps = []
            for pid in PATIENTS:
                # gaps inside the TEST file only (the time-since-gap feature above uses the
                # training tail too; a training gap has no test targets around it to plot)
                t = load_patient(pid, "test").cgm.ts.astype("datetime64[s]").astype(np.int64)
                d = np.diff(t)
                for i in np.flatnonzero(d > LONG_GAP_S):
                    gaps.append((d[i] / 3600.0, pid, t[i], t[i + 1]))
            gaps.sort(reverse=True)
            top = gaps[:10]
            for h in HORIZONS:
                for rank, (hours, pid, g0, g1) in enumerate(top, 1):
                    s = sets[(pid, h)]; iv = _intervals(s, "c3", 0.10)
                    ts = s["ts"].astype("datetime64[s]").astype(np.int64)
                    hit = ((iv.lower <= s["y"]) & (s["y"] <= iv.upper)).astype(float)
                    roll = np.convolve(hit, np.ones(ROLL) / ROLL, mode="same")
                    win = (ts >= g0 - 24 * 3600) & (ts <= g1 + 24 * 3600)
                    for k in np.flatnonzero(win):
                        ts_rows.append(dict(gap_rank=rank, patient=pid, gap_hours=round(hours, 1), horizon_min=h * 5,
                                            hours_from_gap_end=round((ts[k] - g1) / 3600.0, 3), rolling_coverage_1h=round(float(roll[k]), 3),
                                            width=round(float(iv.width[k]), 2) if np.isfinite(iv.width[k]) else "inf"))
            # recovery after each gap (B4): minutes until rolling 1-h coverage is back within 5 points
            # of nominal, and the infinite-interval fraction in the first hour after the gap
            for h in HORIZONS:
                for rank, (hours, pid, g0, g1) in enumerate(top, 1):
                    g = [r for r in ts_rows if r["gap_rank"] == rank and r["horizon_min"] == h * 5 and r["hours_from_gap_end"] >= 0]
                    g.sort(key=lambda r: r["hours_from_gap_end"])
                    rec = next((r["hours_from_gap_end"] * 60 for r in g if r["rolling_coverage_1h"] >= 0.85), float("nan"))
                    first_hour = [r for r in g if r["hours_from_gap_end"] <= 1.0]
                    inf_share = np.mean([r["width"] == "inf" for r in first_hour]) if first_hour else float("nan")
                    fin_w = [r["width"] for r in first_hour if r["width"] != "inf"]
                    rec_rows.append(dict(gap_rank=rank, patient=pid, gap_hours=round(hours, 1), horizon_min=h * 5,
                                         minutes_to_recover_within_5pts=round(rec, 1) if rec == rec else "never within 24 h",
                                         infinite_share_first_hour=round(float(inf_share), 3) if inf_share == inf_share else "",
                                         median_finite_width_first_hour=round(float(np.median(fin_w)), 1) if fin_w else "",
                                         n_targets_first_hour=len(first_hour)))
            # figure: rolling coverage and width around each gap, h=30 min
            fig, axes = plt.subplots(5, 2, figsize=(12, 13), sharex=False)
            for ax, (rank, (hours, pid, g0, g1)) in zip(axes.ravel(), enumerate(top, 1)):
                g = [r for r in ts_rows if r["gap_rank"] == rank and r["horizon_min"] == 30]
                x = [r["hours_from_gap_end"] for r in g]
                ax.plot(x, [r["rolling_coverage_1h"] for r in g], color="#1f77b4", lw=1.2, label="rolling 1-h coverage")
                ax2 = ax.twinx() if False else None  # single axis by rule: width drawn in a second panel row instead
                ax.axvspan(-hours, 0, color="#dddddd", alpha=0.6, label=f"gap {hours:.0f} h")
                ax.axhline(0.9, color="black", lw=0.6, ls="--")
                ax.set_ylim(0, 1.02); ax.set_xlim(-24 - min(hours, 48), 24)
                ax.set_title(f"#{rank} patient {pid}, gap {hours:.1f} h — c3 at 90%, 30 min, conforming", fontsize=8)
                ax.set_xlabel("hours from gap end (gap shaded; left edge clipped at 48 h)", fontsize=7); ax.tick_params(labelsize=7)
            axes[0, 0].legend(fontsize=7, frameon=False)
            fig.suptitle("c3 rolling coverage around the ten largest test-split gaps (widths in results/c3_after_gaps.csv)", fontsize=10)
            fig.tight_layout(); fig.savefig(REPO_ROOT / "figures" / "c3_after_gaps.png", dpi=120)
    res = REPO_ROOT / "results"
    # ---- figures for the README section (D-058): true-glucose bin, rate of change, width vs error
    colours = {"c0": "#1f77b4", "c1": "#ff7f0e", "c2": "#2ca02c", "c3": "#9467bd"}
    def cell(method, h, dim, b):
        return next(r for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["cohort"], r["dimension"], r["bin"]) == (method, 0.9, "bglp", h * 5, "2020", dim, b))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, h in zip(axes, HORIZONS):
        x = np.arange(len(GLUCOSE_NAMES)); wdt = 0.26
        for i, m in enumerate(("c0", "c2", "c3")):
            cells = [cell(m, h, "true_glucose", b) for b in GLUCOSE_NAMES]
            cols = [colours[m] if c["reportable"] else "#cccccc" for c in cells]
            ax.bar(x + (i - 1) * wdt, [c["coverage"] for c in cells], wdt, color=cols, edgecolor="white", label=m)
            for xi, c in zip(x + (i - 1) * wdt, cells):
                ax.text(xi, 0.02, f"n={c['n']}", rotation=90, ha="center", va="bottom", fontsize=6, color="white" if c["reportable"] else "#555555")
        ax.axhline(0.9, color="black", lw=0.8, ls="--"); ax.axhspan(0.85, 0.95, color="black", alpha=0.05)
        ax.set_xticks(x, GLUCOSE_NAMES); ax.set_ylim(0, 1.0); ax.set_title(f"{h * 5} min", fontsize=10); ax.set_xlabel("true glucose (mg/dL)")
    axes[0].set_ylabel("empirical coverage at the 90% level"); axes[1].legend(fontsize=8, frameon=False, loc="upper right")
    fig.suptitle("Coverage by true-glucose bin — c0, c2 (Mondrian by forecast bin), c3; conforming, 2020 cohort. Grey = fewer than 100 targets, not cited.", fontsize=9)
    fig.tight_layout(); fig.savefig(REPO_ROOT / "figures" / "coverage_by_glucose.png", dpi=130)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, h in zip(axes, HORIZONS):
        x = np.arange(len(ROC_NAMES)); wdt = 0.36
        for i, m in enumerate(("c0", "c1")):
            cells = [cell(m, h, "rate_of_change", b) for b in ROC_NAMES]
            cols = [colours[m] if c["reportable"] else "#cccccc" for c in cells]
            ax.bar(x + (i - 0.5) * wdt, [c["coverage"] for c in cells], wdt, color=cols, edgecolor="white", label=m)
            for xi, c in zip(x + (i - 0.5) * wdt, cells):
                ax.text(xi, 0.02, f"n={c['n']}  w={c['median_width']:.0f}  useless {c['useless_fraction']:.0%}", rotation=90, ha="center", va="bottom", fontsize=6,
                        color="white" if c["reportable"] else "#555555")
        ax.axhline(0.9, color="black", lw=0.8, ls="--"); ax.axhspan(0.85, 0.95, color="black", alpha=0.05)
        ax.set_xticks(x, ROC_NAMES); ax.set_ylim(0, 1.0); ax.set_title(f"{h * 5} min", fontsize=10); ax.set_xlabel("|rate of change| at the origin (mg/dL/min)")
    axes[0].set_ylabel("empirical coverage at the 90% level"); axes[1].legend(fontsize=8, frameon=False, loc="upper right")
    fig.suptitle("Coverage by rate of change — constant-width c0 against volatility-scaled c1; conforming, 2020. Bars carry n, median width, useless fraction.", fontsize=9)
    fig.tight_layout(); fig.savefig(REPO_ROOT / "figures" / "coverage_by_rate.png", dpi=130)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, (m, h) in zip(axes.ravel(), [(mm, hh) for hh in HORIZONS for mm in ("c1", "c3")]):
        rhos = []
        for pid, w, e in scatter[(m, h)]:
            fin = np.isfinite(w)
            rho = spearmanr(w[fin], e[fin]).statistic if fin.sum() > 10 and np.ptp(w[fin]) > 1e-6 else float("nan")
            rhos.append(f"{pid} ρ={rho:.2f}")
            ax.scatter(w[fin], e[fin], s=3, alpha=0.25, lw=0, label=None)
        ax.set_xlim(0, 250); ax.set_ylim(0, 150)
        ax.set_title(f"{m}, {h * 5} min\nper-patient Spearman: " + ", ".join(rhos), fontsize=7.5)
        ax.set_xlabel("interval width (mg/dL, finite only)"); ax.set_ylabel("|error| (mg/dL)")
    fig.suptitle("Interval width against absolute error, 90% level, conforming, 2020 cohort (c0 omitted: constant width within patient)", fontsize=9)
    fig.tight_layout(); fig.savefig(REPO_ROOT / "figures" / "width_vs_error.png", dpi=130)
    print("wrote figures/coverage_by_glucose.png, coverage_by_rate.png, width_vs_error.png")
    for data, fn in ((rows, "coverage_breakdown.csv"), (wve, "width_vs_error.csv"), (ts_rows, "c3_after_gaps.csv"), (rec_rows, "c3_recovery.csv")):
        with open(res / fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
        print(f"wrote results/{fn} ({len(data)} rows)")
    print("wrote figures/c3_after_gaps.png")


if __name__ == "__main__":
    main()
