"""Part D: the Week 4 clinical tables, computed for the INTERVALS - c0-c5 at 90%, conforming,
2020 cohort pooled - beside the point forecast and persistence (D-067).

  zone coverage    Clarke and Parkes zones of each bound; the share of intervals whose two
                   bounds both lie in A or B, and the share with either bound in C, D or E
  MARD of bounds   each bound treated as an estimate, beside the point forecast's
  interval TIR     for 70-180: entirely in range, entirely out, or straddling a boundary -
                   the straddle share is the honest indeterminate rate
  hypoglycaemia    onset sensitivity / specificity / PPV / false alarms per day at 70 and 54
                   for the point forecast, the linear AR and each lower bound at 90%

Writes results/interval_clinical.csv. Aggregate only. c4 has no upper bound: its zone, MARD
and TIR rows are lower-bound only and say so.
"""
from __future__ import annotations

import csv

import numpy as np

from src.conformal import c0, c1, c2, c3, c4, c5, choose_gamma, onset_detection
from src.data.loader import COHORT_2020, REPO_ROOT
from src.data.windows import HORIZONS, windows_for
from src.eval.baselines import design, fit_ols
from src.eval.clinical import clarke_zones, mard, parkes_zones
from src.eval.coverage_breakdown import patient_sets

LEVEL = 0.90


def intervals_for(sets, h):
    alpha = 1 - LEVEL
    ys, ps, origins, lins, pers = [], [], [], [], []
    lo = {m: [] for m in ("c0", "c1", "c2", "c3", "c4", "c5")}; up = {m: [] for m in lo}
    for pid in COHORT_2020:
        s = sets[(pid, h)]
        ys.append(s["y"]); ps.append(s["p"]); origins.append(s["origin"])
        ws24 = windows_for(pid, "test", h, 24, variant="bglp"); tr24 = windows_for(pid, "train", h, 24, variant="bglp")
        m24 = np.isin(ws24.target_ts, s["ts"])
        lins.append((design(ws24, True) @ fit_ols([tr24], True))[m24]); pers.append(ws24.history[m24, -1])
        ivs = {"c0": c0(s["cal_y"], s["cal_p"], s["p"], alpha),
               "c1": c1(s["cal_y"], s["cal_p"], s["cal_hist"], s["p"], s["hist"], alpha),
               "c2": c2(s["cal_y"], s["cal_p"], s["p"], alpha),
               "c3": c3(s["cal_y"], s["cal_p"], s["y"], s["p"], alpha, choose_gamma(s["fit_y"], s["fit_p"], s["cal_y"], s["cal_p"], alpha)),
               "c4": c4(s["cal_y"], s["cal_p"], s["p"], alpha),
               "c5": c5(s["cal_y"], s["cal_p"], s["cal_hist"], s["p"], s["hist"], alpha)}
        for m, iv in ivs.items():
            lo[m].append(iv.lower); up[m].append(iv.upper)
    cat = np.concatenate
    return cat(ys), cat(ps), cat(origins), cat(lins), cat(pers), {m: cat(v) for m, v in lo.items()}, {m: cat(v) for m, v in up.items()}


def main() -> None:
    sets = patient_sets("bglp")
    rows = []
    for h in HORIZONS:
        y, p, origin, lin, per, lo, up = intervals_for(sets, h)
        n = len(y)
        # ---- point estimates: zones, MARD, TIR
        for name, est in (("point forecast", p), ("linear AR, point", lin), ("persistence", per)):
            cz, pz = clarke_zones(y, est), parkes_zones(y, est)
            in_range = (est >= 70) & (est <= 180)
            rows.append(dict(section="point estimate", method=name, horizon_min=h * 5, n=n,
                             clarke_AB_pct=round(100 * np.isin(cz, ["A", "B"]).mean(), 2), clarke_CDE_pct=round(100 * np.isin(cz, ["C", "D", "E"]).mean(), 2),
                             parkes_AB_pct=round(100 * np.isin(pz, ["A", "B"]).mean(), 2), parkes_CDE_pct=round(100 * np.isin(pz, ["C", "D", "E"]).mean(), 2),
                             mard_pct=mard(y, est)["mard_pct"], tir_in_pct=round(100 * in_range.mean(), 2), tir_out_pct=round(100 * (~in_range).mean(), 2),
                             tir_straddle_pct=0.0, note=""))
        # ---- intervals
        for m in ("c0", "c1", "c2", "c3", "c4", "c5"):
            l, u = lo[m], up[m]
            fl, fu = np.isfinite(l), np.isfinite(u)
            one_sided = m == "c4"
            zl = np.full(n, "X", dtype="<U1"); zl[fl] = clarke_zones(y[fl], l[fl])
            pl = np.full(n, "X", dtype="<U1"); pl[fl] = parkes_zones(y[fl], l[fl])
            if one_sided:
                rows.append(dict(section="interval, 90%", method=f"{m} (lower bound only)", horizon_min=h * 5, n=n,
                                 clarke_AB_pct=round(100 * np.isin(zl, ["A", "B"]).mean(), 2), clarke_CDE_pct=round(100 * (~np.isin(zl, ["A", "B"])).mean(), 2),
                                 parkes_AB_pct=round(100 * np.isin(pl, ["A", "B"]).mean(), 2), parkes_CDE_pct=round(100 * (~np.isin(pl, ["A", "B"])).mean(), 2),
                                 mard_pct=mard(y[fl], l[fl])["mard_pct"], tir_in_pct="", tir_out_pct="", tir_straddle_pct="",
                                 note="one-sided: zone and MARD of the lower bound; no TIR classification without an upper bound"))
                continue
            zu = np.full(n, "X", dtype="<U1"); zu[fu] = clarke_zones(y[fu], u[fu])
            pu = np.full(n, "X", dtype="<U1"); pu[fu] = parkes_zones(y[fu], u[fu])
            both_ab_c = np.isin(zl, ["A", "B"]) & np.isin(zu, ["A", "B"]); any_cde_c = ~both_ab_c
            both_ab_p = np.isin(pl, ["A", "B"]) & np.isin(pu, ["A", "B"]); any_cde_p = ~both_ab_p
            entirely_in = (l >= 70) & (u <= 180); entirely_out = (u < 70) | (l > 180)
            straddle = ~entirely_in & ~entirely_out
            n_inf = int((~fl | ~fu).sum())
            rows.append(dict(section="interval, 90%", method=m, horizon_min=h * 5, n=n,
                             clarke_AB_pct=round(100 * both_ab_c.mean(), 2), clarke_CDE_pct=round(100 * any_cde_c.mean(), 2),
                             parkes_AB_pct=round(100 * both_ab_p.mean(), 2), parkes_CDE_pct=round(100 * any_cde_p.mean(), 2),
                             mard_pct=f"lower {mard(y[fl], l[fl])['mard_pct']:.1f} / upper {mard(y[fu], u[fu])['mard_pct']:.1f}",
                             tir_in_pct=round(100 * entirely_in.mean(), 2), tir_out_pct=round(100 * entirely_out.mean(), 2), tir_straddle_pct=round(100 * straddle.mean(), 2),
                             note=(f"{n_inf} infinite intervals ({100 * n_inf / n:.1f}%) counted as outside A+B and as straddling" if n_inf else "")))
        # ---- hypoglycaemia onset table at 90%
        for thr in (70.0, 54.0):
            for name, est in [("point forecast", p), ("linear AR, point", lin), ("persistence", per)] + [(f"{m} lower bound, 90%", lo[m]) for m in ("c0", "c1", "c2", "c3", "c4", "c5")]:
                o = onset_detection(y, origin, est, thr)
                rows.append(dict(section=f"hypoglycaemia onset at {thr:.0f}", method=name, horizon_min=h * 5, n=o["n_at_risk"],
                                 clarke_AB_pct="", clarke_CDE_pct="", parkes_AB_pct="", parkes_CDE_pct="", mard_pct="",
                                 tir_in_pct="", tir_out_pct="", tir_straddle_pct="",
                                 note=f"onsets {o['n_onsets']}; sensitivity {o['sensitivity']}; specificity {o['specificity']}; ppv {o['ppv']}; "
                                      f"false alarms/day {o['false_alarms_per_day']}" + (f"; censored: {o['n_onsets_at_sensor_floor']} at floor" if o["censored"] else ""),
                                 **{"sensitivity": o["sensitivity"], "specificity": o["specificity"], "ppv": o["ppv"],
                                    "false_alarms_per_day": o["false_alarms_per_day"], "n_onsets": o["n_onsets"], "censored": o["censored"]}))
    keys = list({k: None for r in rows for k in r})
    with open(REPO_ROOT / "results" / "interval_clinical.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"wrote results/interval_clinical.csv ({len(rows)} rows)\n")
    print(f"{'method':<24}{'h':>4}{'Clarke A+B both':>16}{'any CDE':>9}{'MARD':>22}{'TIR in':>8}{'out':>6}{'straddle':>10}")
    for r in rows:
        if r["section"] in ("point estimate", "interval, 90%"):
            print(f"{r['method']:<24}{r['horizon_min']:>4}{r['clarke_AB_pct']:>16}{r['clarke_CDE_pct']:>9}{str(r['mard_pct']):>22}{str(r['tir_in_pct']):>8}{str(r['tir_out_pct']):>6}{str(r['tir_straddle_pct']):>10}")


if __name__ == "__main__":
    main()
