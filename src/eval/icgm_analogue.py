"""Apply the iCGM ANALOGUES (src/eval/icgm.py, D-059) to the frozen model's point forecast
and to each conformal interval's lower and upper bound at the 90% level, conforming
evaluation, 2020 cohort pooled; the linear AR and persistence point forecasts beside them.

Every row is an analogue of a sensor criterion (docs/icgm_criteria.md §5): a forecast or a
bound stands in for the iCGM reading, the CGM target for the laboratory reference. Writes
results/icgm_analogues.csv. Aggregate only; cells with n < 100 are marked not reportable.
"""
from __future__ import annotations

import csv
import re

import numpy as np

from src.conformal import c0, c1, c2, c3, c4, c5, choose_gamma
from src.data.loader import COHORT_2020, REPO_ROOT
from src.data.windows import HORIZONS, windows_for
from src.eval.baselines import design, fit_ols
from src.eval.coverage_breakdown import patient_sets
from src.eval.icgm import accuracy_analogue, exclusion_analogue, load_criteria, rate_analogue

LEVEL = 0.90


def estimators(sets, h):
    """name -> (reference, estimate, origin) pooled over the 2020 cohort."""
    ys, ps, origins, lins, pers = [], [], [], [], []
    lowers = {m: [] for m in ("c0", "c1", "c2", "c3", "c5")}; uppers = {m: [] for m in lowers}
    alpha = 1 - LEVEL
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
               "c5": c5(s["cal_y"], s["cal_p"], s["cal_hist"], s["p"], s["hist"], alpha)}
        for m, iv in ivs.items():
            lowers[m].append(iv.lower); uppers[m].append(iv.upper)
        lowers.setdefault("c4", []).append(c4(s["cal_y"], s["cal_p"], s["p"], alpha).lower)
    y, o = np.concatenate(ys), np.concatenate(origins)
    out = {"final model, point forecast": (y, np.concatenate(ps), o), "linear AR, point forecast": (y, np.concatenate(lins), o),
           "persistence": (y, np.concatenate(pers), o)}
    for m in lowers:
        out[f"{m} lower bound, 90%"] = (y, np.concatenate(lowers[m]), o)
        if m in uppers:
            out[f"{m} upper bound, 90%"] = (y, np.concatenate(uppers[m]), o)
    return out


def main() -> None:
    crit = load_criteria()
    sets = patient_sets("bglp")
    rows = []
    for h in HORIZONS:
        for name, (ref, est, org) in estimators(sets, h).items():
            for a in crit.accuracy:
                r = accuracy_analogue(ref, est, a)
                rows.append(dict(estimator=name, horizon_min=h * 5, criterion=f"analogue of (b)(1)(v)({a.para})",
                                 estimate_range=f"{'<' if a.lo == -np.inf and np.isfinite(a.hi) else ''}{'' if np.isfinite(a.lo) else ''}"
                                 + (f"< {a.hi:g}" if a.lo == -np.inf and np.isfinite(a.hi) else f"> {a.lo:g}" if a.hi == np.inf and np.isfinite(a.lo) else "whole range" if not np.isfinite(a.lo) else f"{a.lo:g}-{a.hi:g}"),
                                 tolerance=f"±{a.tol_value:g}{'%' if a.tol_relative else ' mg/dL'}", n=r["n"], within=r["within"],
                                 pct=r["pct"], wilson_lower_95=r["wilson_lower_95"], required_pct=a.required_pct,
                                 met="" if r["met"] is None else r["met"], reportable=r["reportable"],
                                 note="" if r["reportable"] else "n < 100: not reportable (D-053)"))
            e = exclusion_analogue(ref, est, crit)
            rows.append(dict(estimator=name, horizon_min=h * 5, criterion="analogue of (b)(1)(v)(H)", estimate_range=f"< {crit.exclusion_low:g}",
                             tolerance=f"reference > {crit.exclusion_high:g}", n=e["n"], within=e["H_violations"], pct=e["H_per_1000"],
                             wilson_lower_95="", required_pct=0, met=e["H_violations"] == 0, reportable=True, note="violations; pct column is per 1,000 estimates"))
            rows.append(dict(estimator=name, horizon_min=h * 5, criterion="analogue of (b)(1)(v)(I)", estimate_range=f"> {crit.exclusion_high:g}",
                             tolerance=f"reference < {crit.exclusion_low:g}", n=e["n"], within=e["I_violations"], pct=e["I_per_1000"],
                             wilson_lower_95="", required_pct=0, met=e["I_violations"] == 0, reportable=True, note="violations; pct column is per 1,000 estimates"))
            rt = rate_analogue(ref, est, org, h * 5, crit)
            for para, key in (("J", "J"), ("K", "K")):
                rows.append(dict(estimator=name, horizon_min=h * 5, criterion=f"analogue of (b)(1)(v)({para})",
                                 estimate_range="all", tolerance=f"indicated {'>' if para == 'J' else '<'} {'+' if para == 'J' else '−'}{crit.rate_indicated:g}, true {'<' if para == 'J' else '>'} {'−' if para == 'J' else '+'}{crit.rate_true:g} mg/dL/min",
                                 n=rt["n"], within=rt[f"{key}_count"], pct=rt[f"{key}_wrong_direction_pct"], wilson_lower_95="",
                                 required_pct=crit.rate_limit_pct, met=rt[f"{key}_met"], reportable=True, note="wrong-direction share, limit is a maximum"))
    with open(REPO_ROOT / "results" / "icgm_analogues.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"wrote results/icgm_analogues.csv ({len(rows)} rows)\n")
    print(f"{'estimator':<30}{'h':>4}" + "".join(f"{p:>14}" for p in ("(A) <70", "(B) 70-180", "(C) >180", "(D) <70/40", "(G) all 20%", "H | I", "J% | K%")))
    for h in (30, 60):
        for name in ("final model, point forecast", "linear AR, point forecast", "persistence", "c0 lower bound, 90%", "c3 lower bound, 90%", "c4 lower bound, 90%", "c5 lower bound, 90%", "c0 upper bound, 90%", "c3 upper bound, 90%"):
            g = {re.search(r"\((\w)\)$", r["criterion"]).group(1): r for r in rows if r["estimator"] == name and r["horizon_min"] == h}
            def cell(p):
                r = g[p]
                return f"{r['wilson_lower_95']:.0f}%/{r['required_pct']:.0f} n{r['n']}" if r["reportable"] else f"n={r['n']} —"
            line = f"{name:<30}{h:>4}" + "".join(f"{cell(p):>14}" for p in ("A", "B", "C", "D", "G"))
            line += f"{str(g['H']['within']) + ' | ' + str(g['I']['within']):>14}{str(g['J']['pct']) + ' | ' + str(g['K']['pct']):>14}"
            print(line)


if __name__ == "__main__":
    main()
