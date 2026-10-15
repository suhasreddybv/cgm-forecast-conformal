"""c4 against c0 and c3 in the hypoglycaemic bins (D-062): one-sided coverage of the lower
bound, the median distance from the forecast to the bound, and the misses below it, for
the true-glucose bins < 54, 54-70 and 70-180, both horizons, both variants, 2020 cohort.
Writes results/c4_hypoglycaemic_bins.csv. Aggregate only.
"""
from __future__ import annotations

import csv

import numpy as np

from src.conformal import c0, c3, c4, choose_gamma
from src.data.loader import COHORT_2020, REPO_ROOT
from src.data.windows import HORIZONS
from src.eval.coverage_breakdown import GLUCOSE_NAMES, MIN_N, glucose_bin, patient_sets

LEVELS = (0.80, 0.90, 0.95)


def main() -> None:
    rows = []
    for variant in ("bglp", "primary"):
        sets = patient_sets(variant)
        for h in HORIZONS:
            for level in LEVELS:
                alpha = 1 - level
                per = {"c0": [], "c3": [], "c4": []}
                ys, ps = [], []
                for pid in COHORT_2020:
                    s = sets[(pid, h)]
                    per["c0"].append(c0(s["cal_y"], s["cal_p"], s["p"], alpha).lower)
                    per["c3"].append(c3(s["cal_y"], s["cal_p"], s["y"], s["p"], alpha,
                                        choose_gamma(s["fit_y"], s["fit_p"], s["cal_y"], s["cal_p"], alpha)).lower)
                    per["c4"].append(c4(s["cal_y"], s["cal_p"], s["p"], alpha).lower)
                    ys.append(s["y"]); ps.append(s["p"])
                y, p = np.concatenate(ys), np.concatenate(ps)
                b = glucose_bin(y)
                for m, parts in per.items():
                    lower = np.concatenate(parts)
                    for bi, name in enumerate(GLUCOSE_NAMES[:3]):
                        msk = b == bi
                        n = int(msk.sum())
                        rows.append(dict(method=m, level=level, variant=variant, horizon_min=h * 5, true_bin=name, n=n,
                                         one_sided_coverage=round(float(np.mean(y[msk] >= lower[msk])), 4) if n else "",
                                         miss_below=int((y[msk] < lower[msk]).sum()),
                                         median_distance_forecast_to_bound=round(float(np.median((p - lower)[msk])), 2) if n else "",
                                         reportable=n >= MIN_N))
    with open(REPO_ROOT / "results" / "c4_hypoglycaemic_bins.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"wrote results/c4_hypoglycaemic_bins.csv ({len(rows)} rows)\n")
    print(f"{'variant':<8}{'h':>4}{'bin':<8}" + "".join(f"{m + ' cov/below/dist':>22}" for m in ("c0", "c3", "c4")) + "   (90% level)")
    for variant in ("bglp", "primary"):
        for h in HORIZONS:
            for name in GLUCOSE_NAMES[:3]:
                line = f"{variant:<8}{h * 5:>4}{name:<8}"
                for m in ("c0", "c3", "c4"):
                    r = next(r for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["true_bin"]) == (m, 0.9, variant, h * 5, name))
                    cell = f"{r['one_sided_coverage']:.3f}/{r['miss_below']}/{r['median_distance_forecast_to_bound']:.0f}" if r["n"] else "—"
                    line += f"{cell + (f' n{r["n"]}' if m == 'c0' else ''):>22}"
                print(line)


if __name__ == "__main__":
    main()
