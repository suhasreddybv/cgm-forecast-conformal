"""The coverage-breakdown bins (D-051): edge values land in the right bin, and every window
lands in exactly one bin per dimension."""
from __future__ import annotations

import numpy as np

from src.eval.coverage_breakdown import (
    GAP_NAMES,
    GLUCOSE_NAMES,
    ROC_NAMES,
    gap_bin,
    glucose_bin,
    minutes_since_long_gap,
    roc_bin,
)


def test_glucose_bins_lower_edge_inclusive():
    y = [40, 53.9, 54, 69.9, 70, 179.9, 180, 249.9, 250, 400]
    names = [GLUCOSE_NAMES[i] for i in glucose_bin(y)]
    assert names == ["<54", "<54", "54-70", "54-70", "70-180", "70-180", "180-250", "180-250", ">250", ">250"]


def test_rate_of_change_bins_in_mg_dl_per_minute():
    # 5-minute step: a 5 mg/dL difference is 1 mg/dL/min
    origin = np.array([100, 105, 110, 115, 116, 100])
    prev = np.array([100, 100, 100, 100, 100, 120])          # 0, 1, 2, 3, 3.2, 4 mg/dL/min
    names = [ROC_NAMES[i] for i in roc_bin(origin, prev)]
    assert names == ["<=1", "<=1", "1-2", "2-3", ">3", ">3"]


def test_gap_bins_and_no_prior_gap():
    m = [0, 29.9, 30, 119.9, 120, 1000, np.nan]
    names = [GAP_NAMES[i] for i in gap_bin(m)]
    assert names == ["<30 min", "<30 min", "30-120 min", "30-120 min", ">120 min (no gap in 2-h window)",
                     ">120 min (no gap in 2-h window)", "no prior gap"]


def test_minutes_since_long_gap_from_a_series():
    t0 = np.datetime64("2021-01-01T00:00:00", "s")
    # readings every 5 min for an hour, then a 2-hour hole, then every 5 min again
    ts = np.concatenate([t0 + np.arange(0, 13) * np.timedelta64(300, "s"),
                         t0 + np.timedelta64(3 * 3600, "s") + np.arange(0, 13) * np.timedelta64(300, "s")])
    origins = [t0 + np.timedelta64(30 * 60, "s"),                 # before the gap: no prior gap
               t0 + np.timedelta64(3 * 3600, "s"),                # exactly at the gap end: 0 min
               t0 + np.timedelta64(3 * 3600 + 45 * 60, "s")]      # 45 min after
    m = minutes_since_long_gap(ts, np.array(origins, dtype="datetime64[s]"))
    assert np.isnan(m[0]) and m[1] == 0.0 and m[2] == 45.0
    assert [GAP_NAMES[i] for i in gap_bin(m)] == ["no prior gap", "<30 min", "30-120 min"]


def test_every_window_lands_in_exactly_one_bin_per_dimension():
    rng = np.random.default_rng(0)
    y = rng.uniform(40, 400, 5000); origin = rng.uniform(40, 400, 5000); prev = origin + rng.normal(0, 10, 5000)
    since = np.where(rng.random(5000) < 0.1, np.nan, rng.uniform(0, 600, 5000))
    for fn, args, k in ((glucose_bin, (y,), len(GLUCOSE_NAMES)), (roc_bin, (origin, prev), len(ROC_NAMES)), (gap_bin, (since,), len(GAP_NAMES))):
        b = fn(*args)
        assert b.shape == (5000,) and b.min() >= 0 and b.max() < k
        assert np.bincount(b, minlength=k).sum() == 5000


def test_breakdown_rows_partition_the_cohort_when_the_results_exist():
    import csv
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "results" / "coverage_breakdown.csv"
    if not path.exists():
        return
    rows = list(csv.DictReader(open(path)))
    for dim in ("true_glucose", "rate_of_change", "time_since_gap"):
        g = [r for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["cohort"], r["dimension"]) == ("c0", "0.9", "bglp", "30", "2020", dim)]
        total = sum(int(r["n"]) for r in g)
        pat = sum(int(r["n"]) for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["dimension"]) == ("c0", "0.9", "bglp", "30", "patient") and r["cohort"] == "2020")
        assert total == pat, (dim, total, pat)


def test_miss_split_sums_to_n_and_reportability_follows_the_rule():
    import csv
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "results" / "coverage_breakdown.csv"
    if not path.exists():
        return
    rows = list(csv.DictReader(open(path)))
    assert "miss_below" in rows[0] and "miss_above" in rows[0] and "reportable" in rows[0]
    for r in rows:
        n = int(r["n"])
        assert int(r["covered"]) + int(r["miss_below"]) + int(r["miss_above"]) == n, r
        assert (r["reportable"] == "True") == (n >= 100), r
        if n:
            assert abs(float(r["coverage"]) - int(r["covered"]) / n) < 1e-3


def test_c0_width_is_constant_within_patient():
    import numpy as np
    from src.conformal import c0
    rng = np.random.default_rng(1)
    cal_p = rng.uniform(80, 250, 400); cal_y = cal_p + rng.normal(0, 10, 400)
    iv = c0(cal_y, cal_p, rng.uniform(60, 300, 1000), 0.1)
    assert np.ptp(iv.width) < 1e-6        # constant up to floating-point: (p+q)-(p-q) can differ by an ulp
