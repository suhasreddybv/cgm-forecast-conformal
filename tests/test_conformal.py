"""Conformal intervals and onset detection (D-045 to D-048), on synthetic data only."""
from __future__ import annotations

import numpy as np
import pytest

from src.conformal import (
    LEVELS,
    bin_of,
    c0,
    c1,
    c2,
    c3,
    choose_gamma,
    conformal_quantile,
    coverage,
    onset_detection,
    sigma_from_history,
    summarize,
)


def test_conformal_quantile_is_the_finite_sample_corrected_order_statistic():
    s = np.arange(1, 101, dtype=float)                  # 1..100
    assert conformal_quantile(s, 0.10) == 91.0           # ceil(101*0.9) = 91 -> 91st smallest
    assert conformal_quantile(s, 0.05) == 96.0
    assert conformal_quantile(s[:5], 0.10) == float("inf")   # n too small: ceil(6*0.9)=6 > 5
    assert conformal_quantile(np.array([]), 0.1) == float("inf")


@pytest.mark.parametrize("alpha", [0.2, 0.1, 0.05])
def test_c0_recovers_nominal_coverage_on_exchangeable_residuals(alpha):
    rng = np.random.default_rng(0)
    cal_pred, test_pred = rng.uniform(80, 250, 2000), rng.uniform(80, 250, 20000)
    cal_y = cal_pred + rng.normal(0, 15, 2000)
    test_y = test_pred + rng.normal(0, 15, 20000)
    iv = c0(cal_y, cal_pred, test_pred, alpha)
    cov = coverage(test_y, iv)
    assert abs(cov - (1 - alpha)) < 0.015, cov
    assert np.allclose(iv.width, iv.width[0]), "constant width"


def test_c1_sigma_uses_no_sample_at_or_after_the_origin():
    rng = np.random.default_rng(1)
    hist = rng.normal(150, 20, (50, 24))
    s = sigma_from_history(hist)
    tampered = hist.copy(); tampered[:, -1] = 10_000.0           # the origin slot
    assert np.array_equal(sigma_from_history(tampered), s), "the origin slot must not enter sigma"
    tampered2 = hist.copy(); tampered2[:, -2] = 10_000.0         # the slot before it does
    assert not np.array_equal(sigma_from_history(tampered2), s)
    assert (s >= 2.0).all(), "floored"
    flat = np.full((3, 24), 100.0)
    assert np.allclose(sigma_from_history(flat), 2.0)


def test_c1_widens_where_the_history_is_volatile_and_holds_coverage():
    rng = np.random.default_rng(2)
    n = 4000
    vol = rng.choice([3.0, 15.0], n)                            # two volatility regimes
    hist = 150 + np.cumsum(rng.normal(0, 1, (n, 24)) * vol[:, None], axis=1)
    pred = hist[:, -1]
    y = pred + rng.normal(0, 1, n) * vol * 1.5
    half = n // 2
    iv = c1(y[:half], pred[:half], hist[:half], pred[half:], hist[half:], 0.1)
    assert abs(coverage(y[half:], iv) - 0.9) < 0.03
    quiet, loud = vol[half:] == 3.0, vol[half:] == 15.0
    assert iv.width[loud].mean() > 3 * iv.width[quiet].mean()


def test_c2_bins_partition_the_forecasts_and_calibrate_per_bin():
    pred = np.array([40, 69.9, 70, 100, 180, 180.1, 300])
    assert bin_of(pred).tolist() == [0, 0, 1, 1, 1, 2, 2]
    rng = np.random.default_rng(3)
    cal_pred = np.concatenate([rng.uniform(45, 69, 300), rng.uniform(71, 179, 300), rng.uniform(181, 300, 300)])
    noise = np.where(cal_pred < 70, 5.0, np.where(cal_pred <= 180, 15.0, 30.0))
    cal_y = cal_pred + rng.normal(0, 1, 900) * noise
    test_pred = np.array([50.0, 120.0, 250.0])
    iv = c2(cal_y, cal_pred, test_pred, 0.1)
    assert iv.width[0] < iv.width[1] < iv.width[2] and iv.note == ""
    sparse = c2(cal_y[300:], cal_pred[300:], test_pred, 0.1)      # no <70 calibration scores
    assert "fallback" in sparse.note and "<70" in sparse.note


def test_c3_converges_to_nominal_on_a_stationary_series_and_tracks_a_step_change():
    """ACI holds nominal coverage on the stationary half, and restores it after the residual
    scale triples - but mostly by emitting infinite intervals (73% of steps here), because the fixed calibration
    scores cannot supply a wide enough finite quantile. Both facts are asserted: the
    recovery is real and it is not free, which is what the useless-fraction column reports."""
    rng = np.random.default_rng(4)
    cal_pred = rng.uniform(100, 200, 1000); cal_y = cal_pred + rng.normal(0, 10, 1000)
    pred = rng.uniform(100, 200, 6000)
    scale = np.where(np.arange(6000) < 3000, 10.0, 30.0)         # residual scale triples half-way
    y = pred + rng.normal(0, 1, 6000) * scale
    fixed = c0(cal_y, cal_pred, pred, 0.1)
    assert coverage(y[3000:], fixed[3000:]) < 0.5, "c0 breaks after the change"
    aci = c3(cal_y, cal_pred, y, pred, 0.1, gamma=0.02)
    first, second = aci[:3000], aci[3000:]
    assert abs(coverage(y[:3000], first) - 0.9) < 0.03
    assert abs(coverage(y[3000:], second) - 0.9) < 0.03, "ACI tracks nominal after the change"
    inf_share = float(np.isinf(second.width).mean())
    assert 0.5 < inf_share < 0.9, f"the recovery is bought mostly with infinite intervals ({inf_share:.2f})"
    assert np.isinf(first.width).mean() < 0.02
    finite = second.width[np.isfinite(second.width)]
    assert finite.mean() > 1.5 * first.width[np.isfinite(first.width)].mean()
    g = choose_gamma(cal_y, cal_pred, y[:3000], pred[:3000], 0.1)
    assert g in (0.002, 0.005, 0.01, 0.02, 0.05)


def test_summary_reports_useless_intervals():
    iv = c0(np.array([0, 0, 0, 0, 0.0]), np.zeros(5), np.zeros(3), 0.5)
    s = summarize(np.zeros(3), iv)
    assert s["n"] == 3 and s["coverage"] == 1.0
    wide = c0(np.array([0, 300.0]), np.zeros(2), np.zeros(2), 0.1)       # n too small -> inf
    assert summarize(np.zeros(2), wide)["useless_fraction"] == 1.0


def test_onset_counting_on_a_hand_built_series():
    #              t:  0    1    2    3    4    5    6    7
    y      = np.array([90,  65,  60,  80,  68,  55,  75,  40.0])
    origin = np.array([95,  90,  65,  60,  80,  68,  55,  75.0])   # reading at the forecast origin
    lower  = np.array([65,  60,  50,  70,  75,  50,  60,  30.0])   # interval lower bound
    r = onset_detection(y, origin, lower, 70.0)
    # at risk (origin >= 70): t = 0, 1, 4, 7. Onsets among them (y < 70): t = 1, 4, 7.
    # Already below at the origin and still below: t = 2, 5 (excluded from everything).
    assert r["n_at_risk"] == 4 and r["n_onsets"] == 3 and r["n_already_below_at_origin"] == 2
    # alarms (lower < 70) among at-risk: t=0 (65: false alarm), t=1 (60: hit), t=4 (75: miss), t=7 (30: hit)
    assert (r["tp"], r["fn"], r["fp"], r["tn"]) == (2, 1, 1, 0)
    assert r["sensitivity"] == round(2 / 3, 4) and r["specificity"] == 0.0 and r["ppv"] == round(2 / 3, 4)
    assert r["censored"] is True and r["n_onsets_at_sensor_floor"] == 1          # t=7 reads 40
    assert r["false_alarms_per_day"] == round(1 / (4 * 5 / 1440), 3)


def test_persistence_detects_no_onset_by_construction():
    rng = np.random.default_rng(5)
    origin = rng.uniform(60, 200, 500); y = origin + rng.normal(0, 20, 500)
    r = onset_detection(y, origin, origin, 70.0)                 # persistence: lower = point = origin
    assert r["tp"] == 0 and r["fp"] == 0 and (r["sensitivity"] == 0.0 or r["n_onsets"] == 0)
    assert len(LEVELS) == 3


def test_operating_curve_reads_sensitivity_at_a_false_alarm_rate():
    from src.eval.operating_curve import DELTAS, LEVELS_FINE, pooled_onset, sensitivity_at
    rows = [dict(threshold=70.0, horizon_min=30, curve="x", parameter=i, false_alarms_per_day=fa, sensitivity=s)
            for i, (fa, s) in enumerate([(0.5, 0.1), (1.5, 0.3), (4.0, 0.8)])]
    assert sensitivity_at(rows, 70.0, 30, "x", 1.0) == pytest.approx(0.2)        # interpolated between 0.5 and 1.5
    assert np.isnan(sensitivity_at(rows, 70.0, 30, "x", 10.0))                  # never reached
    assert DELTAS[0] == 0 and DELTAS[-1] == 60 and len(DELTAS) == 31
    assert LEVELS_FINE[0] == 0.05 and LEVELS_FINE[-1] == 0.95 and len(LEVELS_FINE) == 19
    parts = [dict(tp=2, fn=1, fp=1, tn=10, n_at_risk=14, n_onsets_at_sensor_floor=0),
             dict(tp=0, fn=1, fp=3, tn=20, n_at_risk=24, n_onsets_at_sensor_floor=1)]
    p = pooled_onset(parts)
    assert p["n_onsets"] == 4 and p["sensitivity"] == 0.5 and p["censored"] is True
    assert p["false_alarms_per_day"] == pytest.approx(4 / (38 * 5 / 1440))


def test_a_constant_offset_lower_bound_is_the_point_forecast_sweep():
    """c0's lower bound at any level is the point forecast minus a constant: the control."""
    rng = np.random.default_rng(6)
    cal_pred = rng.uniform(80, 200, 500); cal_y = cal_pred + rng.normal(0, 12, 500)
    pred = rng.uniform(60, 200, 300)
    iv = c0(cal_y, cal_pred, pred, 0.2)
    offset = pred - iv.lower
    assert np.allclose(offset, offset[0]), "constant width per calibration set"
    origin = rng.uniform(70, 200, 300); y = pred + rng.normal(0, 12, 300)
    a = onset_detection(y, origin, iv.lower, 70.0)
    b = onset_detection(y, origin, pred - offset[0], 70.0)
    assert (a["tp"], a["fp"]) == (b["tp"], b["fp"])
