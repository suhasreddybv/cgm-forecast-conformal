"""The three baselines (D-018 to D-021): definitions pinned on synthetic windows, the
fitting rule, NaN-covariate handling, and the sanity checks on the real results."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from src.data import loader
from src.data.loader import PATIENTS
from src.data.windows import WindowSet
from src.eval.baselines import (
    covariate_ok,
    design,
    fit_ols,
    linear_extrapolation,
    persistence,
    sanity_checks,
)

REPO = Path(__file__).resolve().parents[1]
have_data = (loader.DATA_ROOT / "2018" / "train" / "559-ws-training.xml").exists()
needs_data = pytest.mark.skipif(not have_data, reason="OhioT1DM not present; see data/README.md")


def make_ws(history, target, horizon=6, basal=None, bolus=None, carbs=None, pid="000", split="test"):
    history = np.asarray(history, float)
    n, H = history.shape
    z = np.zeros((n, H))
    return WindowSet(pid, split, horizon, H,
                     target_ts=np.arange(n).astype("datetime64[s]"), target=np.asarray(target, float),
                     history=history, max_bracket_s=np.zeros(n, int),
                     basal=z if basal is None else np.asarray(basal, float),
                     bolus=z if bolus is None else np.asarray(bolus, float),
                     carbs=z if carbs is None else np.asarray(carbs, float),
                     n_readings=n, n_candidates=n, n_dropped_long_gap=0, n_dropped_no_history=0)


def test_persistence_is_the_last_history_slot():
    ws = make_ws([[100, 110, 120], [90, 80, 70]], [0, 0])
    assert persistence(ws).tolist() == [120.0, 70.0]


@pytest.mark.parametrize("k", [2, 3, 6])
def test_linear_extrapolation_follows_an_exact_line(k):
    """On a perfectly linear history the extrapolation is exact for every k and horizon,
    and the slope is per step, so a 60-minute horizon extends twice as far as 30."""
    hist = 100 + 2.0 * np.arange(6)[None, :]           # +2 mg/dL per step, last slot = 110
    for h in (6, 12):
        ws = make_ws(hist, [0], horizon=h)
        assert np.allclose(linear_extrapolation(ws, k), 110 + 2.0 * h)


def test_two_point_extrapolation_is_the_last_difference_times_the_horizon():
    ws = make_ws([[0, 0, 0, 0, 100, 103]], [0], horizon=6)
    assert np.isclose(linear_extrapolation(ws, 2)[0], 103 + 3 * 6)


def test_ols_recovers_a_known_linear_generator_from_training_only():
    rng = np.random.default_rng(0)
    H = 6
    hist = rng.normal(150, 40, (500, H))
    beta = np.array([5.0, 0.1, -0.2, 0.3, 0.0, 0.4, 0.5])
    target = design(make_ws(hist, np.zeros(500)), False) @ beta
    train = make_ws(hist, target, split="train")
    b = fit_ols([train], covariates=False)
    assert np.allclose(b, beta, atol=1e-8)
    # with covariates the design widens to 1 + 4H columns
    assert design(train, True).shape == (500, 1 + 4 * H)


def test_nan_covariates_drop_only_the_with_covariate_rows():
    hist = np.tile(np.arange(6.0), (4, 1))
    basal = np.ones((4, 6))
    basal[1, 2] = np.nan
    ws = make_ws(hist, np.ones(4), basal=basal)
    assert covariate_ok(ws).tolist() == [True, False, True, True]
    X = design(ws, True)
    assert np.isnan(X[1]).any() and not np.isnan(X[covariate_ok(ws)]).any()
    # fit_ols with covariates ignores the NaN window; without, it keeps all four
    b_cov = fit_ols([ws], True)
    assert np.isfinite(b_cov).all()
    assert len(design(ws, False)) == 4


# ------------------------------------------------------------------- real data

def _rows(name):
    path = REPO / "results" / name
    return list(csv.DictReader(open(path)))


@needs_data
def test_results_tables_are_consistent_with_each_other_and_carry_no_patient_data():
    pp, agg = _rows("baselines_per_patient.csv"), _rows("baselines.csv")
    assert len(pp) == 12 * 2 * (3 + 3 + 3 * 4)          # p0 per H, l1 per k, l2 per H x fit x cov
    assert len(agg) == len(pp) // 12 * 3
    for col in pp[0]:
        assert not set(col.split("_")) & {"ts", "timestamp", "time", "value", "values", "glucose"}, col
    # the published-comparison columns exist and are blank: Suhas fills them from his reading table
    for r in agg:
        assert r["published_rmse"] == "" and r["published_source"] == "" and r["like_for_like"] == ""
    # the aggregate's n is the sum of its patients' n, and pooled RMSE is n-weighted
    for r in agg:
        g = [p for p in pp if (p["method"], p["variant"], p["fit"], p["covariates"], p["horizon_steps"],
                               p["window_set_H"]) == (r["method"], r["variant"], r["fit"], r["covariates"],
                                                       r["horizon_min"] and str(int(r["horizon_min"]) // 5), r["window_set_H"])
             and (r["cohort"] == "pooled" or p["cohort"] == r["cohort"])]
        assert len(g) == int(r["n_patients"]) == (12 if r["cohort"] == "pooled" else 6)
        assert sum(int(p["n_scored"]) for p in g) == int(r["n_scored"])


@needs_data
def test_scored_n_equals_the_effective_n_table():
    pp = _rows("baselines_per_patient.csv")
    eff = {(r["patient"], r["horizon_steps"], r["history_steps"]): int(r["evaluable_targets"])
           for r in _rows("effective_n.csv") if r["split"] == "test"}
    for r in pp:
        if r["covariates"] == "False":
            assert int(r["n_scored"]) == eff[(r["patient"], r["horizon_steps"], r["window_set_H"])], r
        else:
            assert int(r["n_scored"]) + int(r["n_dropped_nan_covariate"]) == eff[(r["patient"], r["horizon_steps"], r["window_set_H"])]


@needs_data
def test_the_three_sanity_checks_pass_on_the_committed_results():
    pp = _rows("baselines_per_patient.csv")
    for r in pp:
        for k in ("rmse", "mae", "mape"):
            r[k] = float(r[k])
        r["horizon_steps"], r["window_set_H"] = int(r["horizon_steps"]), int(r["window_set_H"])
        r["k"] = int(r["k"]) if r["k"] else ""
        r["covariates"] = r["covariates"] == "True"
    for name, ok, detail in sanity_checks(pp):
        assert ok, f"{name}: {detail}"


@needs_data
def test_nan_covariate_windows_exist_only_in_575s_training_file():
    pp = _rows("baselines_per_patient.csv")
    test_drops = {r["patient"] for r in pp if r["covariates"] == "True" and int(r["n_dropped_nan_covariate"])}
    train_drops = {r["patient"] for r in pp if r["covariates"] == "True" and r["fit"] == "per-patient"
                   and int(r["n_train_dropped_nan"])}
    assert test_drops == set() and train_drops == {"575"}
