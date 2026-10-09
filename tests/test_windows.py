"""The windower's contract with the frozen protocol (docs/protocol.md, D-005, D-013 to D-016).

Synthetic cases pin exact window counts and covariate alignment without the dataset.
The real-data cases assert, over every patient x split x horizon x H: no target is
interpolated, no history crosses more than 30 minutes, history ends strictly before the
target, the test split never reads past the target or outside the patient, the warm-up
is handled as the protocol says, and the counts reconcile with the gap analysis.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from src.data import loader
from src.data.loader import PATIENTS, load_patient
from src.data.windows import (
    DROPPED_LONG_GAP,
    DROPPED_NO_HISTORY,
    GAP_S,
    HISTORY_LENGTHS,
    HORIZONS,
    KEPT,
    MAX_INTERP_GAP_S,
    STEP_S,
    Channel,
    _series_for,
    align_covariates,
    basal_rate_at,
    bolus_delivered_at,
    build_windows,
    windows_for,
)

REPO = Path(__file__).resolve().parents[1]
have_data = (loader.DATA_ROOT / "2018" / "train" / "559-ws-training.xml").exists()
needs_data = pytest.mark.skipif(not have_data, reason="OhioT1DM not present; see data/README.md")

T0 = np.datetime64("2021-01-01T00:00:00", "s")


def series(slots):
    """A CGM series at exact cadence on the given slot indices, values = slot index."""
    slots = np.asarray(slots)
    return T0 + (slots * STEP_S).astype("timedelta64[s]"), slots.astype(float)


def one_gap(n0, g, n1):
    """n0 readings, then g missing slots, then n1 readings."""
    return series(np.r_[np.arange(n0), n0 + g + np.arange(n1)])


def expected_counts(n0, g, n1, h, H):
    """Closed form for a single isolated gap of g missing slots (derivation in D-016).

    Start of series: the first h+H-1 readings have no history. The gap touches the
    history of H-1+min(g,h) targets after it: H-1 whose history span straddles it, plus
    min(g,h) whose history END falls inside it.
    """
    no_hist = h + H - 1
    touched = H - 1 + min(g, h)
    long = touched if g * STEP_S + STEP_S > MAX_INTERP_GAP_S else 0   # interval = (g+1)*300 s
    interp = touched - long
    return dict(n=n0 + n1 - no_hist - long, no_hist=no_hist, long=long, interp=interp)


# --------------------------------------------------------------- synthetic: counts

@pytest.mark.parametrize("h,H", [(6, 6), (6, 12), (12, 6), (12, 24), (6, 24)])
def test_no_gap_keeps_everything_with_history(h, H):
    ts, v = series(np.arange(200))
    b = build_windows(ts, v, h, H)
    assert (b.status == DROPPED_NO_HISTORY).sum() == h + H - 1
    assert (b.status == DROPPED_LONG_GAP).sum() == 0
    assert len(b.target_idx) == 200 - (h + H - 1)
    assert (b.max_bracket_s[b.kept] == 0).all(), "exact cadence: nothing interpolated"


@pytest.mark.parametrize("g", [1, 3, 5])             # intervals of 10, 20, 30 min: interpolated
@pytest.mark.parametrize("h,H", [(6, 6), (12, 6), (6, 24), (12, 24)])
def test_short_gap_is_interpolated_not_dropped(g, h, H):
    ts, v = one_gap(100, g, 100)
    b = build_windows(ts, v, h, H)
    e = expected_counts(100, g, 100, h, H)
    assert (b.status == DROPPED_LONG_GAP).sum() == 0
    assert len(b.target_idx) == e["n"]
    assert int((b.max_bracket_s[b.kept] > GAP_S).sum()) == e["interp"]
    # the interpolated values are linear in time: value == slot index everywhere
    slots_idx = (b.slots_s - T0.astype("datetime64[s]").astype(np.int64)) / STEP_S
    assert np.allclose(b.history, slots_idx)


@pytest.mark.parametrize("g", [6, 7, 20, 100])      # intervals of 35 min and up: dropped
@pytest.mark.parametrize("h,H", [(6, 6), (12, 6), (6, 24), (12, 24)])
def test_long_gap_drops_exactly_the_windows_that_touch_it(g, h, H):
    ts, v = one_gap(100, g, 100)
    b = build_windows(ts, v, h, H)
    e = expected_counts(100, g, 100, h, H)
    assert (b.status == DROPPED_LONG_GAP).sum() == e["long"]
    assert (b.status == DROPPED_NO_HISTORY).sum() == e["no_hist"]
    assert len(b.target_idx) == e["n"]
    assert (b.max_bracket_s[b.kept] <= MAX_INTERP_GAP_S).all()


def test_the_30_minute_boundary_is_inclusive():
    """A 30-minute interval (5 missing slots) is interpolated; 35 minutes is dropped."""
    for g, dropped in ((5, 0), (6, None)):
        ts, v = one_gap(50, g, 50)
        b = build_windows(ts, v, 6, 6)
        n_long = int((b.status == DROPPED_LONG_GAP).sum())
        assert (n_long == 0) if dropped == 0 else (n_long > 0)


def test_cadence_jitter_is_not_reported_as_a_gap():
    """A 301 s interval shifts earlier readings by a second against the target grid.
    They are interpolated across 301 s, which is not a gap, and the value error is tiny."""
    ts, v = series(np.arange(100))
    ts = ts.copy()
    ts[50:] += np.timedelta64(1, "s")
    b = build_windows(ts, v, 6, 6)
    assert (b.status == DROPPED_LONG_GAP).sum() == 0
    assert b.max_bracket_s[b.kept].max() == 301
    assert not (b.max_bracket_s[b.kept] > GAP_S).any()
    assert np.abs(b.history - np.round(b.history)).max() < 1 / 300


def test_non_multiple_gap_is_a_gap_and_is_interpolated():
    ts, v = series(np.arange(100))
    ts = ts.copy()
    ts[50:] += np.timedelta64(783 - 300, "s")      # one interval of 783 s (13 min)
    b = build_windows(ts, v, 6, 6)
    assert (b.status == DROPPED_LONG_GAP).sum() == 0
    assert b.max_bracket_s[b.kept].max() == 783
    assert (b.max_bracket_s[b.kept] > GAP_S).sum() > 0


def test_targets_are_real_readings_and_history_ends_before_them():
    ts, v = one_gap(80, 3, 80)
    for h in HORIZONS:
        b = build_windows(ts, v, h, 12)
        assert np.isin(b.target_idx, np.arange(len(ts))).all()
        tsec = ts.astype(np.int64)
        assert (b.slots_s[:, -1] == tsec[b.target_idx] - h * STEP_S).all()
        assert (b.slots_s[:, -1] < tsec[b.target_idx]).all()
        assert (np.diff(b.slots_s, axis=1) == STEP_S).all()


def test_eval_from_uses_earlier_readings_as_history_only():
    ts, v = series(np.arange(100))
    b = build_windows(ts, v, 6, 12, eval_from=12)
    assert b.target_idx.min() == 17          # 12 + (6+12-1) = 29 would hold without the prior readings
    assert len(b.candidates) == 88
    assert (b.status == DROPPED_NO_HISTORY).sum() == 5


def test_rejects_unordered_timestamps_and_bad_sizes():
    ts, v = series([0, 2, 1])
    with pytest.raises(ValueError, match="strictly increasing"):
        build_windows(ts, v, 6, 6)
    ts, v = series(np.arange(10))
    with pytest.raises(ValueError):
        build_windows(ts, v, 0, 6)


# ----------------------------------------------------------- synthetic: covariates

def _chan(name, ts_list, values, rows):
    ts = np.array([np.datetime64(t, "s") for t in ts_list], dtype="datetime64[s]")
    return Channel(name, ts, np.array(values, float), {"rows": rows})


def _rows_interval(items):
    return [{"ts_begin": b.strftime("%d-%m-%Y %H:%M:%S"), "ts_end": e.strftime("%d-%m-%Y %H:%M:%S"),
             **extra} for b, e, extra in items]


def test_basal_is_a_step_function_with_temp_basal_override():
    import datetime as dt
    d = dt.datetime(2021, 1, 1)
    basal = _chan("basal", ["2021-01-01T00:00:00", "2021-01-01T04:00:00"], [0.8, 1.2], [])
    temp = Channel("temp_basal", np.array([], "datetime64[s]"), np.array([]), {"rows": _rows_interval(
        [(d + dt.timedelta(hours=1), d + dt.timedelta(hours=2), {"value": "0.0"})])})
    at = np.array([dt.datetime(2021, 1, 1, hh, mm) for hh, mm in
                   ((0, 30), (1, 0), (1, 30), (2, 0), (3, 0), (5, 0))], dtype="datetime64[s]").astype(np.int64)
    rate = basal_rate_at(at, basal, temp)
    assert rate.tolist() == [0.8, 0.0, 0.0, 0.8, 0.8, 1.2]        # begin inclusive, end exclusive
    before = basal_rate_at(np.array([at[0] - 7200]), basal, temp)
    assert np.isnan(before).all(), "no rate is known before the first basal event"


def test_bolus_lands_in_its_bin_and_square_bolus_spreads():
    import datetime as dt
    d = dt.datetime(2021, 1, 1, 0, 0, 0)
    rows = _rows_interval([(d + dt.timedelta(minutes=7), d + dt.timedelta(minutes=7), {"dose": "2.0", "type": "normal"}),
                           (d + dt.timedelta(minutes=20), d + dt.timedelta(minutes=50), {"dose": "3.0", "type": "square"})])
    bolus = Channel("bolus", np.array([], "datetime64[s]"), np.array([]), {"rows": rows})
    slots = (np.datetime64(d, "s").astype(np.int64) + np.arange(1, 13) * STEP_S)[None, :]   # bins ending 5..60 min
    delivered = bolus_delivered_at(slots, bolus) - bolus_delivered_at(slots - STEP_S, bolus)
    assert delivered[0, 0] == 0.0 and delivered[0, 1] == 2.0           # minute 7 is in (5, 10]
    assert np.isclose(delivered[0, 3:10].sum(), 3.0)                   # square bolus 20..50 min
    assert np.isclose(delivered[0, 4], 0.5)                            # 5 of 30 minutes per full bin
    assert np.isclose(delivered.sum(), 5.0)


def test_carbs_bin_boundaries_are_left_open_right_closed():
    meal = _chan("meal", ["2021-01-01T00:05:00", "2021-01-01T00:05:01"], [30, 10], [])
    basal = _chan("basal", ["2020-12-31T00:00:00"], [1.0], [])
    empty = Channel("x", np.array([], "datetime64[s]"), np.array([]), {"rows": []})
    slots = (np.datetime64("2021-01-01T00:00:00", "s").astype(np.int64) + np.arange(1, 4) * STEP_S)[None, :]
    _, _, carbs = align_covariates(slots, basal, empty, empty, meal)
    assert carbs.tolist() == [[30.0, 10.0, 0.0]]


# -------------------------------------------------------------------- real data

ALL = [(p, s, h, H) for p in PATIENTS for s in ("train", "test") for h in HORIZONS for H in HISTORY_LENGTHS]


@needs_data
@pytest.mark.parametrize("pid,split,h,H", ALL)
def test_contract_holds_on_every_patient_split_horizon_and_history(pid, split, h, H):
    ws = windows_for(pid, split, h, H)
    rec = load_patient(pid, split)
    raw_ts = rec.cgm.ts
    # 1. no target is interpolated: every target timestamp exists in the raw CGM stream
    assert np.isin(ws.target_ts, raw_ts).all()
    # 2. no history span crosses more than 30 minutes
    assert (ws.max_bracket_s <= MAX_INTERP_GAP_S).all()
    # 3. history end strictly before the target, by exactly the horizon
    assert (ws.history_end_ts < ws.target_ts).all()
    assert (ws.target_ts - ws.history_end_ts == np.timedelta64(h * STEP_S, "s")).all()
    # 4. no history sample later than its target; nothing outside this patient's own data
    hts = ws.history_ts()
    assert (hts < ws.target_ts[:, None]).all()
    first = load_patient(pid, "train").cgm.ts[0]
    assert (hts >= first).all() and (ws.target_ts <= raw_ts[-1]).all()
    # 5. shapes and the covariate arrays are aligned to the same windows
    n = len(ws)
    assert ws.history.shape == ws.basal.shape == ws.bolus.shape == ws.carbs.shape == (n, H)
    assert np.isfinite(ws.history).all() and np.isfinite(ws.target).all()
    assert (ws.bolus >= 0).all() and (ws.carbs >= 0).all()
    assert n + ws.n_dropped_long_gap + ws.n_dropped_no_history == ws.n_candidates


@needs_data
@pytest.mark.parametrize("pid", PATIENTS)
def test_test_split_warm_up_as_the_protocol_specifies(pid):
    """2020: the first 12 test readings serve as history and are never targets; 2018: the
    first test reading is the first target. Early targets draw history from the training
    tail, which ends exactly 300 s before the test file starts."""
    rec = load_patient(pid, "test")
    train = load_patient(pid, "train")
    assert rec.cgm.ts[0] - train.cgm.ts[-1] == np.timedelta64(STEP_S, "s")
    cgm, _, eval_from = _series_for(rec)
    assert eval_from == len(train.cgm) + rec.eval_start_index
    # the scored set starts where the published rule says, whatever history is available
    # earlier: 2020 after the first hour, 2018 at the first test reading (D-017)
    for h in HORIZONS:
        for H in HISTORY_LENGTHS:
            assert windows_for(pid, "test", h, H).target_ts.min() >= rec.cgm.ts[rec.eval_start_index]
    for h in HORIZONS:
        ws = windows_for(pid, "test", h, 24)
        assert ws.n_candidates == len(rec.cgm) - rec.eval_start_index
        assert ws.n_dropped_no_history == 0, "the training tail always supplies history"
        # the first candidate is the first reading after the warm-up. It is the first target
        # unless its history crosses a gap over 30 minutes in the training tail - which is
        # what happens to 552 - and then every candidate before the first target was dropped
        # for exactly that reason, never for lack of readings.
        b = build_windows(cgm.ts, cgm.values, h, 24, eval_from)
        assert b.candidates[0] == eval_from
        before_first = b.status[: np.flatnonzero(b.kept)[0]]
        assert (before_first == DROPPED_LONG_GAP).all()
        assert (pid == "552") == (len(before_first) > 0), pid
        assert ws.target_ts[0] == cgm.ts[b.target_idx[0]]
        # early targets draw history from the training tail, never from before it
        assert ws.history_ts()[0].min() < rec.cgm.ts[0], "first test target's history reaches into training"
        assert ws.history_ts()[0].min() >= train.cgm.ts[0]


@needs_data
def test_history_values_equal_the_readings_where_slots_are_exact():
    rec = load_patient("588", "train")
    ws = windows_for("588", "train", 6, 12)
    t = rec.cgm.ts.astype(np.int64)
    idx = np.searchsorted(t, ws.history_ts().astype(np.int64))
    exact = t[np.minimum(idx, len(t) - 1)] == ws.history_ts().astype(np.int64)
    assert exact.mean() > 0.95
    assert np.array_equal(ws.history[exact], rec.cgm.values[np.minimum(idx, len(t) - 1)][exact])


@needs_data
def test_reconciles_with_the_gap_analysis():
    """The index-space crossing count of gap_analysis.py must equal the number of targets
    the windower finds touched on the history side, up to clustered gaps: when two gaps
    lie within h+H slots of each other the two counts attribute a target to different
    gaps (one sees its history end inside the first gap, the other sees its index-space
    window cross the second). Measured: the windower's history-side count is never above
    the naive one, the per-file shortfall is at most 27 targets, and pooled over all files
    it is below 1.5% of the naive count (D-016)."""
    from src.data.effective_n import naive_crossing
    ga = {(r["patient"], r["split"]): r for r in csv.DictReader(open(REPO / "results" / "gap_analysis.csv"))}
    pooled_naive = pooled_resid = 0
    for pid in PATIENTS:
        for split in ("train", "test"):
            cgm, _, eval_from = _series_for(load_patient(pid, split))
            for h in HORIZONS:
                for H in HISTORY_LENGTHS:
                    naive, n_idx = naive_crossing(cgm.ts, H, h, eval_from)
                    b = build_windows(cgm.ts, cgm.values, h, H, eval_from)
                    has = b.status != DROPPED_NO_HISTORY
                    touched = has & (b.max_bracket_s > GAP_S)
                    end_in = has & (b.end_bracket_s > GAP_S)
                    hist_side = int((touched & ~end_in).sum())
                    assert naive - 30 <= hist_side <= naive, (pid, split, h, H, hist_side, naive)
                    pooled_naive += naive
                    pooled_resid += naive - hist_side
                    # the published naive fraction was computed over the file's own readings
                    if split == "train":
                        assert abs(naive / n_idx - float(ga[(pid, split)][f"windows_crossing_H{H}"])) < 0.003
    assert pooled_resid / pooled_naive < 0.015


@needs_data
def test_effective_n_table_matches_the_windower_and_carries_no_patient_data():
    path = REPO / "results" / "effective_n.csv"
    rows = list(csv.DictReader(open(path)))
    assert len(rows) == len(PATIENTS) * 2 * len(HORIZONS) * len(HISTORY_LENGTHS)
    for col in rows[0]:
        tokens = set(col.split("_"))
        assert not tokens & {"ts", "timestamp", "time", "value", "values", "glucose", "mgdl"}, col
    for r in rows:
        for col, cell in r.items():
            if col not in ("patient", "cohort", "split"):
                float(cell)                       # every cell is a count or a fraction
    for r in rows:
        if r["patient"] in ("552", "588") and r["split"] == "test":
            ws = windows_for(r["patient"], "test", int(r["horizon_steps"]), int(r["history_steps"]))
            assert int(r["evaluable_targets"]) == len(ws)
            assert int(r["readings"]) == ws.n_readings
    # 552's 40% missingness is paid in absolute n, not in the fraction of its own readings:
    # readings inside its 118-hour gap never existed, so they are not targets to lose.
    worst = min(rows, key=lambda r: float(r["fraction_of_expected"]))
    assert (worst["patient"], worst["split"]) == ("552", "test")
    assert float(worst["fraction_of_expected"]) < 0.55
    assert min(float(r["fraction_of_readings"]) for r in rows if r["patient"] == "552") > 0.85


# ------------------------------------------------------ the conforming (bglp) variant, D-027

from src.data.windows import build_windows_zoh  # noqa: E402


def test_zoh_values_every_slot_from_the_last_reading_at_or_before_it():
    """No bracketing reading after a slot is ever used: brute force against the vectorised path."""
    ts, v = one_gap(60, 7, 60)                       # a 40-minute gap: the primary rule would drop
    ts = ts.copy(); ts[30:] += np.timedelta64(2, "s")   # and a jitter, so slots fall between readings
    slots, hist, tidx, maxb, hold, n_nohist = build_windows_zoh(ts, v, 6, 6)
    t = ts.astype(np.int64)
    for i in range(len(tidx)):
        for k in range(6):
            last = np.flatnonzero(t <= slots[i, k])[-1]
            assert hist[i, k] == v[last]
            assert t[last] <= slots[i, k] < (t[last + 1] if last + 1 < len(t) else np.inf)
    assert n_nohist == 6 + 6 - 1 and len(tidx) == 120 - n_nohist, "nothing but the series start is dropped"


def test_zoh_fallback_flags_exactly_the_windows_the_primary_rule_drops():
    ts, v = one_gap(100, 7, 100)                     # interval of 8 steps = 40 min > 30 min
    slots, hist, tidx, maxb, hold, _ = build_windows_zoh(ts, v, 6, 6)
    fb = maxb > MAX_INTERP_GAP_S
    # a window has a slot inside the gap iff the primary rule drops it: H-1+min(g,h) targets (D-016)
    assert int(fb.sum()) == 6 - 1 + min(7, 6)
    b = build_windows(ts, v, 6, 6)
    assert int(fb.sum()) == int((b.status == DROPPED_LONG_GAP).sum())
    assert hold.max() == 7 * STEP_S, "the deepest slot is held for the whole gap"
    ts2, v2 = one_gap(100, 5, 100)                    # 30 min: held, but not a fallback
    _, _, _, maxb2, hold2, _ = build_windows_zoh(ts2, v2, 6, 6)
    assert (maxb2 <= MAX_INTERP_GAP_S).all() and (hold2 > 0).any()


def test_primary_variant_reports_when_it_used_a_reading_after_the_origin():
    """A gap of 3 missing slots just before the target: the history-end slot sits inside it and
    the primary rule interpolates from the reading after the origin (D-027 discrepancy)."""
    ts, v = one_gap(50, 3, 50)
    b = build_windows(ts, v, 6, 6)
    assert b.uses_post_origin.any()
    kept_post = b.uses_post_origin[b.kept]
    assert int(kept_post.sum()) == 3          # the min(g, h) targets whose origin is inside the gap
    ts0, v0 = series(np.arange(100))
    assert not build_windows(ts0, v0, 6, 6).uses_post_origin.any()


@needs_data
@pytest.mark.parametrize("pid", PATIENTS)
def test_bglp_variant_predicts_every_challenge_test_point(pid):
    rec = load_patient(pid, "test")
    expected = len(rec.cgm) - rec.eval_start_index           # Table 2's test count
    for h in HORIZONS:
        for H in HISTORY_LENGTHS:
            ws = windows_for(pid, "test", h, H, variant="bglp")
            assert len(ws) == expected and ws.n_dropped_long_gap == 0 and ws.n_dropped_no_history == 0
            assert np.isin(ws.target_ts, rec.cgm.ts).all()
            assert (ws.target_ts - ws.history_end_ts == np.timedelta64(h * STEP_S, "s")).all()
            assert ws.max_hold_s is not None and (ws.max_hold_s >= 0).all()
            # the fallback set is exactly what the primary variant dropped for a long gap
            prim = windows_for(pid, "test", h, H, variant="primary")
            assert int(ws.fallback.sum()) == prim.n_dropped_long_gap
    # the six 2020 contributors' counts, as the challenge states them
    if rec.cohort == "2020":
        assert expected == {"540": 2884, "544": 2704, "552": 2352, "567": 2377, "584": 2653, "596": 2731}[pid]
