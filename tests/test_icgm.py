"""Forecast analogues of 21 CFR 862.1355(b) (D-059): the code's thresholds are the
document's, and the criterion logic is checked on hand-built data."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from src.eval.icgm import DOC, accuracy_analogue, exclusion_analogue, load_criteria, rate_analogue, wilson_lower_one_sided

CRIT = load_criteria()


def test_every_parsed_threshold_appears_verbatim_in_the_document():
    text = DOC.read_text()
    assert len(CRIT.accuracy) == 7 and [a.para for a in CRIT.accuracy] == list("ABCDEFG")
    for a in CRIT.accuracy:
        row = re.search(rf"^\|\s*\({a.para}\)\s*\|.*$", text, flags=re.M).group(0)
        assert f"{a.required_pct:g}%" in row
        assert (f"±{a.tol_value:g}%" if a.tol_relative else f"±{a.tol_value:g} mg/dL") in row
        for edge in (a.lo, a.hi):
            if np.isfinite(edge):
                assert f"{edge:g}" in row
    assert f"< {CRIT.exclusion_low:g}" in text and f"> {CRIT.exclusion_high:g}" in text
    assert f"No more than {CRIT.rate_limit_pct:g}%" in text
    assert f"> +{CRIT.rate_indicated:g} mg/dL/min" in text and f"< −{CRIT.rate_true:g} mg/dL/min" in text


def test_ranges_are_binned_by_the_estimate_with_the_documented_edges():
    a, b, c, g = CRIT.accuracy[0], CRIT.accuracy[1], CRIT.accuracy[2], CRIT.accuracy[6]
    assert (a.lo, a.hi, a.hi_inclusive) == (-np.inf, 70.0, False)          # < 70
    assert (b.lo, b.hi, b.lo_inclusive, b.hi_inclusive) == (70.0, 180.0, True, True)
    assert (c.lo, c.lo_inclusive, c.hi) == (180.0, False, np.inf)          # > 180
    assert (g.lo, g.hi) == (-np.inf, np.inf)
    assert not a.tol_relative and b.tol_relative


def test_accuracy_analogue_on_hand_built_pairs():
    a = CRIT.accuracy[0]                                   # < 70, within 15 mg/dL
    ref = np.array([60, 60, 60, 100, 60], float)
    est = np.array([50, 74, 80, 60, 69.9], float)          # est 74/80 are not < 70; 50 within 15 of 60; 60 vs 100 out; 69.9 within
    r = accuracy_analogue(ref, est, a)
    assert (r["n"], r["within"]) == (3, 2) and r["pct"] == pytest.approx(66.67, abs=0.01) and r["reportable"] is False
    b = CRIT.accuracy[1]                                   # 70-180, within 15%
    r = accuracy_analogue(np.array([100, 100, 200.0]), np.array([114, 116, 180.0]), b)
    assert (r["n"], r["within"]) == (3, 2)                 # 114 within 15% of 100, 116 not, 180 vs 200 is 10% -> within
    g = CRIT.accuracy[6]
    r = accuracy_analogue(np.full(500, 100.0), np.full(500, 110.0), g)   # all within 20%
    assert r["pct"] == 100.0 and r["wilson_lower_95"] > g.required_pct and r["met"] is True and r["reportable"]


def test_wilson_lower_bound_is_below_the_point_estimate_and_tightens_with_n():
    assert wilson_lower_one_sided(90, 100) < 90.0
    assert wilson_lower_one_sided(900, 1000) > wilson_lower_one_sided(90, 100)
    assert np.isnan(wilson_lower_one_sided(0, 0))


def test_exclusions_and_rate_rules():
    ref = np.array([200, 60, 100, 100.0]); est = np.array([65, 190, 100, np.inf])
    e = exclusion_analogue(ref, est, CRIT)
    assert (e["H_violations"], e["I_violations"], e["n"]) == (1, 1, 3) and e["H_per_1000"] == pytest.approx(333.333)
    origin = np.array([100, 100, 100.0]); est = np.array([140, 60, 110.0]); ref = np.array([30, 170, 120.0])
    r = rate_analogue(ref, est, origin, 30.0, CRIT)
    # indicated +1.33, -1.33, +0.33; true -2.33, +2.33, +0.67 -> J: 1, K: 1
    assert (r["J_count"], r["K_count"]) == (1, 1) and r["J_met"] is False and r["limit_pct"] == CRIT.rate_limit_pct
