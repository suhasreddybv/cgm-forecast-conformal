"""The clinical metrics (D-022, D-023), pinned before any model exists.

Grid boundaries are not adjusted to make a test pass: every expected zone below was
worked out from the published definitions first (Clarke: Boren & Clarke 2010 text and
Figure 1; Parkes: Pfuetzner 2013 Table 1), and the test says which side of a boundary a
point is expected to fall on and why.
"""
from __future__ import annotations

import numpy as np
import pytest

from src.eval.clinical import (
    PARKES_T1,
    clarke,
    clarke_zones,
    hypo_detection,
    mard,
    parkes,
    parkes_zones,
    score,
    score_groups,
    tir_agreement,
)

# Pfuetzner et al. 2013, Table 1, as fetched from PMC3876371: the two type 1 columns of
# each zone block, x/y exactly as printed. The constants in the module must equal this.
PFUETZNER_TABLE_1_TYPE1 = """
B Lower: 50/0 50/30 170/145 385/300 550/450
B Upper: 0/50 30/50 140/170 280/380 430/550
C Lower: 120/0 120/30 260/130 550/250
C Upper: 0/60 30/60 50/80 70/110 260/550
D Lower: 250/0 250/40 550/150
D Upper: 0/100 25/100 50/125 80/215 125/550
E Upper: 0/150 35/155 50/550
"""


def test_parkes_vertices_match_the_published_table_transcription():
    parsed = {}
    for line in PFUETZNER_TABLE_1_TYPE1.strip().splitlines():
        name, pts = line.split(":")
        key = name.strip().replace(" ", "_").replace("Lower", "lower").replace("Upper", "upper")
        parsed[key] = [tuple(int(v) for v in p.split("/")) for p in pts.split()]
    assert parsed == PARKES_T1


# ---------------------------------------------------------------- Clarke zones

CLARKE_CASES = [
    # (reference, estimate, zone, why)
    (100, 110, "A", "within 20%"),
    (100, 120, "A", "exactly 20%: the 20% line is inclusive (within)"),
    (100, 121, "B", "just outside 20%"),
    (60, 65, "A", "both <= 70"),
    (70, 70, "A", "both at the 70 threshold, inclusive"),
    (150, 200, "B", "33% high but below the upper C line (150+110 = 260)"),
    (100, 185, "B", "above 180 but below the upper C line at ref 100 (210): the common 'ref>=70 and est>=180 is C' rule is wrong here"),
    (100, 211, "C", "above the upper C line y = x + 110"),
    (100, 210, "B", "exactly on the upper C line: boundary goes to the less severe zone"),
    (71, 182, "C", "just right of ref 70, just above the line (71+110 = 181)"),
    (70, 182, "E", "ref <= 70 and est > 180 is upper E, not C: the E/C edge is at ref 70"),
    (70, 180, "D", "ref <= 70 and est in (70, 180]: 180 is still D, E starts above it"),
    (60, 100, "D", "upper D: missed hypoglycaemia"),
    (60, 73, "D", "ref <= 70, est above 70 and 13 > 20% of 60 (12): upper D"),
    (60, 72, "A", "exactly 20% of 60 above it: within 20% takes precedence over the upper-D box"),
    (300, 150, "D", "lower D: ref >= 240, est in target range"),
    (240, 180, "D", "lower D corners inclusive"),
    (239, 180, "B", "ref below 240: lower B"),
    (200, 60, "E", "lower E: ref >= 180 and est <= 70"),
    (180, 70, "E", "lower E corner inclusive"),
    (179, 68, "C", "lower C triangle: est < 1.4*(179-130) = 68.6"),
    (150, 27, "C", "lower C: est < 1.4*(150-130) = 28"),
    (150, 28, "B", "exactly on the lower C line: boundary goes to the less severe zone"),
    (130, 0, "B", "the triangle's apex is on the line"),
    (120, 10, "B", "left of the lower C triangle: lower B"),
    (50, 200, "E", "upper E"),
    (300, 250, "A", "within 20% (240 is the bound)"),
    (300, 230, "B", "lower B"),
]


@pytest.mark.parametrize("ref,est,zone,why", CLARKE_CASES, ids=[c[3][:40] for c in CLARKE_CASES])
def test_clarke_hand_placed_points(ref, est, zone, why):
    assert clarke_zones([ref], [est])[0] == zone, why


def test_clarke_every_zone_is_reachable():
    zones = set(clarke_zones([c[0] for c in CLARKE_CASES], [c[1] for c in CLARKE_CASES]))
    assert zones == {"A", "B", "C", "D", "E"}


# ---------------------------------------------------------------- Parkes zones
# Boundary values below are computed by hand from the Table 1 vertices.

PARKES_CASES = [
    (100, 100, "A", "identity"),
    (20, 5, "A", "left of the vertical B-lower segment at x=50: still A"),
    (50, 10, "A", "on the vertical B-lower segment: boundary goes to the less severe zone"),
    (51, 10, "B", "just right of it, below the line"),
    (100, 125, "A", "B-upper at x=100 is y = 50 + 70*(120/110) = 126.4; 125 is below it"),
    (100, 127, "B", "just above the B-upper line"),
    (100, 70, "B", "B-lower at x=100 is y = 30 + 50*(115/120) = 77.9; 70 is below it"),
    (100, 78, "A", "just above the B-lower line"),
    (20, 70, "C", "C-upper: for y in (60, 80] the line is x = 30 + (y - 60) = 40; x=20 is left of it"),
    (40, 70, "B", "exactly on the C-upper line: less severe zone (B-upper at y=70: x = 30 + 20*110/120 = 48.3, so B)"),
    (35, 65, "B", "between B-upper (43.8) and C-upper (35 -> on the line)"),
    (130, 20, "C", "C-lower: x > 120 and y below the line from (120,30)"),
    (120, 10, "B", "on the vertical C-lower segment at x=120: B"),
    (300, 30, "D", "D-lower: x > 250 and y below (250,40)-(550,150) at x=300: 58.3"),
    (300, 59, "C", "just above the D-lower line, below C-lower (x=300: 130 + 40*120/290 = 146.6)"),
    (10, 120, "D", "D-upper: y in (100, 125] -> x = 25 + (y-100) = 45; x=10 is left of it"),
    (45, 120, "C", "exactly on the D-upper line: less severe (C-upper at y=120: 70 + 10*190/440 = 74.3)"),
    (10, 200, "E", "E-upper: y > 155 -> x = 35 + (200-155)*15/395 = 36.7; 10 is left of it"),
    (37, 200, "D", "just right of the E-upper line, left of D-upper at y=200 (50 + 75*30/90 = 75)"),
    (0, 150, "D", "E-upper starts at y=150 exclusive; (0,150) is on the line -> D-upper (y > 100, x < 25)"),
    (0, 151, "E", "just above the E-upper start"),
    (600, 600, "A", "beyond the 550 grid: identity stays A"),
    (600, 300, "C", "beyond 550: C-lower extended, at x=600 y = 250 + 50*120/290 = 270.7; 300 above it -> B? no: B-lower at 600 = 450 + 50*150/165 = 495.5, so below B-lower and above C-lower -> B"),
]
# the last case's note shows the working; the expected zone is corrected to what it derives
PARKES_CASES[-1] = (600, 300, "B", "beyond 550, between the extended C-lower (270.7) and B-lower (495.5) lines")


@pytest.mark.parametrize("ref,est,zone,why", PARKES_CASES, ids=[c[3][:40] for c in PARKES_CASES])
def test_parkes_hand_placed_points(ref, est, zone, why):
    assert parkes_zones([ref], [est])[0] == zone, why


def test_parkes_every_zone_is_reachable():
    zones = set(parkes_zones([c[0] for c in PARKES_CASES], [c[1] for c in PARKES_CASES]))
    assert zones == {"A", "B", "C", "D", "E"}


def test_the_two_grids_disagree_on_known_points():
    """Clarke is harsher at the 70 edge; Parkes is harsher far from the identity line."""
    ref, est = [60, 100, 300], [75, 125, 30]
    assert clarke_zones(ref, est).tolist() == ["D", "B", "E"]
    assert parkes_zones(ref, est).tolist() == ["A", "A", "D"]


# ----------------------------------------------------------------- the rest

def test_perfect_forecast():
    y = np.array([40, 55, 70, 100, 180, 250, 400], float)
    s = score(y, y)
    assert s["clarke_zone_A_pct"] == 100.0 and s["parkes_zone_A_pct"] == 100.0
    assert s["mard_pct"] == 0.0 and s["tir_agreement"] == 1.0
    for t in (70, 54):
        assert s[f"hypo{t}_sensitivity"] == 1.0 and s[f"hypo{t}_specificity"] == 1.0
    assert s["n"] == 7


def test_constant_forecast_at_100_zone_distribution_by_hand():
    """Reference 50, 60, 80, 100, 125, 150, 200, 250, 300, 400 against a constant 100.
    Clarke: 50->D (ref<=70, est in (70,180]); 60->D; 80->A (25%? no: |100-80| = 20 = 25% of
    80 -> B); 100->A; 125->A (20%); 150->B (33%); 200->B (50%, below C line: lower C needs
    est < 1.4*(ref-130) = 98 at ref 200? ref must be < 180, so no); 250->D (ref>=240, est in
    range); 300->D; 400->D."""
    ref = np.array([50, 60, 80, 100, 125, 150, 200, 250, 300, 400], float)
    est = np.full(10, 100.0)
    assert clarke_zones(ref, est).tolist() == ["D", "D", "B", "A", "A", "B", "B", "D", "D", "D"]
    c = clarke(ref, est)
    assert (c["zone_A_n"], c["zone_B_n"], c["zone_C_n"], c["zone_D_n"], c["zone_E_n"]) == (2, 3, 0, 5, 0)
    assert c["zone_AB_pct"] == 50.0
    # Parkes by hand: 50 -> B-upper at x=50 is y = 50 + 20*120/110 = 71.8 -> 100 above -> beyond B;
    # C-upper at y=100: x = 50 + 20*20/30 = 63.3 -> 50 < 63.3 -> beyond C; D-upper at y=100 is
    # the horizontal start (y > 100 needed) -> not D -> C. 60 -> C-upper 63.3 -> 60 < 63.3 -> C.
    # 80 -> B-upper at x=80: 50 + 50*120/110 = 104.5 -> 100 below -> A. 400 -> B-lower at 400:
    # 300 + 15*150/165 = 313.6 -> 100 below; C-lower at 400: 130 + 140*120/290 = 187.9 -> below;
    # D-lower at 400: 40 + 150*110/300 = 95 -> 100 above -> C. 300 -> D-lower 58.3, C-lower
    # 146.6 -> C. 250 -> C-lower 125.9 -> below; D-lower: x > 250 needed -> not D -> C.
    # 200 -> C-lower at 200: 30 + 80*100/140 = 87.1 -> 100 above; B-lower at 200: 145 + 30*155/215
    # = 166.6 -> below -> B. 150 -> B-lower at 150: 30 + 100*115/120 = 125.8 -> below -> B.
    # 125 -> B-lower 101.9 -> 100 below -> B. 100 -> A.
    assert parkes_zones(ref, est).tolist() == ["C", "C", "A", "A", "B", "B", "B", "C", "C", "C"]


def test_mard_uses_the_reference_as_denominator():
    r = mard([100, 200], [110, 180])
    assert r["mard_pct"] == 10.0 and r["denominator"] == "reference" and r["n"] == 2


def test_tir_agreement_counts():
    ref = [60, 100, 180, 200]
    est = [80, 100, 181, 250]
    t = tir_agreement(ref, est)
    assert (t["both_in_range"], t["both_out_of_range"], t["ref_in_pred_out"], t["ref_out_pred_in"]) == (1, 1, 1, 1)
    assert t["agreement"] == 0.5 and t["n"] == 4


def test_hypo_detection_confusion_and_censoring_flag():
    ref = [40, 40, 60, 100, 100]
    est = [50, 80, 55, 60, 120]
    h = hypo_detection(ref, est, 70)
    assert (h["tp"], h["fn"], h["fp"], h["tn"]) == (2, 1, 1, 1)
    assert h["sensitivity"] == round(2 / 3, 4) and h["specificity"] == 0.5 and h["ppv"] == round(2 / 3, 4)
    assert h["censored"] is True and h["n_events_at_sensor_floor"] == 2 and "optimistic" in h["caveat"]
    clean = hypo_detection([60, 100], [50, 110], 70)
    assert clean["censored"] is False and clean["caveat"] == ""


def test_every_reference_event_at_the_floor_returns_the_censoring_flag():
    ref = np.full(5, 40.0)
    est = np.array([30, 45, 50, 60, 80], float)
    h = hypo_detection(ref, est, 54)
    assert h["censored"] is True and h["share_of_events_at_floor"] == 1.0
    assert h["sensitivity"] == 0.6 and np.isnan(h["specificity"])      # 30, 45, 50 are below 54


def test_score_groups_carries_n_through_every_level():
    rng = np.random.default_rng(0)
    groups = []
    for pid, cohort in (("a", "2018"), ("b", "2018"), ("c", "2020")):
        y = rng.uniform(40, 400, 50)
        groups.append(dict(patient=pid, cohort=cohort, horizon_min=30, method="m", y_true=y, y_pred=y * 1.1))
    rows = score_groups(groups)
    levels = {(r["level"], r["cohort"], r["patient"]): r for r in rows}
    assert levels[("pooled", "pooled", "")]["n"] == 150 and levels[("cohort", "2018", "")]["n"] == 100
    assert levels[("patient", "2020", "c")]["n"] == 50
    assert all(r["clarke_zone_A_pct"] == 100.0 for r in rows)      # 10% high is within 20%
    assert len(rows) == 3 + 2 + 1
