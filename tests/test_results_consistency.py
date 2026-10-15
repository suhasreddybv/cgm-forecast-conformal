"""The README's headline sequence-model numbers are read from the committed results files,
not typed, and these tests keep them that way (D-038)."""
from __future__ import annotations

import csv
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _rows(name):
    return list(csv.DictReader(open(REPO / "results" / name)))


def test_final_model_headline_matches_the_steps_table():
    steps = _rows("sequence_steps.csv")
    final = {r["horizon_min"]: float(r["rmse_mean"]) for r in steps
             if r["step"].startswith("7e") and r["variant"] == "bglp" and r["cohort"] == "2020"}
    readme = (REPO / "README.md").read_text()
    assert f"**{final['30']:.2f} / {final['60']:.2f} mg/dL RMSE at 30 / 60 minutes**" in readme
    rank = next(r for r in _rows("sequence_ranking.csv") if r["step"].startswith("7e"))
    assert abs(float(rank["four_score_sum"]) - sum(float(rank[k]) for k in ("rmse30", "mae30", "rmse60", "mae60"))) < 0.01
    assert f"four-score sum **{float(rank['four_score_sum']):.2f}**" in readme
    assert f"**{rank['position_in_official_ranking']}th of nine**" in readme


def test_every_step_has_both_variants_and_the_2020_cohort():
    steps = _rows("sequence_steps.csv")
    labels = {r["step"] for r in steps}
    for label in labels:
        for variant in ("primary", "bglp"):
            assert any(r["step"] == label and r["variant"] == variant and r["cohort"] == "2020" for r in steps), (label, variant)
    assert len(labels) == 12


def test_no_step_scored_below_the_leakage_recheck_threshold():
    steps = _rows("sequence_steps.csv")
    lowest = min(float(r["rmse_mean"]) for r in steps if r["variant"] == "bglp" and r["cohort"] == "2020" and r["horizon_min"] == "30")
    assert lowest >= 18.0, "a step below 18.0 conforming requires the Part G leakage re-check before the README"


def test_sixty_minute_error_exceeds_thirty_for_every_patient_step_and_variant():
    per = _rows("sequence_per_patient.csv")
    key = lambda r: (r["step"], r["variant"], r["patient"])  # noqa: E731
    by = {}
    for r in per:
        by.setdefault(key(r), {})[r["horizon_min"]] = float(r["rmse_mean"])
    bad = [k for k, v in by.items() if not v["60"] > v["30"]]
    assert not bad, bad


def test_selection_was_made_on_validation_only():
    import json
    sel = json.load(open(REPO / "results" / "sequence_selection.json"))
    for step in ("step2", "step3", "step6", "step6fh", "step7"):
        assert step in sel and len(sel[step]) == 24
        for v in sel[step].values():
            assert "val_rmse" in v and not any("test" in k or "bglp" in k or "primary" in k for k in v)


def test_results_files_carry_no_patient_level_timestamps_or_values():
    for name in ("sequence_steps.csv", "sequence_per_patient.csv", "sequence_paired.csv", "clinical_final.csv", "best_epochs.csv"):
        cols = _rows(name)[0].keys()
        assert not any(set(c.split("_")) & {"ts", "timestamp", "glucose"} for c in cols), name
    runs = (REPO / "results" / "sequence_runs.jsonl").read_text().splitlines()
    assert runs and not re.search(r'"\w*(ts|timestamp)\w*":', runs[0])


def test_conformal_first_pass_numbers_are_in_the_files_and_the_readme():
    rows = _rows("coverage.csv")
    pooled = {(r["method"], r["level"], r["horizon_min"]): float(r["coverage"]) for r in rows
              if r["calibration"] == "validation" and r["variant"] == "bglp" and r["cohort"] == "2020" and r["patient"] == ""}
    for m in ("c0", "c1", "c2", "c3"):
        for h in ("30", "60"):
            assert abs(pooled[(m, "0.9", h)] - 0.90) < 0.02, (m, h, pooled[(m, "0.9", h)])
    readme = (REPO / "README.md").read_text()
    assert f"**{pooled[('c0', '0.9', '30')]:.3f} / {pooled[('c0', '0.9', '60')]:.3f} at the 90% level**" in readme
    # the stop condition that fired is recorded, and c3 repairs both sets
    def pat(m, h, pid):
        return float(next(r for r in rows if (r["method"], r["level"], r["calibration"], r["variant"], r["horizon_min"], r["patient"]) == (m, "0.9", "validation", "bglp", h, pid))["coverage"])
    assert pat("c0", "30", "596") < 0.85 and pat("c0", "60", "540") < 0.85
    assert pat("c3", "30", "596") > 0.89 and pat("c3", "60", "540") > 0.89
    onset = _rows("onset_first_pass.csv")
    persist = [r for r in onset if r["detector"] == "persistence"]
    assert persist and all(int(r["tp"]) == 0 and int(r["fp"]) == 0 for r in persist), "persistence detects no onset by construction"
    assert all("sensitivity" in r and "specificity" in r and "ppv" in r and "false_alarms_per_day" in r for r in onset)


def test_where_coverage_breaks_cites_only_reportable_cells_and_matches_the_csv():
    """Every coverage value the section cites by bin resolves to a CSV cell with n >= 100 (D-053)."""
    rows = _rows("coverage_breakdown.csv")
    readme = (REPO / "README.md").read_text()
    sec = readme[readme.index("## Where coverage breaks"):readme.index("## Baselines")]
    def cell(m, h, dim, b):
        return next(r for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["cohort"], r["dimension"], r["bin"]) == (m, "0.9", "bglp", h, "2020", dim, b))
    cited = [("c0", "30", "true_glucose", "54-70"), ("c0", "60", "true_glucose", "54-70"), ("c0", "30", "true_glucose", ">250"),
             ("c0", "60", "true_glucose", ">250"), ("c0", "60", "true_glucose", "180-250"), ("c3", "30", "true_glucose", ">250"),
             ("c3", "60", "true_glucose", ">250"), ("c2", "60", "forecast_glucose", ">250"), ("c2", "60", "true_glucose", ">250"),
             ("c0", "30", "rate_of_change", "2-3"), ("c0", "60", "rate_of_change", "2-3"), ("c0", "30", "rate_of_change", ">3"),
             ("c0", "60", "rate_of_change", ">3"), ("c1", "30", "rate_of_change", ">3"), ("c1", "60", "rate_of_change", ">3")]
    for m, h, dim, b in cited:
        r = cell(m, h, dim, b)
        assert int(r["n"]) >= 100 and r["reportable"] == "True", (m, h, dim, b)
        assert f"{float(r['coverage']):.3f}" in sec, (m, h, dim, b, r["coverage"])
    # the one bin the section says it does not cite
    below54 = cell("c0", "30", "true_glucose", "<54")
    assert int(below54["n"]) < 100 and f"{float(below54['coverage']):.3f}" not in sec
    # the finding rule: a cited finding must exceed 5 points in BOTH variants
    def primary(m, h, dim, b):
        return next(r for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["cohort"], r["dimension"], r["bin"]) == (m, "0.9", "primary", h, "2020", dim, b))
    for m, h, dim, b in (("c0", "30", "true_glucose", ">250"), ("c0", "60", "true_glucose", ">250"), ("c0", "30", "rate_of_change", ">3")):
        assert float(cell(m, h, dim, b)["coverage"]) < 0.85 and float(primary(m, h, dim, b)["coverage"]) < 0.85, (m, h, dim, b)
    assert "n ≥ 100" in sec and "D-053" in sec


def test_icgm_analogues_table_matches_the_readme_and_flags_small_cells():
    rows = _rows("icgm_analogues.csv")
    readme = (REPO / "README.md").read_text()
    assert "### iCGM analogues" in readme and "analogue" in rows[0]["criterion"]
    pf = {(r["horizon_min"], r["criterion"][-2]): r for r in rows if r["estimator"] == "final model, point forecast"}
    assert pf[("60", "A")]["reportable"] == "False" and int(pf[("60", "A")]["n"]) < 100
    assert pf[("30", "A")]["reportable"] == "True" and pf[("30", "A")]["met"] == "False"
    assert f"{float(pf[('30', 'A')]['wilson_lower_95']):.0f}% within 15 mg/dL against the 85% required" in readme
    for r in rows:
        if r["reportable"] == "False":
            assert "not reportable" in r["note"]
    # c4 appears as a lower bound only
    assert any(r["estimator"] == "c4 lower bound, 90%" for r in rows) and not any(r["estimator"].startswith("c4 upper") for r in rows)


def test_c4_is_tighter_and_covers_one_sided_at_nominal_in_the_mid_range():
    rows = _rows("c4_hypoglycaemic_bins.csv")
    def cell(m, v, h, b):
        return next(r for r in rows if (r["method"], r["level"], r["variant"], r["horizon_min"], r["true_bin"]) == (m, "0.9", v, h, b))
    for v in ("bglp", "primary"):
        for h in ("30", "60"):
            mid = {m: cell(m, v, h, "70-180") for m in ("c0", "c3", "c4")}
            assert abs(float(mid["c4"]["one_sided_coverage"]) - 0.9) < 0.025
            assert float(mid["c0"]["one_sided_coverage"]) > float(mid["c4"]["one_sided_coverage"])
            assert float(mid["c4"]["median_distance_forecast_to_bound"]) < float(mid["c0"]["median_distance_forecast_to_bound"])
            low = {m: cell(m, v, h, "54-70") for m in ("c0", "c4")}
            assert int(low["c4"]["miss_below"]) >= int(low["c0"]["miss_below"])
