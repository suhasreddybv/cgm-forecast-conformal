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
