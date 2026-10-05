import numpy as np
import pytest

from src.data import loader
from src.data.loader import (
    BAND_2018,
    BAND_2020,
    COHORT_2018,
    COHORT_2020,
    PATIENTS,
    DataValidationError,
    load_patient,
)

# Table 2, Marling & Bunescu (2020). Training counts, and test counts as the challenge
# scores them: the 2020 files carry 12 extra points (the first hour) that are excluded.
TABLE_2 = {"540": (11947, 2884), "544": (10623, 2704), "552": (9080, 2352),
           "567": (10858, 2377), "584": (12150, 2653), "596": (10877, 2731),
           "559": (10796, 2514), "563": (12124, 2570), "570": (10982, 2745),
           "575": (11866, 2590), "588": (12640, 2791), "591": (10847, 2760)}

have_data = (loader.DATA_ROOT / "2018" / "train" / "559-ws-training.xml").exists()
needs_data = pytest.mark.skipif(not have_data, reason="OhioT1DM not present; see data/README.md")


def test_twelve_patients_in_two_cohorts_of_six():
    assert len(PATIENTS) == 12
    assert len(COHORT_2018) == len(COHORT_2020) == 6
    assert not set(COHORT_2018) & set(COHORT_2020)


def test_unknown_patient_or_split_rejected():
    with pytest.raises(ValueError, match="unknown patient"):
        load_patient("999", "train")
    with pytest.raises(ValueError, match="split must be"):
        load_patient("559", "validation")


@needs_data
@pytest.mark.parametrize("pid", PATIENTS)
def test_cgm_counts_match_the_published_table(pid):
    train, test = TABLE_2[pid]
    assert len(load_patient(pid, "train").cgm) == train
    rec = load_patient(pid, "test")
    # 2020 files carry the first hour, which the challenge excludes from evaluation
    extra = loader.TEST_WARMUP_POINTS_2020 if rec.cohort == "2020" else 0
    assert len(rec.cgm) == test + extra
    assert rec.eval_start_index == extra
    assert len(rec.cgm) - rec.eval_start_index == test


@needs_data
def test_cohorts_carry_different_band_channels():
    """2018 wore the Basis Peak, 2020 the Empatica Embrace; the channel sets differ."""
    a, b = load_patient("559", "train"), load_patient("540", "train")
    for name in BAND_2018:
        assert len(a[name]) > 0, name
    for name in BAND_2020:
        assert len(b[name]) > 0, name
    # heart rate, air temperature and steps exist only for the Basis cohort
    for name in ("basis_heart_rate", "basis_air_temperature", "basis_steps"):
        assert len(a[name]) > 0 and len(b[name]) == 0, name
    # acceleration only for the Empatica cohort
    assert len(a["acceleration"]) == 0 and len(b["acceleration"]) > 0
    # the Empatica band samples faster, so it has more rows than the 5-minute CGM
    assert len(b["basis_gsr"]) > 2 * len(b.cgm)


@needs_data
@pytest.mark.parametrize("pid", ["559", "540"])
def test_cgm_timestamps_are_strictly_increasing(pid):
    for split in ("train", "test"):
        ts = load_patient(pid, split).cgm.ts
        assert np.all(np.diff(ts) > np.timedelta64(0, "s"))


@needs_data
@pytest.mark.parametrize("pid", PATIENTS)
def test_train_and_test_are_disjoint_in_time(pid):
    """The split is chronological: no test reading may precede the end of training."""
    tr, te = load_patient(pid, "train"), load_patient(pid, "test")
    assert te.cgm.ts[0] > tr.cgm.ts[-1], f"{pid}: test starts before training ends"


@needs_data
def test_glucose_values_sit_inside_the_sensor_range():
    rec = load_patient("559", "train")
    assert rec.cgm.values.min() >= loader.GLUCOSE_RANGE[0]
    assert rec.cgm.values.max() <= loader.GLUCOSE_RANGE[1]


@needs_data
def test_short_intervals_are_counted_not_ignored():
    """Three sub-nominal intervals exist cohort-wide; they are recorded, not silently dropped."""
    total = sum(load_patient(p, s).short_intervals for p in PATIENTS for s in ("train", "test"))
    assert total == 3
    assert load_patient("540", "train").short_intervals == 1


@needs_data
def test_every_file_loads_and_validates():
    for pid in PATIENTS:
        for split in ("train", "test"):
            rec = load_patient(pid, split)
            assert rec.patient_id == pid and rec.split == split
            assert rec.provenance.sha256 and rec.insulin_type
            assert rec.cgm_gaps >= 0


@needs_data
def test_validation_rejects_a_cohort_mismatch():
    """A 2020 patient carrying Basis-only channels would mean the cohorts were mixed up."""
    rec = load_patient("540", "train")
    channels = dict(rec.channels)
    channels["basis_heart_rate"] = load_patient("559", "train")["basis_heart_rate"]
    with pytest.raises(DataValidationError, match="should not carry"):
        loader.validate(channels, "540", "train", "2020")
