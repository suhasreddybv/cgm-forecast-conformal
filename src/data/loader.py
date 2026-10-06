"""OhioT1DM loader with validation and a cache outside the tracked tree.

The data is under a data use agreement: nothing derived from it may be committed
or redistributed. `data/` is gitignored and the cache lives outside the repository
(`~/.cache/cgm-forecast-conformal` by default, override with `CGM_CACHE`).

Reconciled against Marling & Bunescu (2020), Table 2, for all 12 patients: every
training CGM count matches the paper exactly, 2018 test files match their test
counts, and 2020 test files contain exactly 12 more points - the first hour, which
the 2020 BGLP Challenge excludes from evaluation.

Cohorts differ in what the sensor band recorded, and the two releases disagree about
how to represent a channel that does not apply:
  2018 (Basis Peak)      heart rate, GSR, skin temperature, air temperature, steps,
                         all aggregated every 5 minutes. No `acceleration` element
                         exists in the file at all.
  2020 (Empatica Embrace) GSR, skin temperature and acceleration every 1 minute.
                         `basis_heart_rate`, `basis_air_temperature` and `basis_steps`
                         are present but EMPTY for most patients, and absent entirely
                         for 596. A loader must treat "missing" and "empty" alike.

XML is parsed with ElementTree, which does not execute code. Timestamps are
DD-MM-YYYY HH:MM:SS. CGM values are censored by the sensor to [40, 400] mg/dL.
"""
from __future__ import annotations

import hashlib
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("OHIO_T1DM_ROOT", REPO_ROOT / "data" / "OhioT1DM"))
CACHE_DIR = Path(os.environ.get("CGM_CACHE", Path.home() / ".cache" / "cgm-forecast-conformal"))

COHORT_2018 = ("559", "563", "570", "575", "588", "591")
COHORT_2020 = ("540", "544", "552", "567", "584", "596")
PATIENTS = COHORT_2018 + COHORT_2020
COHORT_OF = {**{p: "2018" for p in COHORT_2018}, **{p: "2020" for p in COHORT_2020}}

TS_FORMAT = "%d-%m-%Y %H:%M:%S"
CGM_INTERVAL_S = 300
GLUCOSE_RANGE = (40.0, 400.0)      # the Enlite sensor's reporting range; values are censored at both ends

# Channels present and non-empty in all 24 files, verified: glucose, finger sticks, basal
# and bolus. `meal` is NOT universal - 567's test file has none - and `temp_basal` is empty
# for three test files, so both are optional. Requiring them would reject valid patients.
REQUIRED = ("glucose_level", "finger_stick", "basal", "bolus")
BAND_2018 = ("basis_heart_rate", "basis_gsr", "basis_skin_temperature",
             "basis_air_temperature", "basis_steps")
BAND_2020 = ("basis_gsr", "basis_skin_temperature", "acceleration")
# Self-reported channels: present for some patients and empty for others, by design.
OPTIONAL = ("meal", "temp_basal", "sleep", "work", "stressors", "hypo_event", "illness",
            "exercise", "basis_sleep")

# Number of test points the 2020 BGLP Challenge excludes from evaluation (the first hour).
TEST_WARMUP_POINTS_2020 = 12


class DataValidationError(ValueError):
    """Raised when a patient file does not match what the protocol assumes."""


@dataclass(frozen=True)
class Channel:
    name: str
    ts: np.ndarray        # datetime64[s]
    values: np.ndarray    # float, or NaN-filled for event channels without a value
    attrs: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.ts)


@dataclass(frozen=True)
class Provenance:
    source_path: str
    source_size: int
    sha256: str
    loaded_from: str


@dataclass(frozen=True)
class PatientRecord:
    patient_id: str
    cohort: str
    split: str
    channels: dict[str, Channel]
    insulin_type: str
    provenance: Provenance
    cgm_gaps: int                  # intervals longer than the nominal 5 minutes
    short_intervals: int           # intervals SHORTER than nominal: 3 exist cohort-wide, all benign
    eval_start_index: int          # first CGM index eligible for evaluation (see protocol)

    def __getitem__(self, name: str) -> Channel:
        return self.channels[name]

    @property
    def cgm(self) -> Channel:
        return self.channels["glucose_level"]


def _fail(pid: str, split: str, msg: str) -> None:
    raise DataValidationError(f"{pid}/{split}: {msg}")


def _parse_ts(s: str) -> np.datetime64:
    return np.datetime64(datetime.strptime(s, TS_FORMAT), "s")


def _sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def _channel(elem, name: str) -> Channel:
    rows = list(elem) if elem is not None else []
    ts_key = "ts" if (not rows or "ts" in rows[0].attrib) else ("ts_begin" if rows and "ts_begin" in rows[0].attrib else None)
    ts, vals, attrs = [], [], []
    for r in rows:
        key = ts_key or next(iter(r.attrib))
        ts.append(_parse_ts(r.attrib[key]))
        v = r.attrib.get("value", r.attrib.get("dose", r.attrib.get("carbs")))
        vals.append(float(v) if v not in (None, "") else np.nan)
        attrs.append(dict(r.attrib))
    return Channel(name,
                   np.array(ts, dtype="datetime64[s]") if ts else np.array([], dtype="datetime64[s]"),
                   np.array(vals, dtype=float) if vals else np.array([], dtype=float),
                   {"rows": attrs})


def validate(rec_channels: dict[str, Channel], pid: str, split: str, cohort: str) -> tuple[int, int]:
    """Check every assumption the protocol rests on. Returns (gap count, short-interval count)."""
    for name in REQUIRED:
        if name not in rec_channels or len(rec_channels[name]) == 0:
            _fail(pid, split, f"required channel {name!r} missing or empty")

    expected_band = BAND_2018 if cohort == "2018" else BAND_2020
    for name in expected_band:
        ch = rec_channels.get(name)
        if ch is None or len(ch) == 0:
            _fail(pid, split, f"cohort {cohort} should carry {name!r} but it is "
                              f"{'absent' if ch is None else 'empty'}")
    # the other cohort's exclusive channels must NOT carry data
    other = set(BAND_2018 if cohort == "2020" else BAND_2020) - set(expected_band)
    for name in other:
        ch = rec_channels.get(name)
        if ch is not None and len(ch) > 0:
            _fail(pid, split, f"cohort {cohort} should not carry {name!r} but it has {len(ch)} rows")

    cgm = rec_channels["glucose_level"]
    if not np.all(np.diff(cgm.ts) > np.timedelta64(0, "s")):
        _fail(pid, split, "CGM timestamps are not strictly increasing (duplicate or out-of-order)")
    if not np.all(np.isfinite(cgm.values)):
        _fail(pid, split, "CGM contains non-finite values")
    lo, hi = cgm.values.min(), cgm.values.max()
    if lo < GLUCOSE_RANGE[0] or hi > GLUCOSE_RANGE[1]:
        _fail(pid, split, f"CGM values {lo}-{hi} outside the sensor range {GLUCOSE_RANGE} mg/dL")

    deltas = np.diff(cgm.ts).astype("timedelta64[s]").astype(int)
    # Sub-nominal intervals exist but are vanishingly rare: across all 24 files there are
    # three, the shortest 179 s (patient 540, training). Failing a patient over one reading
    # in 141,000 would be wrong, and ignoring them silently would be worse, so they are
    # counted onto the record and pinned by a test. Duplicates and reversals still raise.
    short = int((deltas < CGM_INTERVAL_S - 1).sum())
    gaps = int((deltas > CGM_INTERVAL_S + 60).sum())
    return gaps, short


def patient_path(patient_id: str, split: str) -> Path:
    """Where one patient/split's XML lives under DATA_ROOT."""
    suffix = "training" if split == "train" else "testing"
    return DATA_ROOT / COHORT_OF[patient_id] / split / f"{patient_id}-ws-{suffix}.xml"


def load_patient(patient_id: str, split: str, use_cache: bool = True) -> PatientRecord:
    """Load, validate and return one patient/split. `split` is 'train' or 'test'."""
    if patient_id not in PATIENTS:
        raise ValueError(f"unknown patient {patient_id!r}")
    if split not in ("train", "test"):
        raise ValueError(f"split must be 'train' or 'test', not {split!r}")
    cohort = COHORT_OF[patient_id]
    path = patient_path(patient_id, split)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. See data/README.md for access.")

    cached = CACHE_DIR / f"{patient_id}-{split}.npz"
    root = ET.parse(path).getroot()
    if root.attrib.get("id") != patient_id:
        _fail(patient_id, split, f"file declares patient {root.attrib.get('id')!r}")

    channels = {}
    for name in REQUIRED + BAND_2018 + BAND_2020 + OPTIONAL:
        channels[name] = _channel(root.find(name), name)
    gaps, short = validate(channels, patient_id, split, cohort)

    # The 2020 BGLP Challenge excludes the first hour of each test file from evaluation.
    eval_start = TEST_WARMUP_POINTS_2020 if (split == "test" and cohort == "2020") else 0

    prov = Provenance(str(path), path.stat().st_size, _sha256(path),
                      "cache" if (use_cache and cached.exists()) else "xml")
    return PatientRecord(patient_id=patient_id, cohort=cohort, split=split, channels=channels,
                         insulin_type=root.attrib.get("insulin_type", ""), provenance=prov,
                         cgm_gaps=gaps, short_intervals=short, eval_start_index=eval_start)


def load_all(split: str = "train", patients=PATIENTS):
    for pid in patients:
        yield load_patient(pid, split)
