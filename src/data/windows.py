"""(history, target) windows on the frozen protocol (docs/protocol.md, D-005).

One window per real CGM reading that has enough history. For a target at time T and a
horizon of h steps, the history is H slots at T - (h+k)*5 min, k = H-1 ... 0, so the
history end sits exactly h steps before the target and nothing after it is seen. The
value at a slot is the reading at that instant if one exists, else the linear
interpolation of the two real readings bracketing it. A window is dropped when any
bracketing interval is longer than 30 minutes (protocol rule 2) or when a slot precedes
the first reading available. A target is never interpolated: every target is a reading
that exists in the raw stream, and a test pins that.

Why interpolate rather than snap to the nearest reading: the cadence is not exactly
300 s. Measured over all 24 files, 165,482 intervals are exactly 300 s but 503 are
301-360 s, so after one of those every earlier reading is a few seconds off a grid
anchored on the target. Linear interpolation across a 301 s interval moves the value
by at most 1/301 of one step and is exact whenever the slot coincides with a reading.
The window records the widest bracketing interval it rests on (`max_bracket_s`); a
window counts as interpolated across a gap when that exceeds 360 s, the same gap
definition `gap_analysis.py` uses, so a one-second jitter is not reported as a gap.

Test split: the history series is the training file followed by the test file (the
training tail is exactly 300 s before the first test reading for all 12 patients,
measured), and targets start at `eval_start_index` - the first hour of a 2020 test
file serves as history and is never scored (D-003). History therefore never reaches
into another patient's data and never past the target.

Covariates the primary models may use (D-010) are aligned to the same H slots:
  basal  U/h in effect at the slot instant, with any temporary basal applied
         (NaN before the first basal event of a training file)
  bolus  units delivered in the five minutes ending at the slot; an extended
         ("square") bolus is spread uniformly over its delivery interval
  carbs  grams logged in the five minutes ending at the slot
Nothing after the history end enters any of them.

Caches live outside the tree under CACHE_DIR/windows (data use agreement, D-001).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.data.loader import (
    CACHE_DIR,
    CGM_INTERVAL_S,
    Channel,
    PatientRecord,
    _parse_ts,
    load_patient,
)

HORIZONS = (6, 12)               # steps ahead: 30 and 60 minutes
HISTORY_LENGTHS = (6, 12, 24)    # steps of history: 30 minutes, 1 hour, 2 hours
STEP_S = CGM_INTERVAL_S
GAP_S = CGM_INTERVAL_S + 60      # an interval longer than this is a gap (as in gap_analysis.py)
MAX_INTERP_GAP_S = 30 * 60       # protocol rule 2: history may be interpolated across at most 30 min
CACHE_VERSION = 1


@dataclass(frozen=True)
class WindowSet:
    patient_id: str
    split: str
    horizon: int                  # steps
    history_len: int              # steps
    target_ts: np.ndarray         # datetime64[s], [n]; every one is a real reading
    target: np.ndarray            # [n] mg/dL
    history: np.ndarray           # [n, H] mg/dL, oldest first; column H-1 is the history end
    max_bracket_s: np.ndarray     # [n] widest interval between real readings bracketing any slot
    basal: np.ndarray             # [n, H] U/h
    bolus: np.ndarray             # [n, H] U
    carbs: np.ndarray             # [n, H] g
    n_readings: int               # CGM readings in this split's own file
    n_candidates: int             # readings eligible as targets (after the test warm-up)
    n_dropped_long_gap: int       # candidates whose history crossed a gap over 30 min
    n_dropped_no_history: int     # candidates whose history precedes the first reading

    def __len__(self) -> int:
        return len(self.target)

    @property
    def interpolated(self) -> np.ndarray:
        """Windows whose history was interpolated across a gap (protocol rule 3)."""
        return self.max_bracket_s > GAP_S

    @property
    def history_end_ts(self) -> np.ndarray:
        return self.target_ts - np.timedelta64(self.horizon * STEP_S, "s")

    def history_ts(self) -> np.ndarray:
        """[n, H] slot times, oldest first."""
        offs = slot_offsets_s(self.horizon, self.history_len)
        return self.target_ts[:, None] - offs[None, :].astype("timedelta64[s]")


def slot_offsets_s(horizon: int, history_len: int) -> np.ndarray:
    """Seconds before the target of each history slot, oldest first: (h+H-1 ... h) * 300."""
    return (horizon + np.arange(history_len)[::-1]) * STEP_S


# ----------------------------------------------------------------------------- core

KEPT, DROPPED_LONG_GAP, DROPPED_NO_HISTORY = 0, 1, 2


@dataclass(frozen=True)
class Built:
    """What build_windows found, for every candidate target and for the kept windows."""
    candidates: np.ndarray        # [n_cand] reading indices eligible as targets
    status: np.ndarray            # [n_cand] KEPT / DROPPED_LONG_GAP / DROPPED_NO_HISTORY
    max_bracket_s: np.ndarray     # [n_cand] widest bracketing interval over the history slots (-1 if no history)
    end_bracket_s: np.ndarray     # [n_cand] bracketing interval at the history END slot alone
    slots_s: np.ndarray           # [n_kept, H] epoch seconds of the kept windows' slots, oldest first
    history: np.ndarray           # [n_kept, H]
    target_idx: np.ndarray        # [n_kept] reading index of each kept target

    @property
    def kept(self) -> np.ndarray:
        return self.status == KEPT


def build_windows(ts: np.ndarray, values: np.ndarray, horizon: int, history_len: int,
                  eval_from: int = 0) -> Built:
    """Windows over one strictly increasing CGM series.

    `eval_from` is the index of the first reading eligible as a target; readings before
    it serve as history only.
    """
    if horizon < 1 or history_len < 1:
        raise ValueError("horizon and history_len must be positive step counts")
    t = ts.astype("datetime64[s]").astype(np.int64)
    if len(t) and not np.all(np.diff(t) > 0):
        raise ValueError("timestamps must be strictly increasing")
    cand = np.arange(eval_from, len(t))
    slots = t[cand][:, None] - slot_offsets_s(horizon, history_len)[None, :]        # [n, H]

    r = np.searchsorted(t, slots, side="left")          # first reading at or after each slot
    no_hist = (r == 0) & (t[0] > slots)                 # the slot precedes the first reading
    r = np.minimum(r, len(t) - 1)                       # slot < target <= t[-1], so never past the end
    exact = t[r] == slots
    lft = np.where(exact, r, np.maximum(r - 1, 0))
    bracket = t[r] - t[lft]                             # 0 where the slot is a real reading

    drop_no_hist = no_hist.any(axis=1)
    max_bracket = np.where(drop_no_hist, -1, bracket.max(axis=1))
    status = np.where(drop_no_hist, DROPPED_NO_HISTORY,
                      np.where(max_bracket > MAX_INTERP_GAP_S, DROPPED_LONG_GAP, KEPT))
    keep = status == KEPT

    w = np.where(bracket > 0, (slots - t[lft]) / np.where(bracket > 0, bracket, 1), 0.0)
    hist = values[lft] + (values[r] - values[lft]) * w
    return Built(cand, status, max_bracket, np.where(drop_no_hist, -1, bracket[:, -1]),
                 slots[keep], hist[keep], cand[keep])


# ----------------------------------------------------------------------- covariates

def _interval_rows(ch: Channel):
    """(begin_s, end_s, value) for a channel whose rows carry ts_begin / ts_end."""
    rows = ch.attrs.get("rows", [])
    if not rows:
        z = np.zeros(0)
        return z.astype(np.int64), z.astype(np.int64), z
    b = np.array([_parse_ts(r["ts_begin"]).astype(np.int64) for r in rows])
    e = np.array([_parse_ts(r["ts_end"]).astype(np.int64) for r in rows])
    v = np.array([float(r.get("value", r.get("dose"))) for r in rows])
    order = np.argsort(b, kind="stable")
    return b[order], e[order], v[order]


def _cumulative_at(event_s: np.ndarray, amount: np.ndarray, at_s: np.ndarray) -> np.ndarray:
    """Total amount of instantaneous events at or before each time in `at_s`."""
    if len(event_s) == 0:
        return np.zeros(at_s.shape)
    order = np.argsort(event_s, kind="stable")
    cum = np.concatenate([[0.0], np.cumsum(amount[order])])
    return cum[np.searchsorted(event_s[order], at_s, side="right")]


def basal_rate_at(at_s: np.ndarray, basal: Channel, temp_basal: Channel) -> np.ndarray:
    """U/h in effect at each instant: the scheduled rate, overridden by a temporary basal
    while one is active (begin <= t < end; where two overlap, the later-begun one)."""
    bt = basal.ts.astype("datetime64[s]").astype(np.int64)
    order = np.argsort(bt, kind="stable")
    bt, bv = bt[order], basal.values[order]
    i = np.searchsorted(bt, at_s, side="right") - 1
    rate = np.where(i >= 0, bv[np.maximum(i, 0)], np.nan)
    tb, te, tv = _interval_rows(temp_basal)
    if len(tb):
        j = np.searchsorted(tb, at_s, side="right") - 1
        jj = np.maximum(j, 0)
        active = (j >= 0) & (at_s < te[jj])
        rate = np.where(active, tv[jj], rate)
    return rate


def bolus_delivered_at(at_s: np.ndarray, bolus: Channel) -> np.ndarray:
    """Cumulative insulin (U) delivered at or before each instant. A bolus with
    ts_begin == ts_end lands at once; an extended one is spread uniformly over its interval."""
    b, e, d = _interval_rows(bolus)
    inst = e == b
    total = _cumulative_at(b[inst], d[inst], at_s)
    for bb, ee, dd in zip(b[~inst], e[~inst], d[~inst]):
        total = total + dd * np.clip((at_s - bb) / (ee - bb), 0.0, 1.0)
    return total


def align_covariates(slots_s: np.ndarray, basal: Channel, temp_basal: Channel,
                     bolus: Channel, meal: Channel):
    """basal [n,H] U/h at the slot; bolus and carbs [n,H] delivered in (slot-300 s, slot]."""
    rate = basal_rate_at(slots_s, basal, temp_basal)
    bol = bolus_delivered_at(slots_s, bolus) - bolus_delivered_at(slots_s - STEP_S, bolus)
    mt = meal.ts.astype("datetime64[s]").astype(np.int64)
    carbs = _cumulative_at(mt, meal.values, slots_s) - _cumulative_at(mt, meal.values, slots_s - STEP_S)
    return rate, bol, carbs


# ------------------------------------------------------------------- per patient

def _concat(a: Channel, b: Channel) -> Channel:
    return Channel(a.name, np.concatenate([a.ts, b.ts]), np.concatenate([a.values, b.values]),
                   {"rows": a.attrs.get("rows", []) + b.attrs.get("rows", [])})


def _series_for(rec: PatientRecord):
    """The CGM series and event channels a split's windows may draw history from."""
    if rec.split == "train":
        return rec.cgm, {k: rec[k] for k in ("basal", "temp_basal", "bolus", "meal")}, rec.eval_start_index
    prior = load_patient(rec.patient_id, "train")
    chans = {k: _concat(prior[k], rec[k]) for k in ("basal", "temp_basal", "bolus", "meal")}
    return _concat(prior.cgm, rec.cgm), chans, len(prior.cgm) + rec.eval_start_index


def windows_from_record(rec: PatientRecord, horizon: int, history_len: int) -> WindowSet:
    cgm, chans, eval_from = _series_for(rec)
    b = build_windows(cgm.ts, cgm.values, horizon, history_len, eval_from)
    rate, bol, carbs = align_covariates(b.slots_s, chans["basal"], chans["temp_basal"],
                                        chans["bolus"], chans["meal"])
    return WindowSet(rec.patient_id, rec.split, horizon, history_len,
                     target_ts=cgm.ts[b.target_idx], target=cgm.values[b.target_idx],
                     history=b.history, max_bracket_s=b.max_bracket_s[b.kept],
                     basal=rate, bolus=bol, carbs=carbs,
                     n_readings=len(rec.cgm), n_candidates=len(cgm) - eval_from,
                     n_dropped_long_gap=int((b.status == DROPPED_LONG_GAP).sum()),
                     n_dropped_no_history=int((b.status == DROPPED_NO_HISTORY).sum()))


def cache_path(patient_id: str, split: str, horizon: int, history_len: int) -> Path:
    return CACHE_DIR / "windows" / f"{patient_id}-{split}-h{horizon}-H{history_len}-v{CACHE_VERSION}.npz"


_ARRAYS = ("target_ts", "target", "history", "max_bracket_s", "basal", "bolus", "carbs")
_SCALARS = ("n_readings", "n_candidates", "n_dropped_long_gap", "n_dropped_no_history")


def windows_for(patient_id: str, split: str, horizon: int, history_len: int,
                use_cache: bool = True) -> WindowSet:
    """Build (or read back from the out-of-tree cache) one patient/split/horizon/H."""
    path = cache_path(patient_id, split, horizon, history_len)
    if use_cache and path.exists():
        with np.load(path, allow_pickle=False) as z:
            return WindowSet(patient_id, split, horizon, history_len,
                             **{k: z[k] for k in _ARRAYS}, **{k: int(z[k]) for k in _SCALARS})
    ws = windows_from_record(load_patient(patient_id, split), horizon, history_len)
    if use_cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, **{k: getattr(ws, k) for k in _ARRAYS + _SCALARS})
    return ws
