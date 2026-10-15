"""Conformal prediction intervals around the frozen point forecast (Week 5, D-045 to D-048).

Nothing here touches the model: every method takes calibration residuals and a forecast
and returns an interval. Four methods, all without retraining; conformalised quantile
regression is excluded because the model is frozen as a point forecaster.

  c0  split conformal on the absolute residual |y - yhat|: constant width per
      patient-horizon. Assumes exchangeability of calibration and test scores.
  c1  scaled residual |y - yhat| / sigma(history): sigma from the history ONLY - the SD of
      the first differences of the history slots strictly before the forecast origin,
      floored - so the interval widens when the recent trace is volatile. No sample at or
      after the origin enters sigma (asserted by a sentinel test).
  c2  Mondrian by forecast bin: c0 calibrated separately within yhat < 70, 70-180, > 180.
      Coverage conditional on where the forecast sits; the hypoglycaemia question in its
      simplest form. A bin with fewer than MIN_BIN calibration scores falls back to the
      global quantile and says so.
  c4  one-sided lower bound on the signed residual yhat - y: lower = yhat - q(yhat - y),
      no upper bound. Does a signed score reach down where symmetric ones cannot?
  c3  adaptive conformal inference (Gibbs & Candes 2021): c0's score, with the miscoverage
      level updated online from recent miscoverage, alpha_{t+1} = alpha_t + gamma (alpha -
      err_t). No exchangeability assumption; tracks nominal coverage under drift and after
      gaps. gamma chosen on the validation sequence, never on test.

The finite-sample quantile is the ceil((n+1)(1-alpha))/n empirical quantile of the
calibration scores (Vovk; Lei et al. 2018), conservative by construction.

Onset detection (D-042): an onset event is a target below the threshold whose origin
reading - the history-end slot the forecaster saw - was at or above it. Events already in
progress at the origin are excluded: persistence detects those trivially. A detection is
the lower bound (or the point forecast) below the threshold. Reported with specificity,
PPV and false alarms per day, never sensitivity alone.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

LEVELS = (0.80, 0.90, 0.95)
BINS = ((-np.inf, 70.0), (70.0, 180.0), (180.0, np.inf))
BIN_NAMES = ("<70", "70-180", ">180")
MIN_BIN = 30
SIGMA_FLOOR = 2.0            # mg/dL per step: a flat trace still gets a finite scale
USELESS_WIDTH = 100.0        # an interval wider than this is reported as useless
GAMMAS = (0.002, 0.005, 0.01, 0.02, 0.05)


def conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """The (1 - alpha) split-conformal quantile with the finite-sample correction."""
    s = np.sort(np.asarray(scores, float))
    n = len(s)
    if n == 0:
        return float("inf")
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return float("inf") if k > n else float(s[k - 1])


def sigma_from_history(history: np.ndarray, floor: float = SIGMA_FLOOR) -> np.ndarray:
    """Per-window scale from the history slots strictly before the origin: the SD of their
    first differences. history is [n, H], column -1 the origin slot, which is excluded."""
    h = np.asarray(history, float)[:, :-1]
    d = np.diff(h, axis=1)
    return np.maximum(d.std(axis=1, ddof=1) if d.shape[1] > 1 else np.zeros(len(h)), floor)


@dataclass
class Intervals:
    lower: np.ndarray
    upper: np.ndarray
    note: str = ""

    @property
    def width(self) -> np.ndarray:
        return self.upper - self.lower

    def __getitem__(self, idx) -> "Intervals":
        return Intervals(self.lower[idx], self.upper[idx], self.note)


def c0(cal_y, cal_pred, pred, alpha) -> Intervals:
    q = conformal_quantile(np.abs(cal_y - cal_pred), alpha)
    return Intervals(pred - q, pred + q)


def c1(cal_y, cal_pred, cal_history, pred, history, alpha) -> Intervals:
    s_cal = sigma_from_history(cal_history)
    q = conformal_quantile(np.abs(cal_y - cal_pred) / s_cal, alpha)
    s = sigma_from_history(history)
    return Intervals(pred - q * s, pred + q * s)


def bin_of(pred: np.ndarray) -> np.ndarray:
    p = np.asarray(pred, float)
    return np.where(p < 70, 0, np.where(p <= 180, 1, 2))


def c2(cal_y, cal_pred, pred, alpha) -> Intervals:
    scores = np.abs(cal_y - cal_pred)
    q_global = conformal_quantile(scores, alpha)
    cal_bin, test_bin = bin_of(cal_pred), bin_of(pred)
    q = np.full(len(pred), q_global)
    fell_back = []
    for b in range(3):
        m = cal_bin == b
        if m.sum() >= MIN_BIN:
            q[test_bin == b] = conformal_quantile(scores[m], alpha)
        elif (test_bin == b).any():
            fell_back.append(f"{BIN_NAMES[b]} ({int(m.sum())} calibration scores)")
    return Intervals(pred - q, pred + q, "fallback to global quantile in " + ", ".join(fell_back) if fell_back else "")


def c4(cal_y, cal_pred, pred, alpha) -> Intervals:
    """One-sided lower bound on the SIGNED residual yhat - y (positive when the model
    over-predicts): lower = yhat - q_alpha(yhat - y). No upper bound is produced (upper is
    +inf); coverage is one-sided, P(y >= lower), nominal by construction. Added 15 Oct
    because every miss in the 54-70 bin was below the bound for c0 and c3 in both
    variants: a symmetric score cannot correct direction (D-062)."""
    q = conformal_quantile(cal_pred - cal_y, alpha)
    pred = np.asarray(pred, float)
    return Intervals(pred - q, np.full(len(pred), np.inf), "one-sided lower bound")


def c3(cal_y, cal_pred, y, pred, alpha, gamma) -> Intervals:
    """Adaptive conformal inference over a time-ordered test sequence. y is read only AFTER
    each interval is formed, to update the level - the online protocol."""
    scores = np.sort(np.abs(cal_y - cal_pred))
    a_t = alpha
    lower, upper = np.empty(len(pred)), np.empty(len(pred))
    for t in range(len(pred)):
        a_clip = min(max(a_t, 0.0), 1.0)
        q = conformal_quantile(scores, a_clip) if a_clip > 0 else float("inf")
        lower[t], upper[t] = pred[t] - q, pred[t] + q
        err = 0.0 if lower[t] <= y[t] <= upper[t] else 1.0
        a_t = a_t + gamma * (alpha - err)
    return Intervals(lower, upper, f"gamma={gamma}")


def choose_gamma(cal_y, cal_pred, val_y, val_pred, alpha, gammas=GAMMAS) -> float:
    """gamma that brings ACI's coverage on the validation sequence closest to nominal,
    with calibration scores from the fit portion only."""
    best, best_err = gammas[0], float("inf")
    for g in gammas:
        iv = c3(cal_y, cal_pred, val_y, val_pred, alpha, g)
        err = abs(coverage(val_y, iv) - (1 - alpha))
        if err < best_err:
            best, best_err = g, err
    return best


# ------------------------------------------------------------------ evaluation

def coverage(y, iv: Intervals) -> float:
    return float(np.mean((iv.lower <= y) & (y <= iv.upper)))


def summarize(y, iv: Intervals) -> dict:
    w = iv.width
    finite = np.isfinite(w)
    return dict(n=len(y), coverage=round(coverage(y, iv), 4),
                mean_width=round(float(w[finite].mean()), 2) if finite.any() else float("inf"),
                median_width=round(float(np.median(w[finite])), 2) if finite.any() else float("inf"),
                useless_fraction=round(float(np.mean(~finite | (w > USELESS_WIDTH))), 4), note=iv.note)


def onset_detection(y, origin, lower, threshold: float, minutes_per_target: float = 5.0,
                    sensor_floor: float = 40.0) -> dict:
    """Onset events: y < threshold with origin >= threshold. At-risk set: origin >= threshold.
    Detection: lower < threshold. Specificity and false alarms are over the at-risk
    non-events; false alarms per day use the at-risk set's span at one target per
    `minutes_per_target`."""
    y, origin, lower = (np.asarray(a, float) for a in (y, origin, lower))
    at_risk = origin >= threshold
    event = at_risk & (y < threshold)
    alarm = at_risk & (lower < threshold)
    tp, fn = int((event & alarm).sum()), int((event & ~alarm).sum())
    fp, tn = int((~event & alarm & at_risk).sum()), int((~event & ~alarm & at_risk).sum())
    days = at_risk.sum() * minutes_per_target / 1440.0
    at_floor = int((event & (y <= sensor_floor)).sum())
    return dict(threshold=threshold, n_at_risk=int(at_risk.sum()), n_onsets=tp + fn, tp=tp, fn=fn, fp=fp, tn=tn,
                sensitivity=round(tp / (tp + fn), 4) if tp + fn else float("nan"),
                specificity=round(tn / (tn + fp), 4) if tn + fp else float("nan"),
                ppv=round(tp / (tp + fp), 4) if tp + fp else float("nan"),
                false_alarms_per_day=round(fp / days, 3) if days else float("nan"),
                n_onsets_at_sensor_floor=at_floor, censored=at_floor > 0,
                n_already_below_at_origin=int(((~at_risk) & (y < threshold)).sum()))
