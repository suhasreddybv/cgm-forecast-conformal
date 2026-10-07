"""Error metrics on real targets (mg/dL). Aggregation is the caller's business and is stated."""
from __future__ import annotations

import numpy as np


def rmse(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.sqrt(np.mean((p - y) ** 2)))


def mae(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean(np.abs(p - y)))


def mape(y: np.ndarray, p: np.ndarray) -> float:
    """Percent. Targets are censored to [40, 400], so the denominator is never small."""
    return float(100.0 * np.mean(np.abs(p - y) / np.abs(y)))
