"""Forecast ANALOGUES of the iCGM special-controls accuracy criteria, 21 CFR 862.1355(b).

Every threshold is parsed from docs/icgm_criteria.md at run time - the regulation as
transcribed by Suhas on 15 Oct 2026 - and none is typed here (D-059). The criteria govern a
sensor against a laboratory reference; this repository compares a forecast against a CGM,
so everything computed here is an analogue: substitute "forecast" (or an interval bound)
for "iCGM measurement" and "CGM reading at the target time" for "reference blood glucose".
A pass shows agreement with a sensor, not accuracy against blood; a fail fails something
easier than the regulation asks. Nothing here is a claim that any device meets §862.1355.

  (A)-(G)  bin by the ESTIMATE's value; share within tolerance of the target; the lower
           one-sided 95% Wilson bound is what the regulation requires - reported with the
           plain share, and with the note that serially correlated targets make it optimistic
  (H), (I) zero-tolerance exclusions: counts and rate per 1,000 estimates
  (J), (K) rate-of-change: implied forecast rate = (estimate - origin) / horizon minutes,
           true rate = (target - origin) / horizon minutes; wrong-direction share, limit 1%
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import norm

from src.data.loader import REPO_ROOT

DOC = REPO_ROOT / "docs" / "icgm_criteria.md"
MIN_N = 100     # D-053 rule 1 applies to analogue cells too


@dataclass(frozen=True)
class Accuracy:
    para: str
    lo: float          # estimate range [lo, hi): -inf / inf for open ends; whole range = (-inf, inf)
    hi: float
    lo_inclusive: bool
    hi_inclusive: bool
    tol_value: float
    tol_relative: bool  # True: within tol_value percent; False: within tol_value mg/dL
    required_pct: float


@dataclass(frozen=True)
class Criteria:
    accuracy: tuple
    exclusion_low: float       # (H): estimate < low with reference > high
    exclusion_high: float      # (I): estimate > high with reference < low
    rate_limit_pct: float      # (J), (K): at most this percent
    rate_indicated: float      # |indicated rate| threshold, mg/dL/min
    rate_true: float           # |true rate| threshold, mg/dL/min


def _parse_range(text: str):
    t = text.replace("–", "-").replace("mg/dL", "").strip()
    if "whole" in t:
        return -np.inf, np.inf, True, True
    m = re.match(r"^<\s*([\d.]+)$", t)
    if m:
        return -np.inf, float(m.group(1)), True, False
    m = re.match(r"^>\s*([\d.]+)$", t)
    if m:
        return float(m.group(1)), np.inf, False, True
    m = re.match(r"^([\d.]+)\s*-\s*([\d.]+)$", t)
    if m:
        return float(m.group(1)), float(m.group(2)), True, True
    raise ValueError(f"unparsed range {text!r}")


def load_criteria(doc: Path = DOC) -> Criteria:
    text = doc.read_text()
    acc = []
    for para, rng, tol, req in re.findall(r"^\|\s*\((\w)\)\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*([\d.]+)%\s*\|$", text, flags=re.M):
        lo, hi, loi, hii = _parse_range(rng)
        m = re.search(r"±\s*([\d.]+)\s*(%|mg/dL)", tol)
        if not m:
            raise ValueError(f"unparsed tolerance {tol!r}")
        acc.append(Accuracy(para, lo, hi, loi, hii, float(m.group(1)), m.group(2) == "%", float(req)))
    h = re.search(r"\(H\)\s*\|\s*When the iCGM reads <\s*([\d.]+), no reference value may be >\s*([\d.]+)", text)
    i = re.search(r"\(I\)\s*\|\s*When the iCGM reads >\s*([\d.]+), no reference value may be <\s*([\d.]+)", text)
    j = re.search(r"\(J\)\s*\|\s*No more than\s*([\d.]+)% .*? rate of change >\s*\+?([\d.]+) mg/dL/min when the true rate .*? <\s*[−-]([\d.]+) mg/dL/min", text)
    if not (h and i and j and len(acc) == 7):
        raise ValueError("docs/icgm_criteria.md does not contain the criteria in the expected form")
    assert float(h.group(1)) == float(i.group(2)) and float(h.group(2)) == float(i.group(1))
    return Criteria(tuple(acc), float(h.group(1)), float(h.group(2)), float(j.group(1)), float(j.group(2)), float(j.group(3)))


def wilson_lower_one_sided(k: int, n: int, conf: float = 0.95) -> float:
    """Lower one-sided Wilson bound on a proportion, in percent."""
    if n == 0:
        return float("nan")
    z = norm.ppf(conf)
    p = k / n
    centre = p + z * z / (2 * n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return 100.0 * (centre - half) / (1 + z * z / n)


def accuracy_analogue(reference, estimate, crit: Accuracy) -> dict:
    """Bin by the ESTIMATE (the iCGM reading in the regulation); share of estimates within
    tolerance of the reference (the CGM target here)."""
    ref, est = np.asarray(reference, float), np.asarray(estimate, float)
    fin = np.isfinite(est)
    inb = fin & ((est >= crit.lo) if crit.lo_inclusive else (est > crit.lo)) & ((est <= crit.hi) if crit.hi_inclusive else (est < crit.hi))
    n = int(inb.sum())
    if crit.tol_relative:
        within = np.abs(est[inb] - ref[inb]) <= crit.tol_value / 100.0 * ref[inb]
    else:
        within = np.abs(est[inb] - ref[inb]) <= crit.tol_value
    k = int(within.sum())
    pct = 100.0 * k / n if n else float("nan")
    lower = wilson_lower_one_sided(k, n)
    return dict(para=crit.para, n=n, within=k, pct=round(pct, 2), wilson_lower_95=round(lower, 2) if n else float("nan"),
                required_pct=crit.required_pct, met=bool(lower > crit.required_pct) if n else None, reportable=bool(n >= MIN_N),
                n_non_finite_estimates=int((~fin).sum()))


def exclusion_analogue(reference, estimate, crit: Criteria) -> dict:
    ref, est = np.asarray(reference, float), np.asarray(estimate, float)
    fin = np.isfinite(est)
    h = int((fin & (est < crit.exclusion_low) & (ref > crit.exclusion_high)).sum())
    i = int((fin & (est > crit.exclusion_high) & (ref < crit.exclusion_low)).sum())
    n = int(fin.sum())
    return dict(H_violations=h, H_per_1000=round(1000 * h / n, 3) if n else float("nan"),
                I_violations=i, I_per_1000=round(1000 * i / n, 3) if n else float("nan"), n=n, requirement="zero")


def rate_analogue(reference, estimate, origin, horizon_min: float, crit: Criteria) -> dict:
    ref, est, org = (np.asarray(a, float) for a in (reference, estimate, origin))
    fin = np.isfinite(est)
    indicated = (est - org) / horizon_min
    true = (ref - org) / horizon_min
    j = int((fin & (indicated > crit.rate_indicated) & (true < -crit.rate_true)).sum())
    k = int((fin & (indicated < -crit.rate_indicated) & (true > crit.rate_true)).sum())
    n = int(fin.sum())
    return dict(J_wrong_direction_pct=round(100 * j / n, 3) if n else float("nan"), J_count=j,
                K_wrong_direction_pct=round(100 * k / n, 3) if n else float("nan"), K_count=k, n=n,
                limit_pct=crit.rate_limit_pct, J_met=bool(100 * j / n <= crit.rate_limit_pct) if n else None,
                K_met=bool(100 * k / n <= crit.rate_limit_pct) if n else None)
