"""Clinical accuracy metrics, implemented and frozen before the sequence model exists.

Every function takes (y_true, y_pred) in mg/dL - reference first - and returns a dict
that carries the n it was computed on. Pure functions; the only I/O is the CSV writer.

Sources, and what was and was not accessible (D-022):

Clarke error grid. Clarke et al. (1987), Diabetes Care 10(5):622-628, Figure 1. The 1987
full text is paywalled; the zones are taken from Clarke's own restatement in Boren &
Clarke (2010), J Diabetes Sci Technol 4(1):84-97, whose Figure 1 reproduces the grid and
whose text defines the zones with their thresholds:
    A  within 20% of the reference, and/or <= 70 when the reference is <= 70
    B  deviates by more than 20% but would not lead to a clinically significant
       treatment error (everything not in another zone)
    C  "overcorrection errors": upper C above the line from (70, 180) to (290, 400),
       lower C the triangle (130, 0)-(180, 70)-(180, 0)        [geometry: Figure 1]
    D  "failure to treat": estimate inside the 70-180 target range when the reference
       is low (<= 70) or high (>= 240)
    E  "erroneous": estimate > 180 when the reference is <= 70, or <= 70 when the
       reference is high. The 2010 text says "> 240" here; Clarke's figure draws lower
       E from reference 180 and lower D from 240, and the figure is followed. That is
       the one place text and figure disagree, and it is recorded rather than resolved
       quietly.
Precedence A, E, D, C, B: a point within 20% is A whatever else it satisfies. Points on a
sloped boundary go to the less severe zone; the horizontal/vertical thresholds use the
inclusive comparisons the text states (<= 70, > 180, >= 240).

Parkes consensus error grid, type 1 diabetes. Parkes et al. (2000), Diabetes Care
23(8):1143-1148, Figure 2. The boundaries are polygon vertices; the published coordinates
are Table 1 of Pfuetzner et al. (2013), J Diabetes Sci Technol 7(5):1275-1281, with
"x axis, reference values; y axis, test device results", transcribed below exactly as
printed (x/y pairs). A point on a boundary line belongs to the less severe zone (neither
paper states a convention). Beyond 550 mg/dL the last segment is extended linearly.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ZONES = ("A", "B", "C", "D", "E")
TARGET_RANGE = (70.0, 180.0)        # mg/dL, the error grids' and consensus time-in-range band
HYPO_THRESHOLDS = (70.0, 54.0)      # level 1 and level 2 hypoglycaemia, strict "below"
SENSOR_FLOOR = 40.0                 # Enlite reporting floor; readings below are censored to it

# Pfuetzner et al. 2013, Table 1, type 1 diabetes columns, x/y as printed.
PARKES_T1 = {
    "B_lower": [(50, 0), (50, 30), (170, 145), (385, 300), (550, 450)],
    "B_upper": [(0, 50), (30, 50), (140, 170), (280, 380), (430, 550)],
    "C_lower": [(120, 0), (120, 30), (260, 130), (550, 250)],
    "C_upper": [(0, 60), (30, 60), (50, 80), (70, 110), (260, 550)],
    "D_lower": [(250, 0), (250, 40), (550, 150)],
    "D_upper": [(0, 100), (25, 100), (50, 125), (80, 215), (125, 550)],
    "E_upper": [(0, 150), (35, 155), (50, 550)],
}


def _check(y_true, y_pred):
    y, p = np.asarray(y_true, float).ravel(), np.asarray(y_pred, float).ravel()
    if y.shape != p.shape:
        raise ValueError("y_true and y_pred must have the same length")
    if not (np.isfinite(y).all() and np.isfinite(p).all()):
        raise ValueError("non-finite values")
    return y, p


# ------------------------------------------------------------------ Clarke

def clarke_zones(y_true, y_pred) -> np.ndarray:
    """Zone letter per point."""
    ref, est = _check(y_true, y_pred)
    z = np.full(len(ref), "B", dtype="<U1")
    a = (np.abs(est - ref) <= 0.2 * ref) | ((ref <= 70) & (est <= 70))
    e = ((ref <= 70) & (est > 180)) | ((ref >= 180) & (est <= 70))
    d = ((ref <= 70) | (ref >= 240)) & (est > 70) & (est <= 180)
    c = ((ref > 70) & (est > ref + 110)) | ((ref > 130) & (ref < 180) & (est < 1.4 * (ref - 130)))
    z[c] = "C"
    z[d] = "D"
    z[e] = "E"
    z[a] = "A"
    return z


# ------------------------------------------------------------------ Parkes

def _lower_fn(verts):
    """Boundary below the identity line as y = f(x): the first segment is vertical at
    x0, so a point is beyond the line iff x > x0 and y < f(x)."""
    (x0, _), (x1, y1), *rest = verts
    xs = np.array([x1] + [v[0] for v in rest], float)
    ys = np.array([y1] + [v[1] for v in rest], float)

    def beyond(x, y):
        slope_end = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        f = np.where(x <= xs[-1], np.interp(x, xs, ys), ys[-1] + slope_end * (x - xs[-1]))
        return (x > x0) & (y < f)
    return beyond


def _upper_fn(verts):
    """Boundary above the identity line as x = g(y): beyond iff y > y0 and x < g(y),
    where y0 is the y of the first vertex (the first segment is horizontal or nearly so)."""
    (_, y0), *rest = verts
    pts = [v for v in verts if v[1] > y0] if verts[1][1] == y0 else verts
    # if the first segment is horizontal, g(y) starts at its right-hand end
    if verts[1][1] == y0:
        pts = [verts[1]] + [v for v in verts[2:]]
    xs = np.array([v[0] for v in pts], float)
    ys = np.array([v[1] for v in pts], float)

    def beyond(x, y):
        slope_end = (xs[-1] - xs[-2]) / (ys[-1] - ys[-2])
        g = np.where(y <= ys[-1], np.interp(y, ys, xs), xs[-1] + slope_end * (y - ys[-1]))
        return (y > y0) & (x < g)
    return beyond


_PARKES = {k: (_lower_fn if k.endswith("lower") else _upper_fn)(v) for k, v in PARKES_T1.items()}


def parkes_zones(y_true, y_pred) -> np.ndarray:
    """Zone letter per point, type 1 grid."""
    ref, est = _check(y_true, y_pred)
    z = np.full(len(ref), "A", dtype="<U1")
    for zone in ("B", "C", "D", "E"):            # worst zone wins
        for side in ("lower", "upper"):
            fn = _PARKES.get(f"{zone}_{side}")
            if fn is not None:
                z[fn(ref, est)] = zone
    return z


def zone_table(zones: np.ndarray) -> dict:
    n = len(zones)
    out = {"n": n}
    for zl in ZONES:
        k = int((zones == zl).sum())
        out[f"zone_{zl}_n"] = k
        out[f"zone_{zl}_pct"] = round(100.0 * k / n, 3) if n else float("nan")
    out["zone_AB_pct"] = round(out["zone_A_pct"] + out["zone_B_pct"], 3) if n else float("nan")
    return out


def clarke(y_true, y_pred) -> dict:
    return dict(grid="clarke", **zone_table(clarke_zones(y_true, y_pred)))


def parkes(y_true, y_pred) -> dict:
    return dict(grid="parkes_t1", **zone_table(parkes_zones(y_true, y_pred)))


# ------------------------------------------------------------------ the rest

def mard(y_true, y_pred) -> dict:
    """Mean absolute relative difference, percent, denominator = reference value."""
    ref, est = _check(y_true, y_pred)
    return dict(n=len(ref), mard_pct=round(float(100.0 * np.mean(np.abs(est - ref) / ref)), 3),
                denominator="reference")


def tir_agreement(y_true, y_pred, lo: float = TARGET_RANGE[0], hi: float = TARGET_RANGE[1]) -> dict:
    """Do forecast and reference agree on in-range (lo <= g <= hi) vs out-of-range?"""
    ref, est = _check(y_true, y_pred)
    r_in, e_in = (ref >= lo) & (ref <= hi), (est >= lo) & (est <= hi)
    n = len(ref)
    both_in, both_out = int((r_in & e_in).sum()), int((~r_in & ~e_in).sum())
    return dict(n=n, range_lo=lo, range_hi=hi,
                agreement=round((both_in + both_out) / n, 4) if n else float("nan"),
                both_in_range=both_in, both_out_of_range=both_out,
                ref_in_pred_out=int((r_in & ~e_in).sum()), ref_out_pred_in=int((~r_in & e_in).sum()),
                reference_in_range_pct=round(100.0 * r_in.mean(), 3) if n else float("nan"))


def hypo_detection(y_true, y_pred, threshold: float) -> dict:
    """Sensitivity, specificity and PPV of 'forecast below threshold' against 'reference
    below threshold'. The censoring caveat is part of the result, not only prose: readings
    below the sensor floor are reported as the floor, so events at the floor are of unknown
    true depth, and the reference cannot contain any value below it."""
    ref, est = _check(y_true, y_pred)
    ev, det = ref < threshold, est < threshold
    tp, fn = int((ev & det).sum()), int((ev & ~det).sum())
    fp, tn = int((~ev & det).sum()), int((~ev & ~det).sum())
    at_floor = int((ev & (ref <= SENSOR_FLOOR)).sum())
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    ppv = tp / (tp + fp) if tp + fp else float("nan")
    censored = at_floor > 0
    return dict(threshold=threshold, n=len(ref), n_events=tp + fn, tp=tp, fn=fn, fp=fp, tn=tn,
                sensitivity=round(sens, 4), specificity=round(spec, 4), ppv=round(ppv, 4),
                n_events_at_sensor_floor=at_floor,
                share_of_events_at_floor=round(at_floor / (tp + fn), 4) if tp + fn else float("nan"),
                censored=censored,
                caveat=("" if not censored else
                        f"{at_floor} of {tp + fn} reference events sit at the {SENSOR_FLOOR:.0f} mg/dL "
                        f"sensor floor: their true depth is unknown and errors below the floor are "
                        f"unmeasurable, so this figure is biased optimistic"))


def score(y_true, y_pred) -> dict:
    """Every metric above, flattened, with one n."""
    ref, est = _check(y_true, y_pred)
    out = {"n": len(ref)}
    for grid, fn in (("clarke", clarke), ("parkes", parkes)):
        r = fn(ref, est)
        out.update({f"{grid}_{k}": v for k, v in r.items() if k not in ("n", "grid")})
    out["mard_pct"] = mard(ref, est)["mard_pct"]
    t = tir_agreement(ref, est)
    out.update({f"tir_{k}": v for k, v in t.items() if k not in ("n", "range_lo", "range_hi")})
    for thr in HYPO_THRESHOLDS:
        h = hypo_detection(ref, est, thr)
        out.update({f"hypo{thr:.0f}_{k}": v for k, v in h.items() if k not in ("n", "threshold")})
    return out


def score_groups(groups: list[dict]) -> list[dict]:
    """groups: dicts with patient, cohort, horizon_min, method, y_true, y_pred. Returns rows
    per patient, per cohort and pooled, every one computed on the concatenated points of
    its members (stated as `level`) and carrying its n."""
    rows = []
    keys = sorted({(g["method"], g["horizon_min"]) for g in groups})
    for method, h in keys:
        members = [g for g in groups if g["method"] == method and g["horizon_min"] == h]
        for g in members:
            rows.append(dict(method=method, horizon_min=h, level="patient", patient=g["patient"],
                             cohort=g["cohort"], n_patients=1, **score(g["y_true"], g["y_pred"])))
        for cohort in sorted({g["cohort"] for g in members}) + ["pooled"]:
            sub = [g for g in members if cohort == "pooled" or g["cohort"] == cohort]
            y = np.concatenate([g["y_true"] for g in sub])
            p = np.concatenate([g["y_pred"] for g in sub])
            rows.append(dict(method=method, horizon_min=h, level="cohort" if cohort != "pooled" else "pooled",
                             patient="", cohort=cohort, n_patients=len(sub), **score(y, p)))
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
