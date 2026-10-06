"""The effective-n table, and the reconciliation of the windower against the gap analysis.

results/effective_n.csv  per patient x split x horizon x H: how many real targets are
    evaluable under the frozen protocol, as a count and as a fraction of the raw readings.
    This is the n column that accompanies every result from here on (D-009).

results/window_reconciliation.csv  the windower's gap accounting against the naive
    index-space count of results/gap_analysis.csv, recomputed here over the same candidate
    targets. The two measure different things, and the difference is located exactly:
      naive    the H consecutive *readings* ending h readings before the target contain a
               gap interval. In index space a window never notices a gap on its horizon
               side - the "latest" reading may be hours old.
      windower the H time *slots* ending exactly h steps before the target touch a gap.
    Splitting the windower's count by whether the history END slot itself falls inside a
    gap separates the two: targets touched only on the history side must match the naive
    count, and targets whose history end sits inside a gap are what the naive count misses.
    Both are reported per file; the unexplained remainder is printed and should be small
    (series-start precedence and cadence jitter).

Aggregate counts only - no timestamps, no values (D-001).
"""
from __future__ import annotations

import csv

import numpy as np

from src.data.loader import PATIENTS, REPO_ROOT, load_patient
from src.data.windows import (DROPPED_NO_HISTORY, GAP_S, HISTORY_LENGTHS, HORIZONS, STEP_S,
                              _series_for, build_windows, windows_for)


def naive_crossing(ts: np.ndarray, history_len: int, horizon: int, eval_from: int):
    """Index-space count, as gap_analysis.py does it, over the same candidate targets:
    do the H readings ending h readings before target j contain a gap interval?
    Returns (count, number of candidates with enough readings)."""
    t = ts.astype("datetime64[s]").astype(np.int64)
    gap = np.diff(t) > GAP_S                              # gap[i]: between reading i and i+1
    cand = np.arange(eval_from, len(t))
    first = cand - horizon - history_len + 1               # oldest reading of the index-space history
    ok = first >= 0
    cs = np.concatenate([[0], np.cumsum(gap)])
    crossing = np.zeros(len(cand), bool)
    crossing[ok] = (cs[cand[ok] - horizon] - cs[first[ok]]) > 0
    return int(crossing.sum()), int(ok.sum())


def main() -> None:
    rows, rec_rows = [], []
    print(f"{'pid':<5}{'split':<7}{'h':>4}{'H':>4}{'readings':>9}{'cand':>7}{'evaluable':>10}"
          f"{'frac':>7}{'interp':>8}{'long':>7}{'nohist':>8}  |{'naive':>7}{'touched':>9}{'hist':>7}{'end':>7}")
    for pid in PATIENTS:
        for split in ("train", "test"):
            rec = load_patient(pid, split)
            cgm, _, eval_from = _series_for(rec)
            # readings the sensor should have produced over the split's span (as gap_analysis.py)
            span_s = int((rec.cgm.ts[-1] - rec.cgm.ts[0]).astype("timedelta64[s]").astype(int))
            expected = span_s // STEP_S + 1
            for h in HORIZONS:
                for H in HISTORY_LENGTHS:
                    ws = windows_for(pid, split, h, H)
                    n_interp = int(ws.interpolated.sum())
                    rows.append(dict(
                        patient=pid, cohort=rec.cohort, split=split, horizon_min=h * 5,
                        horizon_steps=h, history_steps=H, readings=ws.n_readings,
                        expected_readings=expected,
                        candidate_targets=ws.n_candidates, evaluable_targets=len(ws),
                        fraction_of_readings=round(len(ws) / ws.n_readings, 4),
                        fraction_of_expected=round(len(ws) / expected, 4),
                        fraction_of_candidates=round(len(ws) / ws.n_candidates, 4),
                        evaluable_interpolated_history=n_interp,
                        evaluable_clean_history=len(ws) - n_interp,
                        dropped_gap_over_30min=ws.n_dropped_long_gap,
                        dropped_no_history=ws.n_dropped_no_history))
                    naive, n_idx = naive_crossing(cgm.ts, H, h, eval_from)
                    b = build_windows(cgm.ts, cgm.values, h, H, eval_from)
                    has_hist = b.status != DROPPED_NO_HISTORY
                    touched = has_hist & (b.max_bracket_s > GAP_S)
                    end_in_gap = has_hist & (b.end_bracket_s > GAP_S)
                    history_side = int((touched & ~end_in_gap).sum())
                    rec_rows.append(dict(
                        patient=pid, split=split, horizon_steps=h, history_steps=H,
                        naive_index_windows_crossing=naive, naive_fraction=round(naive / n_idx, 4),
                        targets_history_touches_gap=int(touched.sum()),
                        touched_fraction=round(touched.sum() / has_hist.sum(), 4),
                        of_which_dropped=ws.n_dropped_long_gap, of_which_interpolated=n_interp,
                        touched_on_history_side_only=history_side,
                        touched_with_history_end_inside_gap=int(end_in_gap.sum()),
                        history_side_minus_naive=history_side - naive))
                    print(f"{pid:<5}{split:<7}{h:>4}{H:>4}{ws.n_readings:>9}{ws.n_candidates:>7}"
                          f"{len(ws):>10}{len(ws) / ws.n_readings:>7.1%}{n_interp:>8}"
                          f"{ws.n_dropped_long_gap:>7}{ws.n_dropped_no_history:>8}  |"
                          f"{naive:>7}{int(touched.sum()):>9}{history_side:>7}{int(end_in_gap.sum()):>7}")

    res = REPO_ROOT / "results"
    for data, fn in ((rows, "effective_n.csv"), (rec_rows, "window_reconciliation.csv")):
        with open(res / fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
        print(f"wrote results/{fn}")

    # summaries the log quotes
    print("\nthe cost of missingness is paid in absolute n: evaluable targets against the readings"
          "\nthe sensor should have produced (test split, h=12, H=24):")
    for r in sorted((r for r in rows if r["split"] == "test" and r["horizon_steps"] == 12
                     and r["history_steps"] == 24), key=lambda r: r["fraction_of_expected"])[:3]:
        print(f"  {r['patient']}: {r['evaluable_targets']} evaluable of {r['readings']} readings "
              f"({r['fraction_of_readings']:.1%}) and of {r['expected_readings']} expected "
              f"({r['fraction_of_expected']:.1%})")
    print("\nevaluable targets as a fraction of raw readings (mean over patients):")
    for split in ("train", "test"):
        for h in HORIZONS:
            for H in HISTORY_LENGTHS:
                v = [r["fraction_of_readings"] for r in rows
                     if r["split"] == split and r["horizon_steps"] == h and r["history_steps"] == H]
                worst = min((r for r in rows if r["split"] == split and r["horizon_steps"] == h
                             and r["history_steps"] == H), key=lambda r: r["fraction_of_readings"])
                print(f"  {split:<6} h={h:>2} H={H:>2}: {np.mean(v):.1%}  "
                      f"(worst {worst['patient']} {worst['fraction_of_readings']:.1%})")
    print("\nreconciliation (all files pooled):")
    for h in HORIZONS:
        for H in HISTORY_LENGTHS:
            g = [r for r in rec_rows if r["horizon_steps"] == h and r["history_steps"] == H]
            naive = sum(r["naive_index_windows_crossing"] for r in g)
            touched = sum(r["targets_history_touches_gap"] for r in g)
            hist = sum(r["touched_on_history_side_only"] for r in g)
            end = sum(r["touched_with_history_end_inside_gap"] for r in g)
            resid = sum(abs(r["history_side_minus_naive"]) for r in g)
            print(f"  h={h:>2} H={H:>2}: naive {naive:>6}  history-side {hist:>6} "
                  f"(|residual| over files {resid:>3})  + history-end-in-gap {end:>6}  = touched {touched:>6}")

if __name__ == "__main__":
    main()
