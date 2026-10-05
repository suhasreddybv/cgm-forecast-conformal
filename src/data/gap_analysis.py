"""B3: how much CGM is missing, and where windowing would cross a gap.

Gaps are to glucose forecasting what the window-label alignment check was to the PPG
repository: handled wrong, they corrupt every result silently. A history window that
spans a two-hour sensor outage is not a history window, and a target interpolated
across a gap is not a measurement.

Aggregate only: counts and fractions per patient. No timestamps, no traces, nothing
that could reconstruct a patient's glucose curve (data use agreement).

Writes results/gap_analysis.csv.
"""
from __future__ import annotations

import csv

import numpy as np

from src.data.loader import CGM_INTERVAL_S, PATIENTS, REPO_ROOT, load_patient

BUCKETS = [(0, 10), (10, 30), (30, 60), (60, 360), (360, 10 ** 9)]
BUCKET_NAMES = ["<=10min", "<=30min", "<=1h", "<=6h", ">6h"]
HISTORY_LENGTHS = (6, 12, 24)


def main() -> None:
    rows = []
    print(f"{'pid':<5}{'split':<7}{'n':>7}{'expected':>10}{'missing':>9}"
          + "".join(f"{b:>9}" for b in BUCKET_NAMES) + f"{'longest':>9}"
          + "".join(f"{'H=' + str(h):>8}" for h in HISTORY_LENGTHS))
    for pid in PATIENTS:
        for split in ("train", "test"):
            rec = load_patient(pid, split)
            ts = rec.cgm.ts
            d = np.diff(ts).astype("timedelta64[s]").astype(int)
            span_s = int((ts[-1] - ts[0]).astype("timedelta64[s]").astype(int))
            expected = span_s // CGM_INTERVAL_S + 1
            missing = 1.0 - len(ts) / expected

            # a "gap" is any interval longer than the nominal cadence; its missing
            # readings are the ones that should have fallen inside it
            gap_mask = d > CGM_INTERVAL_S + 60
            gap_min = d[gap_mask] / 60.0
            counts = [int(((gap_min > lo) & (gap_min <= hi)).sum()) for lo, hi in BUCKETS]

            # how many history windows of length H would span a gap if cut naively
            crossing = {}
            for H in HISTORY_LENGTHS:
                if len(d) < H:
                    crossing[H] = float("nan")
                    continue
                # window i covers readings [i, i+H); it crosses a gap if any interval inside does
                strided = np.lib.stride_tricks.sliding_window_view(gap_mask, H - 1)
                crossing[H] = float(strided.any(axis=1).mean())

            longest = float(gap_min.max()) if gap_mask.any() else 0.0
            print(f"{pid:<5}{split:<7}{len(ts):>7}{expected:>10}{missing:>8.1%}"
                  + "".join(f"{c:>9}" for c in counts) + f"{longest:>8.0f}m"
                  + "".join(f"{crossing[h]:>8.1%}" for h in HISTORY_LENGTHS))
            rows.append(dict(
                patient=pid, cohort=rec.cohort, split=split, readings=len(ts),
                expected_readings=expected, missing_fraction=round(missing, 4),
                gaps_le_10min=counts[0], gaps_le_30min=counts[1], gaps_le_1h=counts[2],
                gaps_le_6h=counts[3], gaps_over_6h=counts[4],
                longest_gap_min=round(longest, 1), total_gaps=int(gap_mask.sum()),
                short_intervals=rec.short_intervals,
                **{f"windows_crossing_H{h}": round(crossing[h], 4) for h in HISTORY_LENGTHS}))

    out = REPO_ROOT / "results" / "gap_analysis.csv"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote results/{out.name}")

    tr = [r for r in rows if r["split"] == "train"]
    te = [r for r in rows if r["split"] == "test"]
    for name, grp in (("train", tr), ("test", te)):
        miss = np.mean([r["missing_fraction"] for r in grp])
        print(f"  {name}: mean missing {miss:.1%}, "
              f"worst {max(r['missing_fraction'] for r in grp):.1%} "
              f"({max(grp, key=lambda r: r['missing_fraction'])['patient']}), "
              f"longest gap {max(r['longest_gap_min'] for r in grp) / 60:.1f} h")
    for h in HISTORY_LENGTHS:
        v = [r[f"windows_crossing_H{h}"] for r in rows]
        print(f"  H={h:>2} ({h * 5:>3} min of history): {np.mean(v):.1%} of naive windows cross a gap "
              f"(worst patient {max(v):.1%})")


if __name__ == "__main__":
    main()
