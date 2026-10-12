"""Sunday addendum §1: where did every fit in the staircase stop, and did the networks
learn beyond their persistence initialisation?

Every model starts exactly at persistence (zero-initialised change head, D-026), so a fit
whose best validation epoch is 1 has learned one epoch's worth and no more. This lists,
for every fit in Steps 1-7, the best epoch and the stop epoch, and draws the histogram of
best epochs per step with patients coloured, so an early-stopping pathology confined to
one patient shows as one colour.

Writes results/best_epochs.csv and figures/best_epochs.png. Aggregate only.
"""
from __future__ import annotations

import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.data.loader import COHORT_2018, COHORT_2020, PATIENTS, REPO_ROOT  # noqa: E402
from src.staircase import RUNS, RUNS_SNAPSHOT  # noqa: E402

STEP_ORDER = ["step1", "step2", "step3", "step4", "step4ens", "step5pre", "step5ft", "step5pop",
              "step6", "step6fh", "step6seeds", "step6ens", "step7pre", "step7", "step7ens", "step7prenocov", "step7nocov"]
# twelve patients, a fixed categorical order (2018 cohort then 2020), never cycled
COLOURS = dict(zip(COHORT_2018 + COHORT_2020,
                   ["#1f77b4", "#aec7e8", "#17becf", "#9edae5", "#2ca02c", "#98df8a",
                    "#d62728", "#ff9896", "#ff7f0e", "#ffbb78", "#9467bd", "#c5b0d5"]))


def load_rows() -> list[dict]:
    src = RUNS if RUNS.exists() else RUNS_SNAPSHOT
    return [json.loads(l) for l in src.read_text().splitlines() if l.strip()]


def main() -> None:
    rows = load_rows()
    out = [dict(step=r["step"], fit=r["fit"], patient=r["patient"] or "all12", horizon_min=r["horizon_steps"] * 5,
                seed=r["seed"], best_epoch=r["best_epoch"], stop_epoch=r["epochs_run"],
                learned_beyond_init=r["best_epoch"] > 1, val_rmse=r["val_rmse"]) for r in rows]
    with open(REPO_ROOT / "results" / "best_epochs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print(f"wrote results/best_epochs.csv ({len(out)} fits)")

    steps = [s for s in STEP_ORDER if any(r["step"] == s for r in out)]
    print(f"\n{'step':<14}{'fits':>6}{'best=1':>8}{'best<=3':>9}{'median best':>12}{'median stop':>12}  best=1 by patient")
    for s in steps:
        g = [r for r in out if r["step"] == s]
        be = np.array([r["best_epoch"] for r in g]); st = np.array([r["stop_epoch"] for r in g])
        early = {}
        for r in g:
            if r["best_epoch"] == 1:
                early[r["patient"]] = early.get(r["patient"], 0) + 1
        print(f"{s:<14}{len(g):>6}{(be == 1).sum():>8}{(be <= 3).sum():>9}{np.median(be):>12.0f}{np.median(st):>12.0f}  "
              + ", ".join(f"{p}×{n}" for p, n in sorted(early.items(), key=lambda x: -x[1])))

    n = len(steps)
    cols = 3
    rows_ = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows_, cols, figsize=(4.6 * cols, 3.2 * rows_), squeeze=False)
    bins = np.arange(0.5, 101.5, 5)
    for ax, s in zip(axes.ravel(), steps):
        g = [r for r in out if r["step"] == s]
        data = [[min(r["best_epoch"], 100) for r in g if r["patient"] == p] for p in PATIENTS]
        pop = [min(r["best_epoch"], 100) for r in g if r["patient"] == "all12"]
        ax.hist(data, bins=bins, stacked=True, color=[COLOURS[p] for p in PATIENTS], label=PATIENTS,
                edgecolor="white", linewidth=0.4)
        if pop:
            ax.hist(pop, bins=bins, color="#555555", label="population", edgecolor="white", linewidth=0.4)
        ax.set_title(f"{s}  (n = {len(g)}; best epoch 1: {sum(r['best_epoch'] == 1 for r in g)})", fontsize=9)
        ax.set_xlabel("best validation epoch (capped at 100)", fontsize=8)
        ax.tick_params(labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=6, ncol=2, frameon=False, title="patient", title_fontsize=7)
    fig.suptitle("Where each fit stopped: best validation epoch per step, patients coloured. Epoch 1 = one epoch past persistence.",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(REPO_ROOT / "figures" / "best_epochs.png", dpi=130)
    print("wrote figures/best_epochs.png")


if __name__ == "__main__":
    main()
