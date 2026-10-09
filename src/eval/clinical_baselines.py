"""Clinical metrics applied to the baselines only: p0 persistence and the best linear AR
(per-patient, H=24, with covariates), at 30 and 60 minutes, test split, scored targets.

Writes results/clinical_baselines.csv (per patient, per cohort, pooled - every row with its
n) and two figures, figures/clarke_baselines.png and figures/parkes_baselines.png, labelled
as baselines. They are not embedded in the README: the hero image is the model-versus-
baseline comparison, which does not exist yet (D-024).

Figures are hexbin densities of (reference, forecast) pairs - aggregate, no trace, no
timestamp (D-001).
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.data.loader import COHORT_OF, PATIENTS, REPO_ROOT  # noqa: E402
from src.data.windows import HORIZONS, windows_for  # noqa: E402
from src.eval.baselines import covariate_ok, design, fit_ols, persistence  # noqa: E402
from src.eval.clinical import PARKES_T1, clarke_zones, parkes_zones, score_groups, write_csv  # noqa: E402

BEST_H = 24
METHODS = ("p0 persistence", "l2 linear AR H=24 per-patient +cov")


def baseline_predictions() -> list[dict]:
    groups = []
    for pid in PATIENTS:
        for h in HORIZONS:
            ws = windows_for(pid, "test", h, BEST_H)
            tr = windows_for(pid, "train", h, BEST_H)
            keep = covariate_ok(ws)
            assert keep.all(), "no test window carries a NaN covariate (D-020)"
            groups.append(dict(patient=pid, cohort=COHORT_OF[pid], horizon_min=h * 5, method=METHODS[0],
                               y_true=ws.target, y_pred=persistence(ws)))
            beta = fit_ols([tr], True)
            groups.append(dict(patient=pid, cohort=COHORT_OF[pid], horizon_min=h * 5, method=METHODS[1],
                               y_true=ws.target, y_pred=design(ws, True) @ beta))
    return groups


def _draw_clarke(ax):
    ax.plot([0, 400], [0, 400], color="black", lw=0.6)
    ax.plot([0, 400 / 1.2], [0, 400], color="black", lw=0.8, ls="--")          # +20%
    ax.plot([0, 400], [0, 320], color="black", lw=0.8, ls="--")                # -20%
    ax.plot([70, 290], [180, 400], color="black", lw=0.8)                      # upper C line
    ax.plot([130, 180], [0, 70], color="black", lw=0.8)                        # lower C
    ax.plot([180, 180], [0, 70], color="black", lw=0.8)
    for x0, x1, y0, y1 in ((0, 70, 70, 180), (0, 70, 180, 400), (240, 400, 70, 180), (180, 400, 0, 70),
                           (0, 70, 0, 70)):
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color="black", lw=0.8)
    for x, y, t in ((35, 290, "E"), (35, 125, "D"), (170, 330, "C"), (320, 125, "D"), (290, 35, "E"),
                    (165, 15, "C"), (300, 290, "A"), (150, 230, "B"), (230, 120, "B")):
        ax.text(x, y, t, ha="center", va="center", fontsize=9)


def _draw_parkes(ax):
    ax.plot([0, 550], [0, 550], color="black", lw=0.6)
    for name, verts in PARKES_T1.items():
        xs, ys = zip(*verts)
        ax.plot(xs, ys, color="black", lw=0.8)
    # label positions checked against the boundary lines at that x or y
    for x, y, t in ((500, 480, "A"), (470, 300, "B"), (470, 170, "C"), (470, 60, "D"),
                    (330, 520, "B"), (185, 520, "C"), (85, 520, "D"), (20, 500, "E")):
        ax.text(x, y, t, ha="center", va="center", fontsize=9)


def figures(groups: list[dict]) -> None:
    for grid, draw, fn, lim in (("Clarke", _draw_clarke, clarke_zones, 400), ("Parkes (type 1)", _draw_parkes, parkes_zones, 550)):
        fig, axes = plt.subplots(2, 2, figsize=(11, 10.5), sharex=True, sharey=True)
        for i, method in enumerate(METHODS):
            for j, h in enumerate((30, 60)):
                ax = axes[i, j]
                sub = [g for g in groups if g["method"] == method and g["horizon_min"] == h]
                y = np.concatenate([g["y_true"] for g in sub])
                p = np.concatenate([g["y_pred"] for g in sub])
                z = fn(y, p)
                ax.hexbin(y, p, gridsize=55, cmap="Blues", bins="log", mincnt=1, extent=(0, lim, 0, lim))
                draw(ax)
                pct = " ".join(f"{zl}:{100 * np.mean(z == zl):.1f}" for zl in "ABCDE")
                ax.set_title(f"BASELINE {method}, {h} min\nn = {len(y):,}   {pct}", fontsize=9)
                ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        for ax in axes[-1]:
            ax.set_xlabel("reference CGM (mg/dL)")
        for ax in axes[:, 0]:
            ax.set_ylabel("forecast (mg/dL)")
        fig.suptitle(f"{grid} error grid — baselines only, test split, all 12 patients pooled. Not the model.", fontsize=11)
        fig.tight_layout()
        out = REPO_ROOT / "figures" / f"{grid.split()[0].lower()}_baselines.png"
        fig.savefig(out, dpi=120)
        print(f"wrote figures/{out.name}")


def main() -> None:
    groups = baseline_predictions()
    rows = score_groups(groups)
    write_csv(rows, REPO_ROOT / "results" / "clinical_baselines.csv")
    print("wrote results/clinical_baselines.csv")
    figures(groups)

    print(f"\n{'method':<36}{'h':>4}{'level':<8}{'n':>7}{'Clarke A':>9}{'A+B':>7}{'Parkes A':>9}{'A+B':>7}"
          f"{'MARD':>7}{'TIR agr':>8}{'sens70':>8}{'sens54':>8}{'at floor':>9}")
    for r in rows:
        if r["level"] == "patient":
            continue
        print(f"{r['method']:<36}{r['horizon_min']:>4}{r['cohort']:<8}{r['n']:>7}{r['clarke_zone_A_pct']:>9.1f}"
              f"{r['clarke_zone_AB_pct']:>7.1f}{r['parkes_zone_A_pct']:>9.1f}{r['parkes_zone_AB_pct']:>7.1f}"
              f"{r['mard_pct']:>7.1f}{r['tir_agreement']:>8.3f}{r['hypo70_sensitivity']:>8.3f}"
              f"{r['hypo54_sensitivity']:>8.3f}{r['hypo54_n_events_at_sensor_floor']:>9}")
    print("\nper-patient Clarke A+B and hypo70 sensitivity at 30 min:")
    for r in rows:
        if r["level"] == "patient" and r["horizon_min"] == 30:
            print(f"  {r['method'][:14]:<15}{r['patient']}  A+B {r['clarke_zone_AB_pct']:>5.1f}  "
                  f"sens70 {r['hypo70_sensitivity']:.3f} (events {r['hypo70_n_events']:>4})  "
                  f"sens54 {r['hypo54_sensitivity']:.3f} (events {r['hypo54_n_events']:>3}, at floor {r['hypo54_n_events_at_sensor_floor']})")


if __name__ == "__main__":
    main()
