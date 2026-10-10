"""Part E: the frozen clinical scorer applied to the final model beside persistence and the
linear AR - conforming variant, 2020 cohort - and the two error-grid figures that replace
the README's placeholder line.

    python -m src.eval.clinical_model [--step step7ens]

Writes results/clinical_final.csv (per patient, per cohort, pooled; every row with its n and
the censoring flag), figures/hero_clarke.png and figures/parkes_final.png. Predictions are
the seed-mean of the step's cached predictions (median over H for an ensemble step).
Aggregate only.
"""
from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.data.loader import COHORT_2020, COHORT_OF, REPO_ROOT  # noqa: E402
from src.eval.clinical import clarke_zones, parkes_zones, score_groups, write_csv  # noqa: E402
from src.eval.clinical_baselines import _draw_clarke, _draw_parkes  # noqa: E402
from src.eval.staircase_eval import STEP_LABELS, baseline_preds, per_patient_seed  # noqa: E402

VARIANT = "bglp"


def groups_for(step: str) -> list[dict]:
    ens = next(e for s, _, e, _ in STEP_LABELS if s == step)
    pp = per_patient_seed(step, VARIANT, ens)
    groups = []
    for h in (6, 12):
        for pid in COHORT_2020:
            seeds = [k[2] for k in pp if k[0] == pid and k[1] == h]
            ts, y, _ = pp[(pid, h, seeds[0])]
            p_model = np.mean([pp[(pid, h, s)][2] for s in seeds], axis=0)
            bts, by, bp0, bl2 = baseline_preds(pid, h, VARIANT)
            m = np.isin(bts, ts)
            assert np.array_equal(bts[m], ts)
            for name, pred in (("persistence", bp0[m]), ("linear AR H=24 +cov", bl2[m]), (f"model ({step})", p_model)):
                groups.append(dict(patient=pid, cohort=COHORT_OF[pid], horizon_min=h * 5, method=name, y_true=y, y_pred=pred))
    return groups


def figures(groups: list[dict], step: str) -> None:
    for grid, draw, fn, lim, out in (("Clarke", _draw_clarke, clarke_zones, 400, "hero_clarke.png"),
                                     ("Parkes (type 1)", _draw_parkes, parkes_zones, 550, "parkes_final.png")):
        methods = ["persistence", f"model ({step})"]
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.2), sharex=True, sharey=True)
        for ax, method in zip(axes, methods):
            sub = [g for g in groups if g["method"] == method and g["horizon_min"] == 30]
            y = np.concatenate([g["y_true"] for g in sub]); p = np.concatenate([g["y_pred"] for g in sub])
            z = fn(y, p)
            ax.hexbin(y, p, gridsize=55, cmap="Blues", bins="log", mincnt=1, extent=(0, lim, 0, lim))
            draw(ax)
            pct = "   ".join(f"{zl} {100 * np.mean(z == zl):.1f}%" for zl in "ABCDE")
            ax.set_title(f"{method}\nn = {len(y):,}   {pct}", fontsize=10)
            ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_xlabel("reference CGM (mg/dL)")
        axes[0].set_ylabel("30-minute forecast (mg/dL)")
        fig.suptitle(f"{grid} error grid, 30 min, conforming evaluation, 2020 cohort pooled — persistence against the final model",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(REPO_ROOT / "figures" / out, dpi=130)
        print(f"wrote figures/{out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", default="step7ens")
    a = ap.parse_args()
    groups = groups_for(a.step)
    rows = score_groups(groups)
    write_csv(rows, REPO_ROOT / "results" / "clinical_final.csv")
    print("wrote results/clinical_final.csv")
    figures(groups, a.step)
    print(f"\n{'method':<24}{'h':>4}{'n':>7}{'Clarke A':>9}{'A+B':>7}{'Parkes A':>9}{'MARD':>7}{'TIR':>7}{'sens70':>8}{'sens54':>8}{'floor':>6}")
    for r in rows:
        if r["level"] == "pooled":
            print(f"{r['method']:<24}{r['horizon_min']:>4}{r['n']:>7}{r['clarke_zone_A_pct']:>9.1f}{r['clarke_zone_AB_pct']:>7.1f}"
                  f"{r['parkes_zone_A_pct']:>9.1f}{r['mard_pct']:>7.1f}{r['tir_agreement']:>7.3f}{r['hypo70_sensitivity']:>8.3f}"
                  f"{r['hypo54_sensitivity']:>8.3f}{r['hypo54_n_events_at_sensor_floor']:>6}")
    print("\n60-min Clarke zone D per patient (Zhu's worst case: 8%) and hypo-70 sensitivity:")
    for r in rows:
        if r["level"] == "patient" and r["horizon_min"] == 60:
            print(f"  {r['method'][:20]:<21}{r['patient']}  D {r['clarke_zone_D_pct']:>5.1f}%  sens70 {r['hypo70_sensitivity']:.3f} ({r['hypo70_n_events']} events)")


if __name__ == "__main__":
    main()
