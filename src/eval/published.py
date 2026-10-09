"""Published figures on the BGLP 2020 test points, transcribed in docs/published_comparison.md
(Suhas, 9 Oct 2026, from the official results page and the system papers). Nothing here is
from memory: every number must appear verbatim in that document, and a test checks it.

Only conforming (bglp) rows for the 2020 cohort receive a published figure, and only from
a system whose inputs match the row's input class (D-012, D-029). Everything else is context.
"""
from __future__ import annotations

import csv
from pathlib import Path

from src.data.loader import REPO_ROOT

DOC = REPO_ROOT / "docs" / "published_comparison.md"

# input class -> officially ranked systems with matching inputs:
#   (system, rmse30, rmse60, mae30, mae60, qualifies, note) - "qualifies" is Table 1's verdict;
#   a note is information, not a disqualifier.
LIKE_FOR_LIKE = {
    "CGM + insulin + meals": [
        ("Zhu (p15), GRU-GAN", 18.34, 32.21, 13.37, 24.20, True, ""),
        ("Yang (p24), multi-scale LSTM", 19.05, 32.03, 13.50, 23.83, True, ""),
        ("Rubin-Falcone (p18), residual LSTM", 18.22, 31.66, 12.83, 23.60, False,
         "partial: pre-trained on the 2018 cohort and ~15M external steps"),
    ],
    "CGM only": [
        ("Hameed (p14), vanilla RNN", 19.21, 31.77, 13.08, 23.09, True, ""),
        ("Bevan (p17), population LSTM", 18.23, 31.10, 14.37, 25.75, True, "population model: the direct LOPO comparator"),
    ],
}

# published classical baselines on the same test points (Joedicke p25; Ma p27), RMSE 30 / 60
PUBLISHED_CLASSICAL = [
    ("linear regression (Joedicke p25)", 21.22, 32.88),
    ("random forest (Joedicke p25)", 22.27, 34.14),
    ("ARIMA (Joedicke p25)", 24.51, 38.66),
    ("AR (Ma p27)", 21.80, None),
    ("ARMA (Ma p27)", 21.44, None),
]

NOTE = ("three of eight officially ranked systems used inputs matching ours; two match CGM-only; "
        "three are online or wearable-based and shown for context")


def input_class(method: str, covariates) -> str | None:
    cov = covariates in (True, "True")
    if method == "l2":
        return "CGM + insulin + meals" if cov else "CGM only"
    if method in ("p0", "l1"):
        return "CGM only"
    return None


def fill(rows: list[dict]) -> list[dict]:
    """Add the published columns to challenge-score rows. Filled only where cohort is 2020,
    the row is conforming (bglp) and the inputs match; blank everywhere else."""
    for r in rows:
        r.setdefault("published_like_for_like", "")
        r.setdefault("published_rmse30", "")
        r.setdefault("published_rmse60", "")
        r.setdefault("published_inputs", "")
        cls = input_class(r["method"], r["covariates"])
        if r["cohort"] == "2020" and r["eval_variant"] == "bglp" and cls and r["method"] == "l2":
            full = [s for s in LIKE_FOR_LIKE[cls] if s[5]]
            r["published_like_for_like"] = "; ".join(f"{s[0]} {s[1]:.2f} / {s[2]:.2f}" for s in full)
            r["published_rmse30"] = min(s[1] for s in full)
            r["published_rmse60"] = min(s[2] for s in full)
            r["published_inputs"] = cls
    return rows


def main() -> None:
    path = REPO_ROOT / "results" / "baselines_challenge_scores.csv"
    rows = fill(list(csv.DictReader(open(path))))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    n = sum(1 for r in rows if r["published_rmse30"] != "")
    print(f"filled {n} rows of results/{path.name} (2020 cohort, conforming, l2 only)")


if __name__ == "__main__":
    main()
