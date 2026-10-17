<!-- Part of cgm-forecast-conformal; the front page is ../README.md -->

# Two evaluation variants, reported side by side

The frozen protocol (`primary`) drops a target whose history crosses a gap over 30 minutes and values history slots by linear interpolation. The Blood Glucose Level Prediction Challenge scores **every** test point after the first hour and permits nothing after the forecast origin. On 9 October, before any sequence model was fitted, a conforming variant (`bglp`) was added beside the primary rules — not in place of them (protocol, "Disclosed addition"; D-027): zero-order hold for every slot, every challenge point predicted (2,884 / 2,704 / 2,352 / 2,377 / 2,653 / 2,731 for the 2020 six, asserted before scoring), windows the primary rule would have dropped predicted from the held reading and counted as *fallback*, and the challenge's own scoring — per-patient RMSE and MAE, mean of six, four-score sum. Every method is reported under both from here on; the gap between them is the measured cost of the primary rule's two departures.

| Baseline | Cohort | Variant | RMSE 30 | MAE 30 | RMSE 60 | MAE 60 | Four-score sum | n (30 min) | Fallback (30 / 60) | Published, like-for-like (RMSE 30 / 60) |
|---|---|---|---|---|---|---|---|---|---|---|
| **p0** persistence | 2020 | primary | 24.16 | 17.62 | 40.27 | 29.92 | **111.98** | 15,239 | — | — |
|  | 2020 | bglp | 25.48 | 18.06 | 41.42 | 30.43 | **115.39** | 15,701 | 462 / 662 | — |
|  | 2018 | primary | 22.61 | 16.41 | 36.81 | 27.49 | **103.31** | 15,673 | — | — |
|  | 2018 | bglp | 23.06 | 16.66 | 37.16 | 27.79 | **104.67** | 15,970 | 297 / 430 | — |
| *published classical baselines, same test points (2020, conforming)* | 2020 | published | linear regression 21.22 · random forest 22.27 · ARIMA 24.51 · AR 21.80 · ARMA 21.44 | | linear regression 32.88 · random forest 34.14 · ARIMA 38.66 | | | 15,701 | — | source: Joedicke p25, Ma p27 |
| l1 extrapolation, k=6 | 2020 | primary | 28.10 | 19.53 | 57.14 | 40.45 | **145.22** | 15,239 | — | — |
|  | 2020 | bglp | 30.92 | 20.50 | 60.95 | 41.69 | **154.06** | 15,701 | 462 / 662 | — |
|  | 2018 | primary | 28.20 | 18.70 | 54.58 | 37.60 | **139.07** | 15,673 | — | — |
|  | 2018 | bglp | 29.25 | 19.24 | 55.93 | 38.30 | **142.72** | 15,970 | 297 / 430 | — |
| l2 linear AR, per-patient, H=24 (CGM only) | 2020 | primary | 19.79 | 14.51 | 34.83 | 26.53 | **95.65** | 14,476 | — | — |
|  | 2020 | bglp | 22.00 | 15.79 | 36.43 | 27.76 | **101.98** | 15,701 | 1225 / 1431 | Hameed (p14) **19.21 / 31.77**; Bevan (p17) **18.23 / 31.10** |
|  | 2018 | primary | 20.10 | 14.09 | 32.74 | 24.33 | **91.26** | 15,202 | — | — |
|  | 2018 | bglp | 20.80 | 14.74 | 33.43 | 25.13 | **94.11** | 15,970 | 768 / 898 | — |
| l2 linear AR, per-patient, H=24, +cov | 2020 | primary | 18.61 | 13.56 | 32.37 | 24.40 | **88.93** | 14,476 | — | — |
|  | 2020 | bglp | 20.77 | 14.79 | 34.12 | 25.71 | **95.41** | 15,701 | 1225 / 1431 | Zhu (p15) **18.34 / 32.21**; Yang (p24) **19.05 / 32.03** |
|  | 2018 | primary | 19.42 | 13.47 | 31.57 | 23.27 | **87.73** | 15,202 | — | — |
|  | 2018 | bglp | 20.12 | 14.07 | 32.23 | 23.98 | **90.40** | 15,970 | 768 / 898 | — |
| l2 linear AR, population LOPO, H=24 (CGM only) | 2020 | primary | 20.13 | 14.77 | 35.05 | 26.68 | **96.62** | 14,476 | — | — |
|  | 2020 | bglp | 22.24 | 15.88 | 36.54 | 27.72 | **102.38** | 15,701 | 1225 / 1431 | Hameed (p14) **19.21 / 31.77**; Bevan (p17) **18.23 / 31.10** |
|  | 2018 | primary | 20.07 | 14.13 | 33.14 | 24.81 | **92.16** | 15,202 | — | — |
|  | 2018 | bglp | 21.02 | 14.96 | 34.25 | 25.91 | **96.14** | 15,970 | 768 / 898 | — |
| l2 linear AR, population LOPO, H=24, +cov | 2020 | primary | 19.55 | 14.35 | 33.98 | 25.91 | **93.78** | 14,476 | — | — |
|  | 2020 | bglp | 21.61 | 15.42 | 35.52 | 26.99 | **99.54** | 15,701 | 1225 / 1431 | Zhu (p15) **18.34 / 32.21**; Yang (p24) **19.05 / 32.03** |
|  | 2018 | primary | 19.69 | 13.79 | 32.66 | 24.29 | **90.43** | 15,202 | — | — |
|  | 2018 | bglp | 20.59 | 14.57 | 33.67 | 25.33 | **94.16** | 15,970 | 768 / 898 | — |

Full table for every baseline, k and H: `results/baselines_challenge_scores.csv`; per-patient rows with fallback counts: `results/baselines_bglp_per_patient.csv`.

**Per patient, conforming, with the fallback count beside the penalty** (2020 cohort; linear AR per-patient H=24 with covariates; `results/baselines_bglp_per_patient.csv`). The cost lands where the gaps are: 584 and 567 carry the most fallback windows and pay four to five mg/dL; 544, with the fewest, pays half a point.

| Patient | n (challenge) | Fallback 30 / 60 (rate at 60) | Persistence RMSE 30, primary → conforming | Linear AR +cov RMSE 30, primary → conforming | Linear AR +cov RMSE 60, primary → conforming |
|---|---|---|---|---|---|
| 540 | 2,884 | 145 / 164 (5.7%) | 28.26 → 29.03 | 20.91 → **22.34** (+1.42) | 38.93 → 40.26 |
| 544 | 2,704 | 87 / 105 (3.9%) | 22.53 → 22.64 | 16.58 → **17.04** (+0.46) | 28.40 → 28.36 |
| 552 | 2,352 | 210 / 252 (10.7%) | 20.94 → 20.86 | 16.68 → **17.50** (+0.82) | 29.91 → 30.40 |
| 567 | 2,377 | 261 / 315 (13.3%) | 27.52 → 30.70 | 20.08 → **25.52** (+5.44) | 35.86 → 40.34 |
| 584 | 2,653 | 348 / 391 (14.7%) | 24.66 → 28.67 | 21.00 → **25.40** (+4.39) | 34.01 → 37.96 |
| 596 | 2,731 | 174 / 204 (7.5%) | 21.21 → 21.01 | 16.42 → **16.84** (+0.42) | 27.09 → 27.43 |

**Published column:** three of eight officially ranked systems used inputs matching ours; two match CGM-only; three are online or wearable-based and shown for context ([docs/published_comparison.md](published_comparison.md), Tables 1 and 3). It is filled only on the 2020-cohort conforming rows and only from systems whose inputs match the row — CGM + insulin + meals against Zhu and Yang, CGM-only against Hameed and Bevan — transcribed from the official results page and the system papers, never from memory (D-012, D-029). Under the challenge's own rules the best linear model here is 2.4 mg/dL behind Zhu and 1.7 behind Yang at 30 minutes; no ranked system reports persistence.

- **Conforming is harder, and more so for the 2020 cohort.** Persistence's four-score sum rises from 111.98 to 115.39 on the 2020 six and the best linear AR's from 88.93 to 95.41: about 2 mg/dL of RMSE at 30 minutes for the fitted model, 1.3 for persistence. The patients with the most fallback windows pay the most — 584 and 567 lose 4–5 mg/dL of 30-minute RMSE on the linear AR; 563, with 29 fallbacks, loses nothing.
- **What the primary rule had been buying, and one thing it should not have.** Dropping post-outage targets is a disclosed choice. The primary interpolation also **leaked a reading from after the forecast origin** into the history-end slot of **2.2% of 30-minute and 4.1% of 60-minute scored windows** — up to 25 minutes ahead, never the target itself. That is a leakage, not a marginal discrepancy. The primary numbers stand as published, with this stated beside them; the conforming variant uses nothing after the origin, and it is the variant the comparisons use.
- **The published comparison column attaches to the `bglp` 2020 rows**, since that is the rule the published figures were computed under; see [docs/published_comparison.md](published_comparison.md) for the eight ranked systems, their inputs, and which qualify.
