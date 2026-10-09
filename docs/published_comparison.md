# Published comparison — OhioT1DM, BGLP Challenge 2020

Literature only. Transcribed by Suhas, 9 Oct 2026, from the official results page and the eight system papers; committed as-is except where marked **(amended)**, with the repository's conforming figures appended as Table 3.

Sources:
- Official rankings: https://webpages.charlotte.edu/rbunescu/data/ohiot1dm/bglp/bglp-results.html
- Challenge rules: https://webpages.charlotte.edu/rbunescu/data/ohiot1dm/bglp/bglp-rules.html
- Proceedings: CEUR-WS Vol-2675, https://ceur-ws.org/Vol-2675/

The official ranking is the April 2020 submissions: 8 of 16 systems conformed to the rules and were scored on the same test points for contributors 540, 544, 552, 567, 584 and 596 (2884 / 2704 / 2352 / 2377 / 2653 / 2731 points). The overall score is the sum of the four metrics. Mean RMSE and MAE are the mean of the six per-patient values.

## Table 1 — the eight officially ranked systems

| Rank | System (first author, paper) | Off/on | Inputs | Input class | RMSE 30 | MAE 30 | RMSE 60 | MAE 60 | Overall | Per-patient | Baseline reported | Qualifies | Why / why not |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Rubin-Falcone (p18) — residual LSTM blocks, N-BEATS-style | Offline, fine-tuned per patient | CGM, finger-stick, bolus, carbs, time of day, missingness flag. **No activity.** Pre-trained on 2018 cohort + Tidepool (~100 people, ~15M steps) | CGM + insulin + meals, **+ external pre-training** | 18.22 | 12.83 | 31.66 | 23.60 | 86.31 | Yes (Table 1) | N-BEATS 21.2 at 30 min; no persistence | **Partial** | Inputs match; external data no other entrant had. Compare with the caveat stated. |
| 2 | Hameed (p14) — vanilla RNN | Offline | CGM only (univariate; the ranked result is Approach I). 2018 six used as validation | CGM only | 19.21 | 13.08 | 31.77 | 23.09 | 87.15 | Yes, 10 repeats | None | **Yes — CGM-only rows** | Note: the paper's preprocessing section lists four variables but the ranked result and conclusion are univariate. |
| 3 | Zhu (p15) — GRU generator with CNN discriminator (GAN) | Offline, personalised | CGM + bolus + carbs; 1.5 h history; test gaps extrapolated only. Pre-training tried and dropped (+0.5 worse) | CGM + insulin + meals | 18.34 | 13.37 | 32.21* | 24.20 | 88.12 | Yes, 10 runs, + Clarke zones per patient | None | **Yes** | Closest match to our primary models. 567 had no meal records, so its input was reduced. *Table says 32.21, abstract 32.31; 32.21 sums to the official overall. |
| 4 | Yang (p24) — multi-scale, multi-lag LSTM | Offline | CGM + basal + bolus + meal + timestamp; GSR/temp/accel explicitly excluded; model-free fallback after ≥12 missing | CGM + insulin + meals | 19.05 | 13.50 | 32.03 | 23.83 | 88.41 | Yes (Table 5) | None | **Yes** | Their post-gap fallback ("remain prediction") is what our conforming variant needs. |
| 5 | Bevan (p17) — single population LSTM, 128 units, 30-min history | Offline, population | CGM only; trained on all 12 (2018 train+test, 2020 train); missing → zero after standardisation | CGM only | 18.23 | 14.37 | 31.10 | 25.75 | 89.45 | **Yes (Table 2)** | Linear/feed-forward in a validation figure only | **Yes — CGM-only rows; direct LOPO comparator** | Patient-excluded 18.32 vs all-patients 18.23 vs patient-only 19.21 (Table 3). 30-min history found optimal. |
| 6 | Sun (p20) — recursive latent-variable (PCA) linear model | **Online** | CGM + basal/bolus as insulin-on-board; no meals; 2 h history | CGM + insulin | 19.37 | 13.76 | 32.59 | 24.64 | 90.36 | Yes, + Clarke zones | None | **Partial** | Online, and no meals. But it is a *linear* method — context for our linear AR. |
| 7 | Joedicke (p25) — MOGE (official entry; numbers match MOGE) | Offline | CGM, basal, bolus type/dose, **GSR, skin temperature**; 2 h history at 15-min steps | + wearables | 19.60 | 14.25 | 34.12 | 25.99 | 93.96 | Yes | **LR 21.22 / 32.88; RF 22.27 / 34.14; ARIMA 24.51 / 38.66 (RMSE 30/60)** | **No** | Wearable inputs. Its classical baselines are the most useful published context on these test points. |
| 8 | Ma (p27) — online ARMA + residual NN | **Online** | CGM only | CGM only | 20.03 | 14.52 | 34.89 | 26.41 | 95.85 | Yes | **AR 21.80; ARMA 21.44; BP-NN 33.45 (RMSE 30 only)** | **Partial** | Online. Skipped test points with null inputs (their rule 2). |

**No system reports a persistence baseline.**

### The like-for-like set for the README column
- Against our **CGM + insulin + meals** rows (per-patient and LOPO l2 with covariates, and Saturday's model): **Zhu 18.34 / 32.21** and **Yang 19.05 / 32.03**; Rubin-Falcone 18.22 / 31.66 with the external-pre-training caveat.
- Against our **CGM-only** rows (l2 without covariates): **Hameed 19.21 / 31.77** and **Bevan 18.23 / 31.10**.
- Published classical baselines on the same test points: LR 21.22 / 32.88, ARIMA 24.51 / 38.66, AR 21.80, ARMA 21.44, N-BEATS 21.2.
- Context only (online or wearable inputs): Sun, Ma, Joedicke.

Three of eight match our primary inputs; two match CGM-only; three are context. Say so in the README.

### The camera-ready (July 2020) ranking
Not transcribed. 13 of 16 systems; informational only per the organisers. Add if wanted; the official ranking is the comparison.

## Table 2 — conformity audit of our protocol against the challenge rules

| Rule | Our frozen protocol | Conforms |
|---|---|---|
| Report on 540, 544, 552, 567, 584, 596 only | We report all 12 with cohort rows; **the 2020-cohort rows are the comparable ones** (persistence 24.16 at 30 min, not the 12-patient 23.39) | ✓ using 2020 rows |
| 2018 six used for pre-training only, or not at all | Per-patient fits don't use them; the LOPO population fit trains on all other patients — permitted | ✓ |
| Evaluation begins 60 min after test start (2020) | Verified Day 3 Part A; scored candidates match challenge counts | ✓ |
| Every challenge test point predicted | **No — targets whose history crosses a >30-min gap are dropped** (540: 2,829 of 2,884 scored; 552: 2,275 of 2,352) | ✗ |
| No interpolation of missing data; extrapolation permitted | **No — history interpolated across gaps ≤ 30 min** (targets never) | ✗ |
| No data after the prediction time | **Marginal as drafted ("up to 60 s after T−h"); measured larger: the history-end slot was valued from a reading after T−h in 2.2% of 30-min and 4.1% of 60-min scored windows, up to 1,500 s after the origin, never the target itself (D-027). Verdict as drafted was "minor"; revised by the author on 10 Oct to leakage once the measurement was in.** | ✗ **leakage** |
| RMSE/MAE per person, then mean of six | Mean-of-patients is our primary aggregation; pooled in brackets | ✓ |
| Offline: one model per patient, fixed at test time | Per-patient l2 and LOPO are offline | ✓ |

**Consequence.** Our primary numbers score a subset that excludes the hardest points — those immediately after a sensor outage — so they are optimistic relative to every number in Table 1 by an unmeasured amount. A conforming evaluation variant (no interpolation, zero-order hold or extrapolation, all challenge points predicted with a stated fallback, 2020 cohort, mean of six) is added as a disclosed addition before any sequence model is fitted. Only conforming numbers are placed beside Table 1.

## Table 3 — this repository's conforming figures beside Table 1 (added 10 Oct 2026)

Every row is the `bglp` variant on the 2020 six — every challenge test point predicted, zero-order hold, nothing after the forecast origin, mean of six — from `results/baselines_challenge_scores.csv` (D-027, D-028). Primary-protocol figures are deliberately not placed here. Published column filled only where the inputs match (Table 1, "Qualifies"); Rubin-Falcone is listed as partial because of its external pre-training and is not used as a comparator figure.

| This repository (conforming, 2020 six, mean of six) | Inputs | RMSE 30 | MAE 30 | RMSE 60 | MAE 60 | Overall | Fallback 30 / 60 | Like-for-like published (RMSE 30 / 60) |
|---|---|---|---|---|---|---|---|---|
| p0 persistence | CGM only | 25.48 | 18.06 | 41.42 | 30.43 | **115.39** | 462 / 662 | — |
| l1 extrapolation k=6 | CGM only | 30.92 | 20.50 | 60.95 | 41.69 | **154.06** | 462 / 662 | — |
| l2 linear AR, per-patient, H=24 (CGM only) | CGM only | 22.00 | 15.79 | 36.43 | 27.76 | **101.98** | 1225 / 1431 | Hameed (p14) **19.21 / 31.77**; Bevan (p17) **18.23 / 31.10** |
| l2 linear AR, per-patient, H=24, +insulin+meals | CGM + insulin + meals | 20.77 | 14.79 | 34.12 | 25.71 | **95.41** | 1225 / 1431 | Zhu (p15) **18.34 / 32.21**; Yang (p24) **19.05 / 32.03** |
| l2 linear AR, population LOPO, H=24 (CGM only) | CGM only | 22.24 | 15.88 | 36.54 | 27.72 | **102.38** | 1225 / 1431 | Hameed (p14) **19.21 / 31.77**; Bevan (p17) **18.23 / 31.10** |
| l2 linear AR, population LOPO, H=24, +insulin+meals | CGM + insulin + meals | 21.61 | 15.42 | 35.52 | 26.99 | **99.54** | 1225 / 1431 | Zhu (p15) **18.34 / 32.21**; Yang (p24) **19.05 / 32.03** |

Published classical baselines on the same test points (RMSE 30 / 60): linear regression (Joedicke p25) 21.22 / 32.88 · random forest (Joedicke p25) 22.27 / 34.14 · ARIMA (Joedicke p25) 24.51 / 38.66 · AR (Ma p27) 21.80 / — · ARMA (Ma p27) 21.44 / —.

**Reading it:** under the challenge's own rules the best linear model here is 2.4 mg/dL behind Zhu and 1.7 behind Yang at 30 minutes, and 1.9–2.1 behind at 60; the CGM-only linear AR is 2.8 behind Hameed and 3.8 behind Bevan at 30. The published linear regression (21.22 / 32.88) sits between this repository's persistence and its linear AR at 30 minutes and is better than both at 60. No ranked system reports persistence; this repository's conforming persistence on these points is 25.48 / 41.42.

---

## For Claude Code

1. Commit this file as `docs/published_comparison.md`. Literature, not data.
2. In the README's Baselines table, fill the published column **only for the 2020-cohort rows and only from the like-for-like set above**, with a one-line note: "three of eight officially ranked systems used inputs matching ours; two match CGM-only; three are online or wearable-based and shown for context."
3. Add the published classical baselines (LR, ARIMA, AR, ARMA) as a separate context row beneath persistence, labelled as published, same test points.
4. Fill the audit's blanks from `docs/protocol.md` if any wording differs; do not change the verdicts.
5. Build the conforming variant per `bglp-comparison-table.md` and rerun the baselines under it; place the conforming 2020-cohort figures beside Table 1. Report primary and conforming side by side — the gap is the measured cost of the subset.
6. Decision-log entries: the transcription and its two corrections (Rubin-Falcone inputs; Bevan per-patient), the like-for-like rule as applied, the four conformity findings, and the variant's definition.
