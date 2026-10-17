<!-- Part of cgm-forecast-conformal; the front page is ../README.md -->

# The sequence model, built as a staircase

**A sequence model with the same inputs as a linear autoregression does not beat it until it is trained on the full horizon, ensembled over history lengths and pre-trained on the other patients; those steps were worth 0.25, 0.25 and 0.34 mg/dL respectively at 30 minutes (conforming, 2020), the choice of family a further 0.12 and the ensemble of the final model 0.40, of which pre-training and both ensembles exceeded seed variability (0.15) and configuration selection (-0.01) and the family choice did not.** At 60 minutes the full-horizon loss alone was worth 1.36.

Every step keeps what the step before chose; every choice is made on the frozen validation tail (D-025), never on test; every step is scored under both evaluation variants on every patient (`results/sequence_steps.csv`, `results/sequence_per_patient.csv`, `results/sequence_paired.csv`, seeds mean ± SD). Inputs are frozen at CGM history plus basal (with temporary overrides), boluses (spread over delivery) and carbohydrates — no wearable channels, no finger-stick, no time of day (D-032).

| Step | Conforming 30 min | Conforming 60 min | Primary 30 min | Primary 60 min | Four-score (conforming) | Position |
|---|---|---|---|---|---|---|
| Step 1 — single configuration: GRU 64, 1 layer, dropout 0.2, H = 24, 3 seeds | 20.57 ± 0.15 | 34.17 ± 0.03 | 18.79 ± 0.05 | 33.01 ± 0.08 | 94.50 | 8 of 9 |
| Step 2 — validation-selected configuration (coordinate search, 1 seed) | 20.58 | 34.84 | 18.70 | 33.19 | 95.74 | 8 of 9 |
| Step 3 — full-horizon loss, kept where it won on validation (1 seed) | 20.33 | 33.49 | 18.45 | 31.90 | 92.75 | 7 of 9 |
| Step 4 — chosen configuration, 5 seeds | 20.39 ± 0.15 | 34.07 ± 0.21 | 18.56 ± 0.10 | 32.45 ± 0.15 | 93.86 | 7 of 9 |
| Step 4e — history-length ensemble, median over H ∈ {12…42} | 20.13 | 33.66 | 18.32 | 32.10 | 92.75 | 7 of 9 |
| Step 5 — pre-trained on the twelve training files + 2018 test files, fine-tuned per patient, 5 seeds | 20.05 ± 0.14 | 33.29 ± 0.52 | 18.21 ± 0.09 | 31.59 ± 0.54 | 91.48 | 7 of 9 |
| Step 5p — population model alone, leave-one-patient-out (1 seed) | 20.59 | 34.94 | 18.90 | 33.56 | 95.75 | 8 of 9 |
| Step 6 — dilated causal CNN, validation-selected, 5 seeds | 20.35 ± 0.12 | 33.76 ± 0.20 | 18.67 ± 0.08 | 32.55 ± 0.20 | 93.28 | 7 of 9 |
| Step 6e — TCN history-length ensemble | 20.12 | 33.27 | 18.35 | 31.83 | 91.90 | 7 of 9 |
| Step 7 — final model: per patient the family that won on validation, pre-trained and fine-tuned, 5 seeds | 19.93 ± 0.12 | 33.42 ± 0.21 | 18.18 ± 0.07 | 31.91 ± 0.18 | 91.44 | 7 of 9 |
| **Step 7e — final model, history-length ensemble (frozen)** | **19.53** | **32.43** | 17.64 | 30.82 | 88.96 | 5 of 9 |
| Step 7a — final configuration without covariates, 5 seeds | 20.85 ± 0.06 | 34.96 ± 0.29 | 19.12 ± 0.07 | 33.62 ± 0.30 | 95.86 | 9 of 9 |

*2020 cohort, mean of six patients, RMSE in mg/dL; ± is the SD over seeds where more than one seed ran. Position is in the official BGLP 2020 ranking of eight systems plus this one (`results/sequence_ranking.csv`). Linear AR for reference: 20.77 / 34.12 conforming, 18.61 / 32.37 primary.*

**Did the networks learn beyond their persistence initialisation?** Every model starts exactly at persistence (zero-initialised change head), so a fit whose best validation epoch is 1 has learned one epoch's worth. In Steps 1–4 they did learn: the median best epoch in Step 1 is 22, and 64 of 72 fits improved past epoch 1. The epoch-1 stops in those steps are **50 fits, 49 of them for two patients — 591 (28) and 570 (21), plus one for 588** — whose validation tails rise from the first epoch; per-patient early stopping on those two splits is a limitation, and the split is frozen and not changed. **The fine-tuning steps are different**: in Step 7, 74 of 120 fine-tunes and 78 of 120 ensemble fine-tunes stopped at epoch 1 (62 of 120 in Step 5), while the pre-trained population models they start from peaked around epoch 6. Fine-tuning from pre-trained weights barely moved on validation. **So the final model is, in practice, a population model trained on all twelve training files and the 2018 test files — the patient's own training file included — with one epoch of patient adaptation, ensembled over history lengths.** That is Bevan's all-patients scheme (18.23 / 31.10), not Rubin-Falcone's per-patient fine-tuning, and the comparison is drawn accordingly. The leave-one-patient-out population model, which has never seen the evaluated patient, scores 20.59; the "fine-tuned" one 20.05: most of that 0.54 is having seen the patient's training data, not the fine-tuning epochs. The fine-tuner was never offered zero epochs — the pre-trained state itself was not a candidate — which is a limitation of the harness, listed below. Full distribution: `results/best_epochs.csv`, [figures/best_epochs.png](../figures/best_epochs.png).

**Step 1, stated plainly.** A paired gain of +0.36 [0.01, 0.67], 5 of 6 at 30 minutes, with nothing at 60 and −0.03 [−0.32, 0.27] under the primary variant, is marginal, not a gain: *a single-configuration GRU with the same inputs as the linear autoregression does not meaningfully beat it.* The pre-registered expectation of 18.5–19.5 for this model (D-030) was not met; the check was recorded before any later step ran (D-034).

**The second family.** A dilated causal CNN, searched with the same discipline, won on validation for 9 of 24 patient-horizons and is the final family for those; its own ensemble scored 20.12 / 33.27, within seed variability of the recurrent one.

**Per patient, conforming, with the fallback count beside the penalty** — does the model lose on 567 and 584 the way the linear AR did? Yes, by almost the same amount:

| Patient | n | Fallback 30 / 60 | Final, primary → conforming, 30 min | penalty | Final, primary → conforming, 60 min | Linear AR penalty, 30 min |
|---|---|---|---|---|---|---|
| 540 | 2,884 | 145 / 164 | 19.35 → **20.92** | +1.57 | 36.38 → 38.08 | +1.42 |
| 544 | 2,704 | 87 / 105 | 15.94 → **16.60** | +0.66 | 26.76 → 27.28 | +0.46 |
| 552 | 2,352 | 210 / 252 | 15.17 → **15.58** | +0.41 | 28.27 → 28.94 | +0.82 |
| 567 | 2,377 | 261 / 315 | 19.22 → **23.53** | +4.31 | 33.83 → 37.98 | +5.44 |
| 584 | 2,653 | 348 / 391 | 20.34 → **24.78** | +4.44 | 33.50 → 36.32 | +4.39 |
| 596 | 2,731 | 174 / 204 | 15.84 → **15.80** | -0.04 | 26.18 → 25.99 | +0.42 |

*Final model (7e). The linear AR's conforming penalties at 30 minutes were +5.44 (567) and +4.39 (584); the model's are +4.31 and +4.44. The primary-to-conforming gap for the model is 1.89 mg/dL at 30 minutes against the linear AR's 2.16: the cost of the outages is paid by every method alike. 596 is the one patient marginally easier under the conforming rule (-0.04); the cohort-level check that conforming is never easier passed.*

**Covariate ablation.** The final configuration without insulin and meals scores 20.85 ± 0.06 / 34.96 ± 0.29 against 19.93 / 33.42 with them: **the covariates are worth -0.92 mg/dL at 30 minutes and -1.54 at 60** to the sequence model — the same ~0.9 they were worth to the linear model — and without them the sequence model does not beat the linear AR that has them (+0.11 [-0.41, 0.67], 3 of 6).

**Clinical scoring of the final model** (`results/clinical_final.csv`; frozen scorer, D-022–024):

| Method (conforming, 2020 pooled) | h | n | Clarke A / A+B | Parkes A / A+B | MARD | TIR agreement | Sens 70 (events) | Sens 54 (events, at floor) |
|---|---|---|---|---|---|---|---|---|
| persistence | 30 | 15,701 | 80.7 / 98.6 | 83.0 / 99.4 | 12.5% | 0.867 | **0.546** (507) | 0.051 (39, 2) |
| persistence | 60 | 15,701 | 60.2 / 96.1 | 63.4 / 96.9 | 21.2% | 0.780 | **0.353** (507) | 0.000 (39, 2) |
| linear AR H=24 +cov | 30 | 15,701 | 86.3 / 97.9 | 89.8 / 99.6 | 10.6% | 0.887 | **0.233** (507) | 0.051 (39, 2) |
| linear AR H=24 +cov | 60 | 15,701 | 64.8 / 95.5 | 69.6 / 97.5 | 18.7% | 0.818 | **0.047** (507) | 0.000 (39, 2) |
| **final model (7e)** | 30 | 15,701 | 88.9 / 98.4 | 90.7 / 99.5 | 9.6% | 0.899 | **0.318** (507) | 0.000 (39, 2) |
| **final model (7e)** | 60 | 15,701 | 70.2 / 95.9 | 73.3 / 97.8 | 16.9% | 0.832 | **0.012** (507) | 0.000 (39, 2) |

The error grids separate the three methods by a few points of zone A and not at all on A+B. **Hypoglycaemia sensitivity separates them the other way**: at 30 minutes the model detects 31.8% of reference readings below 70 mg/dL, the linear AR 23.3%, persistence 54.6%; at 60 minutes the model detects **1.2%** against persistence's 35.3%, and at 54 mg/dL it detects none. A model trained on squared error forecasts the conditional mean, and the conditional mean an hour ahead is almost never below 70. Zhu's reported worst case of 8% of 60-minute forecasts in Clarke zone D is not exceeded here — the worst patient is 567 at 7.7% — but the over-prediction-into-hypoglycaemia failure Zhu describes appears as its mirror: the model under-reads the lows. The 54 mg/dL figures rest on 39 events across the cohort, 2 at the sensor floor, and are flagged as censored in the file.

**Not done, and logged as not done (D-038):** no hyperparameter chosen on test; no external data; no test-time adaptation; no wearable inputs; no finger-stick or time-of-day inputs. The remaining gap to Rubin-Falcone is, on the honest reading, attributable to its ~15M steps of external pre-training (Tidepool) and its two extra inputs; the gap to Zhu (1.2 at 30 minutes, 0.2 at 60) is not explained by inputs and belongs to this implementation.
