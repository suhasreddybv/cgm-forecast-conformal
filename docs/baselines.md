<!-- Part of cgm-forecast-conformal; the front page is ../README.md -->

# Baselines

Written before any sequence model exists. Test split only, scored targets only (D-017), RMSE and MAE in mg/dL as mean of per-patient values with the pooled figure in brackets, MAPE beside them, and the number of real targets under every figure. Per-patient figures, including those for each cohort, are in `results/baselines_per_patient.csv`; the cohort rows are in `results/baselines.csv`.

**The persistence floor is 23.4 mg/dL RMSE at 30 minutes and 38.5 at 60** — the error of forecasting the most recent reading, which every other method also sees. That is the floor every later method is measured from: the linear AR beats it by 4.4 and 6.6 mg/dL here, and the final sequence model by 5.9 and 9.0 under the conforming rules on the 2020 cohort.

| Baseline | Fit | H | 30 min RMSE (pooled) | MAE | MAPE | n | SD / worst | 60 min RMSE (pooled) | MAE | MAPE | n | SD / worst | Published |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **p0** persistence | none | 6 | **23.39** (23.55) | 17.01 | 11.4% | 30,912 | 2.88 / 540 28.4 | **38.54** (38.79) | 28.70 | 19.4% | 30,579 | 4.73 / 540 47.5 | — |
| l1 extrapolation, k=2 | none | 6 | 29.86 (30.33) | 19.14 | 12.8% | 30,912 | 5.33 / 575 39.9 | 63.35 (64.07) | 41.93 | 28.2% | 30,579 | 9.52 / 575 80.3 | — |
| l1 extrapolation, k=3 | none | 6 | 28.23 (28.61) | 18.58 | 12.4% | 30,912 | 4.74 / 575 37.9 | 59.12 (59.72) | 40.11 | 27.0% | 30,579 | 8.57 / 575 74.4 | — |
| l1 extrapolation, k=6 | none | 6 | 28.15 (28.43) | 19.12 | 12.8% | 30,912 | 4.17 / 575 35.7 | 55.86 (56.37) | 39.02 | 26.4% | 30,579 | 8.02 / 540 66.8 | — |
| l2 linear AR | per-patient | 6 | 19.95 (20.09) | 14.33 | 9.8% | 30,912 | 2.39 / 575 24.4 | 33.89 (34.09) | 25.55 | 17.9% | 30,579 | 3.85 / 540 41.2 | — |
| l2 linear AR + covariates | per-patient | 6 | 19.42 (19.56) | 13.91 | 9.5% | 30,912 | 2.42 / 575 23.9 | 32.72 (32.90) | 24.57 | 17.2% | 30,579 | 3.77 / 540 39.4 | — |
| l2 linear AR | population, LOPO | 6 | 20.10 (20.23) | 14.48 | 9.9% | 30,912 | 2.27 / 575 23.9 | 34.20 (34.37) | 25.88 | 18.1% | 30,579 | 3.56 / 540 40.8 | — |
| l2 linear AR + covariates | population, LOPO | 6 | 19.79 (19.92) | 14.27 | 9.8% | 30,912 | 2.29 / 575 23.5 | 33.63 (33.78) | 25.46 | 17.8% | 30,579 | 3.42 / 540 40.3 | — |
| l2 linear AR | per-patient | 12 | 19.94 (20.09) | 14.30 | 9.8% | 30,501 | 2.42 / 575 24.4 | 33.89 (34.09) | 25.52 | 17.8% | 30,170 | 3.94 / 540 41.2 | — |
| l2 linear AR + covariates | per-patient | 12 | 19.09 (19.23) | 13.60 | 9.3% | 30,501 | 2.43 / 575 23.6 | 32.23 (32.41) | 24.10 | 16.8% | 30,170 | 3.70 / 540 38.9 | — |
| l2 linear AR | population, LOPO | 12 | 20.11 (20.23) | 14.47 | 9.9% | 30,501 | 2.30 / 575 23.9 | 34.20 (34.35) | 25.86 | 18.0% | 30,170 | 3.63 / 540 40.9 | — |
| l2 linear AR + covariates | population, LOPO | 12 | 19.61 (19.73) | 14.08 | 9.6% | 30,501 | 2.27 / 575 23.3 | 33.36 (33.51) | 25.18 | 17.6% | 30,170 | 3.26 / 540 40.0 | — |
| l2 linear AR | per-patient | 24 | 19.94 (20.09) | 14.30 | 9.8% | 29,678 | 2.50 / 575 24.7 | 33.78 (33.97) | 25.43 | 17.7% | 29,342 | 3.99 / 540 41.0 | — |
| l2 linear AR + covariates | per-patient | 24 | 19.02 (19.17) | 13.51 | 9.2% | 29,678 | 2.46 / 575 23.6 | 31.97 (32.18) | 23.83 | 16.5% | 29,342 | 3.79 / 540 38.9 | — |
| l2 linear AR | population, LOPO | 24 | 20.10 (20.23) | 14.45 | 9.9% | 29,678 | 2.38 / 575 24.0 | 34.09 (34.23) | 25.75 | 17.9% | 29,342 | 3.61 / 540 40.8 | — |
| l2 linear AR + covariates | population, LOPO | 24 | 19.62 (19.75) | 14.07 | 9.6% | 29,678 | 2.34 / 575 23.5 | 33.32 (33.47) | 25.10 | 17.4% | 29,342 | 3.20 / 540 40.0 | — |

- **p0** forecasts the value at the last history slot, T − h. The slot value is the raw reading wherever a reading sits on the slot and within 1/301 of a step of it otherwise (D-013); it is the same "most recent reading" every other method sees.
- **l1** is worse than persistence at every k and both horizons, and two-point extrapolation an hour out is worse than persistence for **12 of 12 patients**: local trend is not the task.
- **l2** is ordinary least squares on the H history slots, coefficients from training windows only. The best configuration is per-patient, H=24, with covariates: **19.02 / 31.97**. Covariates (basal with temporary overrides, boluses spread over delivery, carbohydrates) are worth about 0.9 mg/dL at 30 minutes and 1.8 at 60 for the per-patient fit; history beyond an hour is worth nothing without them (19.94 at H=12 and H=24 alike). The leave-one-patient-out population fit is within 0.2 mg/dL of the per-patient fit at 30 minutes and 0.3–0.5 worse at 60.
- **Published column:** blank here by design — these are 12-patient, primary-protocol figures and no published number was computed that way. The comparison lives in the conforming 2020-cohort table below, where cohort, inputs and evaluation rule match (D-012, D-029).
- **Sanity checks passed** (D-021): persistence at 30 minutes sits between 19.2 and 28.5 mg/dL for every patient; 60-minute error exceeds 30-minute error in all 216 baseline × patient comparisons; k=2 extrapolation is noisier than persistence at 60 minutes for every patient.
- Nothing was selected on the test split. All H and all k are reported; no row is a "best".

**Clinical metrics on the same two rows** (`results/clinical_baselines.csv`; Clarke and Parkes grids, MARD, time-in-range agreement, hypoglycaemia sensitivity at 70 and 54 mg/dL, all with n). **Persistence puts 99.0% of 30-minute forecasts in Clarke zones A+B** — an error grid alone does not separate a trivial method from a model. The linear AR improves every aggregate (Clarke A 83.7% → 89.6%, MARD 11.3% → 9.2%) **and detects hypoglycaemia worse**: sensitivity at 70 mg/dL falls from 0.578 to 0.413 at 30 minutes and from 0.356 to 0.104 at 60, because least squares shrinks toward a mean that is not hypoglycaemic. Sensitivity is reported beside every grid from here on, and the 54 mg/dL figures carry their censoring flag in the row (D-023, D-024). The baseline-only grid figures are in `figures/` and are not embedded; the hero image above is the model-versus-persistence comparison on the same scorer.
