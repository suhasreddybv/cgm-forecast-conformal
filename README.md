# cgm-forecast-conformal

Glucose forecasting at 30 and 60 minutes with calibrated conformal prediction intervals, evaluated by clinical error grids rather than RMSE alone.

**Result so far.** Forecasting the most recent reading — persistence — scores **23.4 mg/dL RMSE at 30 minutes and 38.5 at 60** on the test split (mean of 12 patients; n = 30,912 and 30,579 real targets). A per-patient linear autoregression on two hours of history plus insulin and meals brings that to **19.0 and 32.0**, a margin of 4.4 and 6.6 mg/dL, and a leave-one-patient-out population fit is within 0.2 of it at 30 minutes. **The sequence model and the conformal intervals are pending**; the evaluation protocol was frozen before any of this was fitted ([docs/protocol.md](docs/protocol.md)), and the bar the sequence model has to clear is the linear model at 19.0, not persistence at 23.4. Under the challenge's own rules — every test point predicted, nothing after the forecast origin — the same linear model scores 20.8 on the 2020 cohort and persistence 25.5; both variants are reported for every method. No published comparison appears until the inputs and evaluation rule are matched (D-012).

## The data layer

12 patients, two cohorts of six, eight weeks each, CGM every five minutes. **Every patient reconciles exactly with Table 2 of the dataset paper**: all 12 training counts match, the 2018 test files match their published counts, and every 2020 test file carries exactly 12 more readings — the first hour, which the 2020 challenge excludes from evaluation and which this repository therefore excludes too.

The cohorts do not share a sensor-band channel set. The 2018 cohort wore the Basis Peak and has heart rate, GSR, skin temperature, air temperature and steps at a 5-minute aggregation; the 2020 cohort wore the Empatica Embrace and has GSR, skin temperature and acceleration at a 1-minute aggregation. The absent channels are represented inconsistently — present-but-empty for five 2020 patients, absent entirely for 596 — so the loader treats missing and empty alike and asserts the difference in both directions.

**Gaps are the thing that would corrupt this silently** (`results/gap_analysis.csv`). The CGM series is **11.9% incomplete on average in training and 10.7% in test**; the worst patient/split, 552/test, is **40.2% missing**; and the longest single gap is **118 hours**. Cut naively, 1.6% of 30-minute history windows cross a gap, 3.4% of 1-hour windows and **7.1% of 2-hour windows** — 13.9% for the worst patient. The frozen rule: history may be interpolated across at most 30 minutes, longer gaps drop the window, and **a target is never interpolated** — a forecast is scored only against a real CGM reading.

**Windows are built on that rule** (`src/data/windows.py`, D-013 to D-016). One window per real reading with enough history; history slots sit at exactly 30 or 60 minutes before the target and are valued by interpolation in time, because the cadence is not exactly 300 s (503 of 166,000 intervals are 301–360 s). Covariates the primary models may use — basal with temporary overrides, boluses with extended ones spread over their delivery, carbohydrates — are aligned to the same slots and see nothing after the history end. **`results/effective_n.csv` is the n column that accompanies every result from here on**: 92–97% of raw readings are evaluable depending on horizon and history length, and 552's test split, 40% missing, keeps 89% of its own readings but only **53% of the readings its sensor should have produced**. Reconciling the windower against the naive gap count found that index-space windows miss every target whose latest reading is older than it looks — the true cost of gaps at 30 minutes is about twice the naive figure (D-016).

_The hero figure — Clarke and Parkes error grids at both horizons — is Week 5 work and is not embedded until it exists._

## Baselines

Written before any sequence model exists. Test split only, scored targets only (D-017), RMSE and MAE in mg/dL as mean of per-patient values with the pooled figure in brackets, MAPE beside them, and the number of real targets under every figure. Per-patient figures, including those for each cohort, are in `results/baselines_per_patient.csv`; the cohort rows are in `results/baselines.csv`.

**The persistence floor is 23.4 mg/dL RMSE at 30 minutes and 38.5 at 60** — the error of forecasting the most recent reading, which every other method also sees. That is the number Saturday's model has to beat, and the first fitted model beats it by 4.4 and 6.6 mg/dL.

| Baseline | Fit | H | 30 min RMSE (pooled) | MAE | MAPE | n | SD / worst | 60 min RMSE (pooled) | MAE | MAPE | n | SD / worst | Published (30 / 60) |
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
- **Published column:** left blank by design. It is filled from the reading table, only where cohort, inputs and evaluation rule match (D-012) — never from memory.
- **Sanity checks passed** (D-021): persistence at 30 minutes sits between 19.2 and 28.5 mg/dL for every patient; 60-minute error exceeds 30-minute error in all 216 baseline × patient comparisons; k=2 extrapolation is noisier than persistence at 60 minutes for every patient.
- Nothing was selected on the test split. All H and all k are reported; no row is a "best".

**Clinical metrics on the same two rows** (`results/clinical_baselines.csv`; Clarke and Parkes grids, MARD, time-in-range agreement, hypoglycaemia sensitivity at 70 and 54 mg/dL, all with n). **Persistence puts 99.0% of 30-minute forecasts in Clarke zones A+B** — an error grid alone does not separate a trivial method from a model. The linear AR improves every aggregate (Clarke A 83.7% → 89.6%, MARD 11.3% → 9.2%) **and detects hypoglycaemia worse**: sensitivity at 70 mg/dL falls from 0.578 to 0.413 at 30 minutes and from 0.356 to 0.104 at 60, because least squares shrinks toward a mean that is not hypoglycaemic. Sensitivity is reported beside every grid from here on, and the 54 mg/dL figures carry their censoring flag in the row (D-023, D-024). The grid figures are in `figures/` labelled as baselines and are not embedded: the hero image is the model-versus-baseline comparison, which does not exist yet.

## Two evaluation variants, reported side by side

The frozen protocol (`primary`) drops a target whose history crosses a gap over 30 minutes and values history slots by linear interpolation. The Blood Glucose Level Prediction Challenge scores **every** test point after the first hour and permits nothing after the forecast origin. On 9 October, before any sequence model was fitted, a conforming variant (`bglp`) was added beside the primary rules — not in place of them (protocol, "Disclosed addition"; D-027): zero-order hold for every slot, every challenge point predicted (2,884 / 2,704 / 2,352 / 2,377 / 2,653 / 2,731 for the 2020 six, asserted before scoring), windows the primary rule would have dropped predicted from the held reading and counted as *fallback*, and the challenge's own scoring — per-patient RMSE and MAE, mean of six, four-score sum. Every method is reported under both from here on; the gap between them is the measured cost of the primary rule's two departures.

| Baseline | Cohort | Variant | RMSE 30 | MAE 30 | RMSE 60 | MAE 60 | Four-score sum | n (30 min) | Fallback (30 / 60) |
|---|---|---|---|---|---|---|---|---|---|
| **p0** persistence | 2020 | primary | 24.16 | 17.62 | 40.27 | 29.92 | **111.98** | 15,239 | — |
|  | 2020 | bglp | 25.48 | 18.06 | 41.42 | 30.43 | **115.39** | 15,701 | 462 / 662 |
|  | 2018 | primary | 22.61 | 16.41 | 36.81 | 27.49 | **103.31** | 15,673 | — |
|  | 2018 | bglp | 23.06 | 16.66 | 37.16 | 27.79 | **104.67** | 15,970 | 297 / 430 |
| l1 extrapolation, k=6 | 2020 | primary | 28.10 | 19.53 | 57.14 | 40.45 | **145.22** | 15,239 | — |
|  | 2020 | bglp | 30.92 | 20.50 | 60.95 | 41.69 | **154.06** | 15,701 | 462 / 662 |
|  | 2018 | primary | 28.20 | 18.70 | 54.58 | 37.60 | **139.07** | 15,673 | — |
|  | 2018 | bglp | 29.25 | 19.24 | 55.93 | 38.30 | **142.72** | 15,970 | 297 / 430 |
| l2 linear AR, per-patient, H=24 | 2020 | primary | 19.79 | 14.51 | 34.83 | 26.53 | **95.65** | 14,476 | — |
|  | 2020 | bglp | 22.00 | 15.79 | 36.43 | 27.76 | **101.98** | 15,701 | 1225 / 1431 |
|  | 2018 | primary | 20.10 | 14.09 | 32.74 | 24.33 | **91.26** | 15,202 | — |
|  | 2018 | bglp | 20.80 | 14.74 | 33.43 | 25.13 | **94.11** | 15,970 | 768 / 898 |
| l2 linear AR, per-patient, H=24, +cov | 2020 | primary | 18.61 | 13.56 | 32.37 | 24.40 | **88.93** | 14,476 | — |
|  | 2020 | bglp | 20.77 | 14.79 | 34.12 | 25.71 | **95.41** | 15,701 | 1225 / 1431 |
|  | 2018 | primary | 19.42 | 13.47 | 31.57 | 23.27 | **87.73** | 15,202 | — |
|  | 2018 | bglp | 20.12 | 14.07 | 32.23 | 23.98 | **90.40** | 15,970 | 768 / 898 |
| l2 linear AR, population LOPO, H=24, +cov | 2020 | primary | 19.55 | 14.35 | 33.98 | 25.91 | **93.78** | 14,476 | — |
|  | 2020 | bglp | 21.61 | 15.42 | 35.52 | 26.99 | **99.54** | 15,701 | 1225 / 1431 |
|  | 2018 | primary | 19.69 | 13.79 | 32.66 | 24.29 | **90.43** | 15,202 | — |
|  | 2018 | bglp | 20.59 | 14.57 | 33.67 | 25.33 | **94.16** | 15,970 | 768 / 898 |

Full table for every baseline, k and H: `results/baselines_challenge_scores.csv`; per-patient rows with fallback counts: `results/baselines_bglp_per_patient.csv`.

- **Conforming is harder, and more so for the 2020 cohort.** Persistence's four-score sum rises from 111.98 to 115.39 on the 2020 six and the best linear AR's from 88.93 to 95.41: about 2 mg/dL of RMSE at 30 minutes for the fitted model, 1.3 for persistence. The patients with the most fallback windows pay the most — 584 and 567 lose 4–5 mg/dL of 30-minute RMSE on the linear AR; 563, with 29 fallbacks, loses nothing.
- **What the primary rule had been buying, and one thing it should not have.** Dropping post-outage targets is a disclosed choice. But the primary interpolation also valued the history-end slot from the reading *after* the forecast origin in **2.2% of 30-minute and 4.1% of 60-minute scored windows** (at most 25 minutes ahead, never the target itself). The primary numbers stand as published, with this stated; the conforming variant uses nothing after the origin.
- **The published comparison column attaches to the `bglp` rows**, since that is the rule the published figures were computed under. It stays blank until filled from the reading table (D-012).

## Method

_To be written (5–18 Oct 2026)._ Baselines and their clinical metrics above; the clinical scorer and the model harness are frozen before the model (D-022 to D-026). Planned next: sequence model (Saturday 10 Oct); split-conformal intervals with empirical vs nominal coverage; MARD, Clarke and Parkes zones, time-in-range agreement, hypoglycaemia sensitivity at 70 and 54 mg/dL. Error-grid figure is the hero image.

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# obtain OhioT1DM first: see data/README.md
pytest                              # 337 tests; 207 need the dataset and skip without it
python -m src.data.gap_analysis     # results/gap_analysis.csv
python -m src.data.effective_n      # results/effective_n.csv and window_reconciliation.csv (~3 min first run)
python -m src.eval.baselines        # baselines under both variants: results/baselines*.csv and baselines_challenge_scores.csv (~6 min)
python -m src.eval.clinical_baselines  # results/clinical_baselines.csv and the two error-grid figures
```

Data: see [data/README.md](data/README.md). No data is included in this repository.

## Limitations and failure cases

Known before any model was fitted. Horizon-dependent coverage failures, hypoglycaemia misses and per-patient variation are added once results exist.

- **The reference is a CGM, not a laboratory measurement.** Ground truth here is a Medtronic Enlite sensor, so every result — this repository's and every published figure on OhioT1DM — measures agreement with a device that itself differs from the patient's true blood glucose. A forecast matching the CGM exactly would still be wrong by the sensor's own error against a laboratory assay. The Enlite was an **adjunctive** device, not cleared for treatment decisions without a confirmatory fingerstick; that is specific to this sensor and era, and current CGMs (Dexcom G6/G7, FreeStyle Libre 2/3) carry non-adjunctive labelling under the iCGM special controls at 21 CFR 862.1355. Results here are described as agreement with the CGM reference, never as accuracy against blood glucose.
- **Readings at the sensor floor are censored, and the bias lands where it matters most.** CGM values are clipped to [40, 400] mg/dL. Both hypoglycaemia thresholds (70 and 54) sit inside that range and can be evaluated; what cannot be measured is error *below* 40, where a true 30 is reported as 40 and a forecast of 40 scores as exact while being 10 mg/dL high. 206 readings (0.124%) sit exactly at the floor — but **17.8% of all readings below 54 mg/dL do**, so low-end error-grid zones and any metric on floored readings are biased optimistic in the severe-hypoglycaemia region.
- **Per-patient figures are not comparable without their target counts.** 552's test split is 40.2% missing with a 118-hour gap, so its error rests on far fewer evaluable targets than any other patient's. Every per-patient number is reported with its effective n.
- **Wearable channels are not uniformly available**, so they are restricted to a secondary stratified analysis: heart rate is empty for five 2020 patients and absent for 596, and acceleration exists only in the 2020 cohort. Primary models use CGM, insulin and meals only — the channels all 12 patients have and the published comparisons use.
- **12 patients, one sensor, one pump family, eight weeks each.** Nothing here supports a claim about other devices or populations.

## References

- Marling, C., & Bunescu, R. (2020). The OhioT1DM dataset for blood glucose level prediction: Update 2020. *CEUR Workshop Proceedings*, 2675, 71–74.
- Clarke, W. L., et al. (1987). Evaluating clinical accuracy of systems for self-monitoring of blood glucose. *Diabetes Care*, 10(5), 622–628.
- Parkes, J. L., et al. (2000). A new consensus error grid to evaluate the clinical significance of inaccuracies in the measurement of blood glucose. *Diabetes Care*, 23(8), 1143–1148.

## Licence

MIT
