<!-- Part of cgm-forecast-conformal; the front page is ../README.md -->

# Week 5 — intervals around the frozen model (first pass, 12 Oct)

The model is frozen at `75e0817`; nothing here retrains, re-selects or refits it. Conformal prediction wraps it. **Calibration set:** the frozen 20% validation tail of each training file — which early stopping consumed, so its residuals are mildly optimistic. Rather than argue about the contamination, it is measured: empirical test coverage against nominal, with an **in-sample bracket** (calibration on the fit portion's residuals, known optimistic) beside it. The frozen model's predictions on its training file came from regenerating all 144 ensemble members from their recorded seeds and starting weights and asserting, member by member, that the regenerated test predictions equal the cached ones to the last bit (D-045). Four methods, none retraining (D-046): split conformal (c0), residuals scaled by the volatility of the history strictly before the forecast origin (c1), Mondrian by forecast bin — below 70, 70–180, above 180 (c2), and adaptive conformal inference, online (c3); a fifth, a one-sided lower bound on the signed residual (c4), was added on 15 October after the breakdown (D-062). Conformalised quantile regression is excluded: the model is a frozen point forecaster.

![Empirical coverage against nominal for the four methods, conforming evaluation, 2020 cohort](../figures/coverage_calibration.png)

| Method | h | 80% coverage (bracket) | 90% coverage (bracket) | 95% coverage (bracket) | 90% median width | 90% useless |
|---|---|---|---|---|---|---|
| c0 split conformal | 30 | 0.814 (0.791) | 0.909 (0.893) | 0.956 (0.949) | 63 mg/dL | 16.9% |
| c1 scaled by history volatility | 30 | 0.806 (0.783) | 0.898 (0.887) | 0.950 (0.943) | 56 mg/dL | 16.4% |
| c2 Mondrian by forecast bin | 30 | 0.815 (0.792) | 0.910 (0.894) | 0.954 (0.947) | 61 mg/dL | 16.9% |
| c3 adaptive (online) | 30 | 0.799 (0.799) | 0.900 (0.900) | 0.950 (0.950) | 51 mg/dL | 5.6% |
| c0 split conformal | 60 | 0.812 (0.772) | 0.900 (0.872) | 0.948 (0.930) | 107 mg/dL | 65.4% |
| c1 scaled by history volatility | 60 | 0.787 (0.754) | 0.897 (0.876) | 0.950 (0.937) | 99 mg/dL | 49.6% |
| c2 Mondrian by forecast bin | 60 | 0.814 (0.770) | 0.898 (0.869) | 0.947 (0.930) | 104 mg/dL | 57.9% |
| c3 adaptive (online) | 60 | 0.801 (0.801) | 0.900 (0.900) | 0.950 (0.950) | 85 mg/dL | 30.5% |

*Conforming evaluation, 2020 cohort pooled over points (n = 15,701 at 30 min, 15,701 at 60). Bracket = the same method calibrated on the fit portion. "Useless" = width over 100 mg/dL or infinite. Every patient, both variants, all levels: `results/coverage.csv`.*

**What the first pass says.** Pooled, split conformal is on nominal: **0.909 / 0.900 at the 90% level** against the in-sample bracket's 0.893 / 0.872 — the validation tail's contamination by early stopping cost almost nothing at 30 minutes and about three points at 60, which the bracket shows and the validation calibration removes. **But a 90% interval an hour ahead is 107 mg/dL wide at the median and 65% of them are wider than 100**: the honest 60-minute interval is mostly too wide to act on. The adaptive method narrows that to 85 with 31% useless, at the cost of some infinite intervals after outages. **The stop condition fired on two patient sets** (c0 at 90% more than five points under nominal): 596 at 30 minutes (0.846) and 540 at 60 (0.842). For 596 the in-sample bracket covers *better* (0.912) — its test residuals are simply larger than its training residuals, an exchangeability failure, not an early-stopping artefact. c3 restores both to nominal (0.901 at 44 mg/dL; 0.899 at 93 mg/dL, 45% useless), which is what it is for (D-048).

![Interval widths by method at the 90% level](../figures/interval_widths.png)

**Hypoglycaemia onset detection, defined correctly** (D-047). Week 4's sensitivity counted every target below 70 mg/dL, including those already below it at the forecast origin — which persistence "detects" for free, so its 0.353 was not a detection rate. An **onset** is a target below the threshold whose origin reading was at or above it; persistence's onset sensitivity is zero by construction. Reported with specificity, PPV and false alarms per day, never alone:

| Detector (conforming, 2020 pooled) | h | Onsets | Sensitivity | Specificity | PPV | False alarms / day |
|---|---|---|---|---|---|---|
| persistence (zero by construction) | 30 | 230 | 0.000 | 1.000 | — | 0.00 |
| linear AR, point forecast | 30 | 230 | 0.083 | 0.998 | 0.396 | 0.55 |
| final model, point forecast | 30 | 230 | 0.070 | 0.999 | 0.432 | 0.40 |
| **c0 lower bound at 90%** | 30 | 230 | 0.904 | 0.855 | 0.087 | 41.15 |
| persistence (zero by construction) | 60 | 328 | 0.000 | 1.000 | — | 0.00 |
| linear AR, point forecast | 60 | 328 | 0.012 | 0.998 | 0.121 | 0.55 |
| final model, point forecast | 60 | 328 | 0.000 | 0.999 | 0.000 | 0.19 |
| **c0 lower bound at 90%** | 60 | 328 | 0.875 | 0.719 | 0.064 | 79.28 |

*Conforming, 2020 cohort, 70 mg/dL, 90% level for the lower bound; the at-risk set is every target whose origin reading was ≥ 70. At 54 mg/dL: 37 onsets at 30 minutes, the point forecast detects 0, 2 at the sensor floor. `results/onset_first_pass.csv`.*

**The point forecast detects 7.0% of 30-minute onsets and 0.0% of 60-minute ones.** The 90% lower bound detects 90.4% and 87.5% — by alarming on **41 and 79 non-events per day**, with a positive predictive value of 0.09 and 0.06. At this level the lower bound is not an alarm; it is the bottom of a wide interval. Whether any level or method sits on a usable point of the sensitivity–false-alarm trade-off is Saturday's operating curve, which is now the figure this repository is building toward.

## Operating points at two false alarms per day

| Method | h | Setting at ≤ 2/day | Achieved false alarms/day | Sensitivity | Specificity | PPV | Interpolated sensitivity at 2/day | Verdict vs sweep (D-064) |
|---|---|---|---|---|---|---|---|---|
| **point forecast − δ (reference)** | 30 | δ = 4 mg/dL | 1.18 | 0.170 | 0.996 | 0.386 | 0.266 | — |
| linear AR − δ | 30 | δ = 6 mg/dL | 1.40 | 0.196 | 0.995 | 0.378 | 0.245 (-0.021 [-0.076, 0.025]) | within variability |
| c0 split | 30 | level 0.30 | 1.88 | 0.248 | 0.993 | 0.365 | 0.262 (-0.004 [-0.037, 0.037]) | within variability |
| c1 volatility-scaled | 30 | level 0.30 | 1.82 | 0.222 | 0.994 | 0.347 | 0.233 (-0.033 [-0.111, 0.032]) | within variability |
| c2 Mondrian by forecast | 30 | level 0.30 | 1.63 | 0.230 | 0.994 | 0.381 | 0.267 (+0.001 [-0.036, 0.040]) | within variability |
| **c3 adaptive** | 30 | level 0.25 | 1.55 | 0.261 | 0.995 | 0.422 | 0.296 (+0.030 [-0.026, 0.093]) | within variability |
| c4 signed, one-sided | 30 | level 0.60 | 1.16 | 0.183 | 0.996 | 0.408 | 0.260 (-0.006 [-0.036, 0.023]) | within variability |
| c5 Mondrian by rate | 30 | level 0.30 | 1.71 | 0.243 | 0.994 | 0.384 | 0.271 (+0.005 [-0.060, 0.039]) | within variability |
| **point forecast − δ (reference)** | 60 | δ = 8 mg/dL | 1.32 | 0.037 | 0.995 | 0.146 | 0.076 | — |
| linear AR − δ | 60 | δ = 10 mg/dL | 1.72 | 0.064 | 0.994 | 0.188 | 0.069 (-0.007 [-0.111, 0.079]) | within variability |
| c0 split | 60 | level 0.30 | 1.95 | 0.079 | 0.993 | 0.202 | 0.081 (+0.005 [-0.012, 0.041]) | within variability |
| c1 volatility-scaled | 60 | level 0.25 | 1.70 | 0.018 | 0.994 | 0.062 | 0.021 (-0.055 [-0.102, -0.006]) | below the sweep |
| c2 Mondrian by forecast | 60 | level 0.30 | 1.69 | 0.067 | 0.994 | 0.198 | 0.080 (+0.005 [-0.011, 0.025]) | within variability |
| **c3 adaptive** | 60 | level 0.25 | 1.99 | 0.192 | 0.993 | 0.375 | 0.193 (+0.117 [0.019, 0.179]) | above the sweep |
| c4 signed, one-sided | 60 | level 0.60 | 1.21 | 0.021 | 0.996 | 0.099 | 0.056 (-0.020 [-0.043, 0.017]) | within variability |
| c5 Mondrian by rate | 60 | level 0.25 | 1.42 | 0.037 | 0.995 | 0.138 | 0.075 (-0.001 [-0.046, 0.021]) | within variability |

*`results/operating_points.csv`, `results/operating_curve_claim.csv`; 70 mg/dL, conforming, 2020 pooled. "Setting" is the last point on each curve at or below 2 per day (the grid is discrete, so the achieved rate sits below the mark); the interpolated column is what the claim rule reads, with its difference from the sweep and the six-patient bootstrap interval. Specificity exceeds 0.99 for every method at these rates — false alarms per day is the axis that discriminates. The regulation's (b)(7)(iv) quartet reads off these rows: true alert rate = sensitivity, missed = 1 − sensitivity, false alert rate = false alarms per day, correct alert rate = PPV (D-066).*

## The clinical tables, computed for the intervals

| Estimate (90%, conforming, 2020) | h | Both bounds in Clarke A+B | Any bound in C/D/E | MARD, lower / upper | Entirely in 70–180 | Entirely out | **Straddles** | Onset sens / FA per day at 70 (this setting) |
|---|---|---|---|---|---|---|---|---|
| point forecast | 30 | 98.4% | 1.6% | 9.5 | 70.68% | 29.32% | **0.0%** | 0.070 / 0.4 |
| persistence | 30 | 98.6% | 1.4% | 12.5 | 68.8% | 31.2% | **0.0%** | 0.000 / 0.0 |
| c0 | 30 | 95.6% | 4.4% | lower 23.8 / upper 27.2 | 32.7% | 13.99% | **53.31%** | 0.904 / 41.1 |
| c1 | 30 | 92.6% | 7.4% | lower 24.3 / upper 27.9 | 37.27% | 14.04% | **48.68%** | 0.913 / 39.6 |
| c2 | 30 | 95.8% | 4.2% | lower 23.0 / upper 26.5 | 34.52% | 13.04% | **52.44%** | 0.900 / 38.5 |
| c3 | 30 | 94.0% | 6.0% | lower 19.1 / upper 22.6 | 42.43% | 15.21% | **42.36%** | 0.917 / 29.2 |
| c4 (lower bound only) | 30 | 99.3% | 0.7% | 15.8 | — | — | **—** | 0.787 / 18.2 |
| c5 | 30 | 95.0% | 5.0% | lower 23.0 / upper 26.4 | 35.04% | 14.08% | **50.88%** | 0.913 / 37.9 |
| point forecast | 60 | 95.9% | 4.1% | 16.9 | 73.56% | 26.44% | **0.0%** | 0.000 / 0.2 |
| persistence | 60 | 96.1% | 3.9% | 21.2 | 68.88% | 31.12% | **0.0%** | 0.000 / 0.0 |
| c0 | 60 | 88.1% | 11.9% | lower 36.8 / upper 45.8 | 7.95% | 7.07% | **84.98%** | 0.875 / 79.3 |
| c1 | 60 | 76.2% | 23.8% | lower 42.8 / upper 52.1 | 18.46% | 6.88% | **74.66%** | 0.939 / 88.8 |
| c2 | 60 | 87.9% | 12.1% | lower 36.0 / upper 45.1 | 10.2% | 6.6% | **83.2%** | 0.854 / 76.0 |
| c3 | 60 | 85.9% | 14.1% | lower 28.6 / upper 38.0 | 24.86% | 7.67% | **67.47%** | 0.927 / 57.2 |
| c4 (lower bound only) | 60 | 96.6% | 3.4% | 24.8 | — | — | **—** | 0.658 / 36.3 |
| c5 | 60 | 87.5% | 12.5% | lower 36.0 / upper 45.0 | 9.43% | 6.8% | **83.77%** | 0.875 / 76.1 |

*`results/interval_clinical.csv`, 90% level. "Both bounds in A+B" is the share of intervals whose lower and upper bound both fall in Clarke zone A or B; infinite c3 intervals count as outside. A bound is not an estimate, which its MARD shows. **Straddles** is the honest indeterminate rate: at 60 minutes 85% of c0's intervals say neither "in range" nor "out" (D-067).*
