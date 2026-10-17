<!-- Part of cgm-forecast-conformal; the front page is ../README.md -->

# Where coverage breaks

Pooled, every method is on nominal at the 90% level (c0 0.909, c1 0.898, c2 0.910, c3 0.900 at 30 minutes, conforming, 2020 — the Week 5 first pass above); this section is where it stops being. The rules for calling a cell a finding were logged before the tables were read (D-053): **n ≥ 100** or the cell is greyed and never cited; **more than 5 points from nominal at 90% in both evaluation variants**, or it is "variant-specific", not a finding; true-glucose bin at 90%, conforming, 2020 is the headline stratification; and every miss is split into *below* the lower bound (the lower bound too high — the dangerous direction for a hypoglycaemia target) and *above* the upper bound. All 2,304 cells are in `results/coverage_breakdown.csv`; the prose cites by the rule.

![Coverage by true-glucose bin for c0, c2 and c3 at the 90% level](../figures/coverage_by_glucose.png)

| Method | h | <54 | 54-70 | 70-180 | 180-250 | >250 |
|---|---|---|---|---|---|---|
| c0 split conformal | 30 | — (n = 39) | 0.923 · 100% below · n 468 · 63 mg/dL · within 5 | 0.930 · 55% below · n 10,715 · 63 mg/dL · within 5 | 0.880 · 19% below · n 3,386 · 63 mg/dL · within 5 | 0.793 · 8% below · n 1,093 · 63 mg/dL · **finding** |
| c2 Mondrian, by **true** bin | 30 | — (n = 39) | 0.906 · 100% below · n 468 · 59 mg/dL · within 5 | 0.924 · 51% below · n 10,715 · 59 mg/dL · within 5 | 0.893 · 15% below · n 3,386 · 67 mg/dL · within 5 | 0.841 · 5% below · n 1,093 · 67 mg/dL · **finding** |
| c2 Mondrian, by **forecast** bin (what it guarantees) | 30 | — (n = 3) | 0.868 · 16% below · n 242 · 41 mg/dL · within 5 | 0.911 · 36% below · n 11,097 · 59 mg/dL · within 5 | 0.921 · 45% below · n 3,411 · 67 mg/dL · within 5 | 0.883 · 45% below · n 948 · 67 mg/dL · within 5 |
| c3 adaptive | 30 | — (n = 39) | 0.910 · 100% below · n 468 · 56 mg/dL · within 5 | 0.906 · 54% below · n 10,715 · 49 mg/dL · within 5 | 0.887 · 23% below · n 3,386 · 55 mg/dL · within 5 | 0.881 · 13% below · n 1,093 · 71 mg/dL · within 5 |
| c0 split conformal | 60 | — (n = 39) | 0.857 · 100% below · n 468 · 108 mg/dL · within 5 | 0.943 · 68% below · n 10,715 · 107 mg/dL · within 5 | 0.848 · 12% below · n 3,386 · 108 mg/dL · **finding** | 0.668 · 0% below · n 1,093 · 108 mg/dL · **finding** |
| c2 Mondrian, by **true** bin | 60 | — (n = 39) | 0.850 · 100% below · n 468 · 108 mg/dL · within 5 | 0.940 · 65% below · n 10,715 · 99 mg/dL · within 5 | 0.846 · 8% below · n 3,386 · 109 mg/dL · **finding** | 0.678 · 1% below · n 1,093 · 109 mg/dL · **finding** |
| c2 Mondrian, by **forecast** bin (what it guarantees) | 60 | — (n = 0) | — (n = 28) | 0.907 · 29% below · n 11,549 · 99 mg/dL · within 5 | 0.885 · 43% below · n 3,440 · 109 mg/dL · within 5 | 0.822 · 57% below · n 684 · 115 mg/dL · **finding** |
| c3 adaptive | 60 | — (n = 39) | 0.885 · 100% below · n 468 · 77 mg/dL · within 5 | 0.914 · 59% below · n 10,715 · 77 mg/dL · within 5 | 0.864 · 13% below · n 3,386 · 88 mg/dL · within 5 | 0.883 · 7% below · n 1,093 · 121 mg/dL · within 5 |

*Each cell: coverage · share of misses that are below the lower bound · n · median width · verdict under D-053. Conforming evaluation, 2020 cohort, 90% level. The < 54 bin has 39 targets and is not cited.*

**By true glucose, the expectation did not come true where it was expected and came true where it was not.** The constant-width c0 was expected to under-cover below 70 with the misses below. In the 54–70 bin (n = 468) it covers **0.923 at 30 minutes and 0.857 at 60** — within the 5-point rule at both horizons, so *not a finding* — and every miss it does make there is below the lower bound (36 of 36 misses at 30 minutes, 67 of 67 at 60): when the interval fails on a low reading, it fails in the dangerous direction, every time. **The finding under the rule is at the top of the range:** above 250 mg/dL (n = 1,093) c0 covers **0.793 at 30 minutes and 0.668 at 60**, and 180–250 at 60 minutes 0.848 — in both variants, with the misses *above* the upper bound (100% of them at 60 minutes). The model shrinks toward the mean at the top of the range as much as at the bottom, and a constant-width interval cannot follow it there. In the 54–70 bin c3 holds coverage with *less* width than c0 (56 vs 63 mg/dL at 30 minutes, 77 vs 108 at 60) and, like every method, every miss it makes there is below the lower bound (42 / 0, 54 / 0): it is not buying hypoglycaemia coverage with symmetric width, and the misses that remain are all in the dangerous direction. **c3 holds nominal in every reportable bin** (0.881 / 0.883 above 250) by widening online — its median width above 250 at 60 minutes is 121 mg/dL.

**c2, both ways.** Mondrian by forecast bin guarantees coverage per *forecast* bin, which is not the bin a clinician asks about. By forecast bin it holds within the rule at 30 minutes everywhere reportable and fails above 250 at 60 minutes (0.822, n = 684); by *true* bin it fails above 250 at both horizons (0.841 / 0.678) and 180–250 at 60. In the 54–70 true bin it covers 0.906 / 0.850 — within the rule. The forecast bin below 70 is nearly empty (242 forecasts at 30 minutes, 28 at 60): **the model almost never forecasts below 70, which is why conditioning on its forecast cannot protect the readings that end up there.** That is a statement about the model's bias (Week 4's sensitivity result), not about conformal.

![Coverage by rate of change at the origin, c0 against c1](../figures/coverage_by_rate.png)

**Rate of change: fast-moving glucose breaks the constant-width methods, and c1 fixes it at a price.**

| Rate at origin (mg/dL/min) | n (30 / 60) | c0 coverage 30 / 60 | c1 coverage 30 / 60 | c1 median width 30 / 60 | c1 useless 30 / 60 |
|---|---|---|---|---|---|
| <=1 | 12,265 / 12,307 | 0.928 (within 5) / 0.916 (within 5) | 0.903 (within 5) / 0.892 (within 5) | 47 / 86 | 12% / 42% |
| 1-2 | 2,521 / 2,506 | 0.864 (variant-specific) / 0.863 (variant-specific) | 0.887 (within 5) / 0.921 (within 5) | 77 / 140 | 28% / 74% |
| 2-3 | 675 / 661 | 0.809 (**finding**) / 0.782 (**finding**) | 0.874 (within 5) / 0.894 (within 5) | 88 / 161 | 36% / 87% |
| >3 | 240 / 227 | 0.713 (**finding**) / 0.811 (**finding**) | 0.838 (**finding**) / 0.890 (within 5) | 94 / 167 | 46% / 79% |

c0 falls to 0.809 / 0.782 at 2–3 mg/dL/min and 0.713 / 0.811 above 3 — findings in both variants. c1's volatility-scaled width brings every one of those cells back within the rule except above 3 mg/dL/min at 30 minutes (0.838), **and it does so with intervals of 94 / 167 mg/dL at the median, 46% / 79% of them wider than 100.** On fast-moving glucose c1 holds coverage largely by saying it does not know.

![Interval width against absolute error for c1 and c3](../figures/width_vs_error.png)

**Width against error: does the uncertainty inform?** c0's width is constant within a patient, so its within-patient correlation is zero by construction; across the six patients the correlation between a patient's c0 width and its RMSE is 0.66 at 30 minutes and 0.83 at 60 — the patient with the wider interval is the patient with the larger errors. Within patient, **c3's width correlates with the error at a Spearman median of 0.34 (range 0.24–0.41) at 30 minutes and 0.36 at 60; c1's at 0.28 and 0.18; c2's at 0.08 and 0.09.** The online method's intervals are wide where the errors are; the volatility-scaled ones less so; Mondrian's barely at all.

![c3 rolling coverage around the ten largest test-file gaps](../figures/c3_after_gaps.png)

**Recovery after outages.** For the ten largest test-file gaps, 552's 118 hours first, the minutes after the gap closes until c3's rolling one-hour coverage is back within 5 points of nominal:

| # | Patient | Gap (h) | Minutes to recover, 30 / 60 min | Infinite share, first hour, 30 / 60 | Median finite width, first hour, 30 / 60 |
|---|---|---|---|---|---|
| 1 | 552 | 118.0 | 0.0 / 0.0 | 0% / 0% | 47.3 / 76.5 |
| 2 | 544 | 28.9 | 55.0 / 0.0 | 0% / 0% | 43.4 / 64.5 |
| 3 | 584 | 16.9 | 55.0 / 45.0 | 54% / 54% | 92.0 / 269.2 |
| 4 | 567 | 16.0 | 60.0 / 70.0 | 0% / 15% | 68.2 / 101.1 |
| 5 | 567 | 11.3 | 20.0 / 0.0 | 0% / 0% | 140.3 / 140.6 |
| 6 | 596 | 10.7 | 40.0 / 45.0 | 0% / 0% | 71.5 / 80.1 |
| 7 | 563 | 9.9 | 0.0 / 0.0 | 0% / 0% | 55.5 / 51.9 |
| 8 | 540 | 9.6 | 55.0 / 65.0 | 38% / 23% | 66.5 / 94.5 |
| 9 | 559 | 9.0 | 55.0 / 40.0 | 54% / 69% | 68.7 / 248.4 |
| 10 | 559 | 7.4 | 40.0 / 40.0 | 8% / 0% | 99.2 / 115.5 |

*`results/c3_recovery.csv`; "0" means the first hour after the gap was already within 5 points. Widths in mg/dL.*

Every gap recovered within the 24 hours examined; the median is **48 minutes at 30 minutes and 40 at 60**, the longest 60 / 70. After 552's 118-hour gap c3 was within 5 points from the first hour with no infinite intervals and a median width of 47.3 / 76.5 mg/dL. On 4 of the ten gaps at 30 minutes (4 at 60) the first hour contained infinite intervals, up to 54% / 69% of its targets. **For a deployed system: a user coming back from a sensor outage gets intervals that track nominal coverage again within about an hour, and in that hour they are 43–140 mg/dL wide at 30 minutes, 52–269 at 60, with no finite bound at all on some readings after some outages.**

**What this does not show.** Cells with fewer than 100 targets — the whole < 54 mg/dL bin (39 targets at 90%), the forecast bins below 70 — are in the CSV and greyed in the figures, and nothing is claimed from them; the level-2 hypoglycaemia question is answered on Saturday by the onset metric, not by a coverage cell. Mondrian c2 conditions on the forecast bin and is reported by the true bin as well because the two are different questions; a method conditioned on the *true* bin would condition on the target and is not permitted (D-053). None of this is a claim about sensor accuracy: the reference is a CGM, every interval is an interval around agreement with a device ([Limitations](../README.md#limitations-and-failure-cases)). 

## iCGM analogues

These criteria govern a **sensor against a laboratory reference**; this repository compares a **forecast against a CGM** (Medtronic Enlite, adjunctive-labelled), which is itself off from blood glucose by its own error — so a forecast that *passes* an analogue has shown agreement with a sensor, not accuracy against blood, and one that *fails* has failed something easier than the regulation asks ([docs/icgm_criteria.md](icgm_criteria.md) §5, Suhas's transcription of 21 CFR 862.1355(b), verified against the eCFR on (A), (H)/(I) and (J); D-059). Every criterion is binned by the *estimate's* value — "when the device says X, how far can the truth be" — and the regulation requires the lower one-sided 95% bound, not the point estimate; here that bound is optimistic because consecutive five-minute targets are serially correlated. The thresholds are parsed from the document at run time; none is typed into code, and a test checks the parse against the document's rows. `results/icgm_analogues.csv`.

| Estimate (analogue; conforming, 2020, 90% bounds) | h | (A) est < 70, ±15 mg/dL, ≥ 85% | (B) 70–180, ±15%, ≥ 70% | (C) > 180, ±15%, ≥ 80% | (G) all, ±20%, ≥ 87% | (H) / (I) violations (per 1,000) | (J) / (K) wrong direction, ≤ 1% |
|---|---|---|---|---|---|---|---|
| final model, point forecast | 30 | 75% (80%) vs 85% · not met · n 245 | 77% (78%) vs 70% · met · n 11,097 | 88% (89%) vs 80% · met · n 4,359 | 89% (89%) vs 87% · met · n 15,701 | 0 (0.00) / 0 (0.00) | 0.00% / 0.00% |
| final model, point forecast | 60 | not reportable, n = 28 | 55% (55%) vs 70% · not met · n 11,549 | 66% (68%) vs 80% · not met · n 4,124 | 70% (70%) vs 87% · not met · n 15,701 | 0 (0.00) / 0 (0.00) | 0.00% / 0.00% |
| linear AR, point | 30 | 55% (60%) vs 85% · not met · n 224 | 73% (74%) vs 70% · met · n 11,116 | 85% (86%) vs 80% · met · n 4,361 | 86% (86%) vs 87% · not met · n 15,701 | 0 (0.00) / 0 (0.00) | 0.00% / 0.00% |
| linear AR, point | 60 | not reportable, n = 87 | 48% (49%) vs 70% · not met · n 11,533 | 62% (63%) vs 80% · not met · n 4,081 | 64% (65%) vs 87% · not met · n 15,701 | 1 (0.06) / 0 (0.00) | 0.00% / 0.00% |
| persistence | 30 | 63% (67%) vs 85% · not met · n 498 | 67% (68%) vs 70% · not met · n 10,802 | 75% (76%) vs 80% · not met · n 4,401 | 80% (81%) vs 87% · not met · n 15,701 | 0 (0.00) / 0 (0.00) | 0.00% / 0.00% |
| persistence | 60 | 40% (43%) vs 85% · not met · n 490 | 47% (47%) vs 70% · not met · n 10,815 | 53% (54%) vs 80% · not met · n 4,396 | 60% (60%) vs 87% · not met · n 15,701 | 7 (0.45) / 2 (0.13) | 0.00% / 0.00% |
| c0 lower bound | 30 | 6% (7%) vs 85% · not met · n 2,869 | 27% (28%) vs 70% · not met · n 10,636 | 53% (55%) vs 80% · not met · n 2,196 | 43% (44%) vs 87% · not met · n 15,701 | 10 (0.64) / 0 (0.00) | 0.00% / 0.32% |
| c0 lower bound | 60 | 3% (4%) vs 85% · not met · n 4,931 | 12% (12%) vs 70% · not met · n 9,660 | 27% (29%) vs 80% · not met · n 1,110 | 17% (17%) vs 87% · not met · n 15,701 | 122 (7.77) / 0 (0.00) | 0.00% / 0.07% |
| c3 lower bound | 30 | 24% (26%) vs 85% · not met · n 2,004 | 45% (45%) vs 70% · not met · n 11,085 | 64% (65%) vs 80% · not met · n 2,388 | 62% (63%) vs 87% · not met · n 15,477 | 21 (1.36) / 0 (0.00) | 0.00% / 0.38% |
| c3 lower bound | 60 | 11% (12%) vs 85% · not met · n 3,374 | 28% (29%) vs 70% · not met · n 10,735 | 48% (50%) vs 80% · not met · n 1,205 | 36% (37%) vs 87% · not met · n 15,314 | 181 (11.82) / 0 (0.00) | 0.00% / 0.27% |
| c4 lower bound (one-sided) | 30 | 40% (42%) vs 85% · not met · n 1,601 | 52% (53%) vs 70% · not met · n 11,228 | 77% (78%) vs 80% · not met · n 2,872 | 72% (72%) vs 87% · not met · n 15,701 | 6 (0.38) / 0 (0.00) | 0.00% / 0.05% |
| c4 lower bound (one-sided) | 60 | 16% (17%) vs 85% · not met · n 2,523 | 28% (29%) vs 70% · not met · n 11,413 | 44% (46%) vs 80% · not met · n 1,765 | 41% (41%) vs 87% · not met · n 15,701 | 43 (2.74) / 0 (0.00) | 0.00% / 0.01% |
| c0 upper bound | 30 | not reportable, n = 0 | 20% (21%) vs 70% · not met · n 8,004 | 39% (40%) vs 80% · not met · n 7,697 | 42% (43%) vs 87% · not met · n 15,701 | 0 (0.00) / 3 (0.19) | 0.15% / 0.00% |
| c0 upper bound | 60 | not reportable, n = 0 | 12% (13%) vs 70% · not met · n 5,380 | 19% (19%) vs 80% · not met · n 10,321 | 23% (24%) vs 87% · not met · n 15,701 | 0 (0.00) / 50 (3.19) | 0.10% / 0.00% |
| c3 upper bound | 30 | not reportable, n = 0 | 34% (34%) vs 70% · not met · n 8,608 | 50% (51%) vs 80% · not met · n 6,869 | 56% (57%) vs 87% · not met · n 15,477 | 0 (0.00) / 5 (0.32) | 0.31% / 0.00% |
| c3 upper bound | 60 | not reportable, n = 0 | 21% (22%) vs 70% · not met · n 6,673 | 27% (28%) vs 80% · not met · n 8,641 | 34% (35%) vs 87% · not met · n 15,314 | 0 (0.00) / 55 (3.59) | 0.14% / 0.00% |

*Wilson lower one-sided 95% bound (plain share in brackets) against the required percentage; (H) and (I) are zero-tolerance; (J) and (K) are the wrong-trend rules. Cells with n < 100 say so instead of a percentage (D-053). A lower bound is expected to sit below the target by construction and fails the accuracy analogues; its row answers a different question — whether a bound that reads below 70 is ever paired with a reading above 180 (H).*

The contrast the transcription predicted is half there. At 60 minutes the point forecast reads below 70 only 28 times, so (A) and (D) are not reportable for it; at 30 minutes it reads below 70 245 times — just reportable — and **fails (A): 75% within 15 mg/dL against the 85% required**, while meeting (B), (C) and the whole-range (G) (89% against 87%). Persistence fails (A) at 63%. The point forecasts violate neither (H) nor (I) and trip neither trend rule. The c0 lower bound violates (H) 122 times at 60 minutes (a bound below 70 against a reading above 180), c3's 181, c4's 43; the upper bounds violate (I) up to 50 times.

The labelling ranges the regulation requires in (b)(7)(i) — < 54, 54 to < 70, 70–180, > 180 to 250, > 250 — are exactly the true-glucose bins of the coverage breakdown above, so the breakdown table is that labelling, computed on forecasts. The alert-rate quartet of (b)(7)(iv) — true, missed, false and correct alert rates — is exactly the onset-detection table (sensitivity, 1 − sensitivity, false alarms, PPV), which Saturday reports at every nominal level under the regulation's words. Nothing here is a claim that any device meets, or could meet, §862.1355.

**A one-sided bound on the signed residual, c4** (added 15 October, D-062). Every miss in the 54–70 bin was below the bound for every symmetric method in both variants (D-061); a symmetric score cannot correct direction, so c4 calibrates the quantile of the *signed* residual ŷ − y and produces a lower bound only, ŷ − q(ŷ − y), with one-sided coverage nominal by construction. The question was whether a signed score reaches down where the symmetric ones cannot. **It cannot, and the algebra says why: |ŷ − y| ≥ ŷ − y pointwise, so the symmetric quantile is never smaller than the signed one and c4's bound is never below c0's.** What c4 gives is a *tighter* bound at exact one-sided coverage: in the 70–180 bin its one-sided coverage is 0.893 against c0's lower side at 0.962, and in the 54–70 bin — where the point forecast sits above the reading — the tighter bound costs it: 0.793 against c0's 0.923 and c3's 0.910 at 30 minutes, 0.675 against 0.857 / 0.885 at 60, with 97 / 152 misses below against c0's 36 / 67. The symmetric methods protect the low readings by *over*-covering on the lower side everywhere else. Saturday's operating curve, where c4 appears beside c0–c3 and the point-forecast sweep, is where a tighter bound's fewer false alarms are weighed against those misses (`results/c4_hypoglycaemic_bins.csv`).

| True bin | h | n | c0: one-sided coverage · misses below · forecast − bound | c3 | c4 (signed score) |
|---|---|---|---|---|---|
| 54-70 | 30 | 468 | 0.923 · 36 · 32 mg/dL | 0.910 · 42 · 28 mg/dL | 0.793 · 97 · 21 mg/dL |
| 70-180 | 30 | 10,715 | 0.962 · 410 · 32 mg/dL | 0.949 · 545 · 24 mg/dL | 0.893 · 1151 · 19 mg/dL |
| 54-70 | 60 | 468 | 0.857 · 67 · 54 mg/dL | 0.885 · 54 · 40 mg/dL | 0.675 · 152 · 34 mg/dL |
| 70-180 | 60 | 10,715 | 0.962 · 412 · 54 mg/dL | 0.949 · 544 · 39 mg/dL | 0.889 · 1193 · 31 mg/dL |

*Conforming, 2020, 90% level; one-sided coverage is P(reading ≥ lower bound); "forecast − bound" is the median distance from the point forecast to the lower bound.*
