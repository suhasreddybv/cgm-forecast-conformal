# Evaluation protocol — frozen 5 October 2026, before any model was fitted

Nothing in this file may change after a model is fitted. If it has to change, the change is
a new dated section with the reason, not an edit, and every number produced under the old
protocol is re-run or struck.

Every rule below is taken from Marling & Bunescu (2020), *The OhioT1DM Dataset for Blood
Glucose Level Prediction: Update 2020*, or measured from the data and cited to a committed
artifact.

## The data

12 patients in two cohorts of six. The 2018 release (559, 563, 570, 575, 588, 591) wore the
Basis Peak band; the 2020 release (540, 544, 552, 567, 584, 596) wore the Empatica Embrace.
Eight weeks per patient, CGM every five minutes from a Medtronic Enlite sensor.

**Verified against the paper's Table 2 for all 12 patients** (`tests/test_loader.py`): every
training CGM count matches exactly, every 2018 test file matches its published count, and
every 2020 test file contains exactly 12 more readings — the first hour.

**CGM values are censored to [40, 400] mg/dL** by the sensor. This matters for the Week 5
hypoglycaemia work: a reading of 40 means "40 or below", so sensitivity at the 54 mg/dL
threshold cannot be evaluated against a sensor that does not resolve below 40.

## Splits

The provided training and testing files, per patient, used as the published work uses them.
**No random splits, no re-splitting, no pooling across patients for the split.** Training and
testing are disjoint in time and the test period follows the training period for every
patient (asserted in `tests/test_loader.py`).

**The start of the test period.** Quoting the paper: for the 2018 challenge "the number of
test examples was equal to the number of data points in the XML testing file", while for the
2020 challenge "the first hour of data in each XML testing file is excluded from the set of
points used for evaluation … as the first test points would otherwise be too close
chronologically to the training data."

So: **2020 cohort — skip the first 12 CGM readings of the test file. 2018 cohort — use all
of them.** The loader exposes this as `eval_start_index`. History for the earliest scored
test points is taken from the readings preceding them, including the end of the training
file where needed; no model may see a target before predicting it.

## Horizons

30 and 60 minutes: 6 and 12 steps ahead at the 5-minute cadence.

## Metrics

- **RMSE and MAE in mg/dL** — primary, and the metrics the published comparisons report.
- **MAPE** alongside, as in the PPG repository, because a 20 mg/dL error means something
  different at 70 and at 300 mg/dL.
- Reported per patient and aggregated, with the aggregation stated (mean of per-patient
  values, with the pooled figure beside it).
- **Error-grid (Clarke, Parkes), time-in-range agreement and hypoglycaemia sensitivity are
  Week 5 work and are not computed now.** They are listed as pending, not estimated.

## Gaps and targets

Measured today (`results/gap_analysis.csv`): the CGM series is **11.9% incomplete on average
in training and 10.7% in test**, the worst patient/split is 552/test at **40.2% missing**, and
the longest single gap is **118 hours**. Cut naively, **7.1% of 2-hour history windows cross a
gap** (13.9% for the worst patient).

**The rules, frozen:**

1. **A target is never interpolated.** A forecast is scored only against a real CGM reading.
   Any evaluation on an interpolated target inflates accuracy and is not permitted.
2. **History may be interpolated across gaps of at most 30 minutes** (6 readings), linearly,
   and a window containing a longer gap is dropped rather than repaired. 30 minutes is chosen
   because it is the shorter forecast horizon: a method may not lean on more imputed history
   than the distance it is being asked to predict.
3. **Every window records whether it was interpolated**, and results are reported with and
   without interpolated-history windows, so the cost of rule 2 is visible rather than assumed.
4. **Gap handling is applied identically to training and test.** No rule may make the test set
   easier than the training set.

## Fitting

Every scaler, imputer, hyperparameter and model parameter is fitted on training data only.

- **Patient-specific models**: fitted on that patient's training file, evaluated on that
  patient's test file.
- **Population models**: fitted on training patients only, with the evaluated patient held
  out, leave-one-patient-out.
- **Both families are reported.** Neither is dropped because the other looks better.

No hyperparameter is chosen on test data, including by looking at a test curve once.

## Comparison targets

The per-patient RMSE and MAE at 30 and 60 minutes from the dataset paper and from the Blood
Glucose Level Prediction Challenge results on this dataset. **Those numbers are extracted
from the sources in the Tuesday reading block and cited there — they are not written from
memory, and this file does not quote them until they are.**

## Stop conditions

- A result that beats the published figures by a wide margin is treated as a leak until
  proven otherwise, not as a finding.
- Any evaluation that scores an interpolated target is a bug, not a variant.
- If a model's advantage disappears when hyperparameters are chosen inside the fold rather
  than globally, it is reported as no gain — as in the PPG repository.
