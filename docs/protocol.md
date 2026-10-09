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

**CGM values are censored to [40, 400] mg/dL** by the sensor. Both hypoglycaemia thresholds
(70 and 54 mg/dL) sit *inside* that range, so sensitivity at either can be evaluated; what
censoring prevents is measuring error **below** 40. A true value of 30 is reported as 40, so a
forecast of 40 scores as exact when it is 10 mg/dL high.

Consequences, measured over all 166,533 readings: **206 readings (0.124%) sit exactly at the
40 mg/dL floor** and 335 (0.201%) at the 400 ceiling. The floor readings are not spread evenly
— **17.8% of all readings below 54 mg/dL sit exactly at the floor**, against 3.8% of those
below 70. So error-grid zones at the low end, and any metric computed on floored readings, are
**biased optimistic**, and that bias concentrates in the severe-hypoglycaemia region where
accuracy matters most. Errors below 40 are unmeasurable on this dataset. Any low-end result is
reported with the count of floored targets it rests on.

## The reference is a CGM, not a laboratory measurement

**The ground truth in this dataset is a Medtronic Enlite continuous glucose monitor, not a
laboratory reference.** Every result computed here — this repository's and every published
figure on OhioT1DM — measures agreement with a sensor that itself differs from the patient's
true blood glucose. A forecast that matched the CGM perfectly would still be wrong by the
sensor's own error against a laboratory assay.

This is not a caveat about precision at the margins. It bounds what any accuracy figure on
this dataset can mean: the numbers describe agreement with a device. **The Enlite was an
adjunctive device — not cleared for treatment decisions without a confirmatory fingerstick.**

That is specific to this sensor and this era, and should not be generalised to CGMs as a
class. Current sensors including the Dexcom G6 and G7 and the FreeStyle Libre 2 and 3 carry
non-adjunctive labelling, and the FDA's iCGM special controls (21 CFR 862.1355) set out what
a sensor must demonstrate to be dosed from directly. A result on OhioT1DM is measured against
a reference from the adjunctive generation.

No MARD figure for the Enlite is quoted here, because none has been sourced; if one is cited
later it comes with its reference, not from memory.

The practical rule: results are described as agreement with the CGM reference, never as
accuracy against blood glucose.

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
- **Every per-patient figure carries the number of real targets it was computed on, and every
  pooled figure states the total.** Patients are not comparable without it: 552's test split is
  40.2% missing with a 118-hour gap, so its per-patient error rests on far fewer evaluable
  targets than any other patient's. This is the same convention as the fold counts in the PPG
  repository, where S6 lacked three activities and the per-activity tables said so in every row.
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

## Covariate scope

Frozen before any fit, because the wearable channels are not uniformly available (D-004):
heart rate is present-but-empty for five 2020 patients and absent entirely for 596,
acceleration exists only in the 2020 cohort at a 1-minute cadence, and the 2018 band carries
channels the 2020 band does not.

- **Primary models use CGM, insulin (basal, temporary basal, bolus) and meals only.** These
  are present for all 12 patients, and they are what the published comparisons use, so the
  comparison stays like-for-like.
- **Wearable channels are a secondary, stratified analysis**, evaluated only on the patients
  where they exist and reported separately. They are never mixed into the primary comparison,
  and a result that depends on them is never presented as a result on the dataset.
- **"Missing" and "empty" are handled identically**, as established in D-004.

A model that quietly used heart rate would be evaluated on six patients and compared against
published numbers computed on twelve. That is the specific error this rule prevents.

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

**Each published figure is recorded with the inputs the method used.** The primary models here
are frozen at CGM, insulin and meals, so only methods using those same inputs are a like-for-
like comparison. A method that also consumed wearable channels, or that used only CGM, is a
different experiment and is labelled as such in the comparison table rather than ranked beside
ours. The table carries an `inputs` column and states, per row, whether the figure qualifies as
a like-for-like comparison.

## Stop conditions

- A result that beats the published figures by a wide margin is treated as a leak until
  proven otherwise, not as a finding.
- Any evaluation that scores an interpolated target is a bug, not a variant.
- If a model's advantage disappears when hyperparameters are chosen inside the fold rather
  than globally, it is reported as no gain — as in the PPG repository.

## Clarification, 9 October 2026 — validation hold-out for fitted models

Added before any sequence model is fitted, as the last edit to this file before a model
exists. It clarifies the "Fitting" section; it changes no rule above.

- **Early stopping, checkpoint selection and any hyperparameter choice use a temporal
  hold-out taken from the end of each training file: the last 20% of that patient's
  training windows, by time.** The first 80% fits; the hold-out is never fitted on and the
  test file is never looked at.
- **Patient-specific models** hold out that patient's training tail. **Population models**
  hold out every training patient's tail, pooled, with the evaluated patient's whole
  training file excluded as before (leave-one-patient-out).
- The split is by window index in time order, so the hold-out follows the fitting period
  and a window never straddles the split. Windows whose history reaches into the fitting
  period are allowed, as the first test windows draw history from the training tail.
- **Seeds are fixed and logged.** Every fitted model is run with at least three seeds and
  the per-seed figures are reported beside the mean; no seed is dropped.
- The 20% is a fixed fraction, not a tuned one. With roughly 9,500–11,600 training windows
  per patient it leaves 1,900–2,300 hold-out windows, enough to stop on.
