# Design decisions

An append-only record of the design choices in this repository, the measured evidence behind
each, and what was rejected. Entries are dated by when the decision landed in git, with the
commit cited.

Every number here is copied from a committed results file, test or script output, and the
source is named. Nothing is quoted from memory. Entries are never rewritten: a changed
decision gets a new entry, and the old entry's **Status** line alone is updated to point at
it. Corrections stay visible beside what they correct.

Same format and rules as the decision log in `wrist-ppg-motion-robust-hr`.

---

### D-001 · The data use agreement governs everything in this repository
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted
- **Decision:** nothing derived from OhioT1DM enters git — not raw XML, not caches, not
  intermediate arrays, not per-patient statistics keyed to timestamps, and no figure that
  reproduces any patient's glucose trace. Aggregate metrics, error grids and distribution
  plots are permitted; a plotted time series is not.
- **Evidence that the protections are active, not merely configured:** the pre-commit hooks
  were tested by attempting four commits. A random (non-allowlisted) AWS-style credential was
  **blocked** by gitleaks; a private-key header was **blocked** by two hooks; a 6 MB file was
  **blocked** by the large-file hook; and a file under `data/` **could not be staged at all**
  because the directory is gitignored. An earlier attempt using AWS's *documented* example key
  passed, because gitleaks allowlists it by design — that test was inconclusive, not a gap, and
  was repeated with a random key.
- **Cache:** outside the tracked tree entirely, at `~/.cache/cgm-forecast-conformal`
  (`CGM_CACHE` to override), so no cache file can be added even by accident.
- **Rejected:** keeping a cache under `data/` and relying on `.gitignore`. One `git add -f`
  would defeat it.

### D-002 · The evaluation protocol is frozen before any model is fitted
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted
- **Decision:** `docs/protocol.md` fixes the splits, horizons, metrics, gap rules, fitting
  rules and stop conditions before a single baseline is fitted. **SHA-256 prefix `34bc943a95d39159`.**
  Any later change is a new dated section with a reason, never an edit, and results produced
  under the old protocol are re-run or struck.
- **Why now:** the PPG repository froze its evaluation design stage by stage and still had to
  correct a claim that came from a brief rather than a CSV (D-042 there). Writing the protocol
  down before the first fit removes the opportunity to choose the rule that flatters the result.
- **Rejected:** fitting a baseline "just to see" today. A baseline fitted before the protocol
  is frozen is evidence about the protocol, not about the baseline.

### D-003 · The test period starts where the published challenge says it does
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted
- **Decision:** for the 2020 cohort the **first 12 CGM readings of each test file are excluded
  from evaluation**; for the 2018 cohort every test reading is scored. Exposed as
  `eval_start_index` on the record.
- **Evidence:** the paper states that for the 2020 challenge "the first hour of data in each
  XML testing file is excluded from the set of points used for evaluation", and that the files
  therefore hold 12 more points than Table 2's test counts. **Measured on the data: every one
  of the six 2020 test files is exactly +12 readings against Table 2, and every 2018 test file
  is exactly +0** (`tests/test_loader.py::test_cgm_counts_match_the_published_table`).
- **Why it matters:** scoring those 12 points would make the first predictions trivially easy —
  they sit minutes after the training data ends — and would silently inflate every 2020 result
  relative to published work.

### D-004 · The two cohorts do not share a channel set
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted
- **Decision:** band channels are validated per cohort, and any model using them must either
  restrict itself to the intersection (GSR, skin temperature) or be reported per cohort.
- **Measured** (all 12 patients, training files): the 2018 cohort carries heart rate, GSR, skin
  temperature, air temperature and steps at a 5-minute aggregation and has **no `acceleration`
  element at all**. The 2020 cohort carries GSR, skin temperature and acceleration at a
  1-minute aggregation — roughly 28,000–53,000 rows against ~11,000 CGM readings — and its
  `basis_heart_rate`, `basis_air_temperature` and `basis_steps` elements are **present but
  empty** for five patients and **absent entirely** for 596.
- **Consequence for code:** "channel missing" and "channel present but empty" mean the same
  thing and must be handled identically. The loader does so and asserts the cohort difference
  in both directions — a 2020 patient carrying Basis-only data would mean the cohorts had been
  mixed up, and that raises.
- **Also measured:** `meal` is **not** universal — 567's test file contains none — and
  `temp_basal` is empty in three test files. Both are optional; requiring them would reject
  valid patients. Only glucose, finger sticks, basal and bolus are present everywhere.

### D-005 · Targets are never interpolated
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted
- **Decision:** a forecast is scored only against a real CGM reading. History may be
  interpolated across gaps of **at most 30 minutes**; a window containing a longer gap is
  dropped rather than repaired; every window records whether it was interpolated, and results
  are reported with and without those windows.
- **Evidence** (`results/gap_analysis.csv`): CGM is **11.9% incomplete on average in training
  and 10.7% in test**. The worst patient/split is **552/test at 40.2% missing**, and the longest
  single gap in the dataset is **118 hours**. Cut naively, **1.6% of 30-minute history windows,
  3.4% of 1-hour and 7.1% of 2-hour windows cross a gap** — 13.9% for the worst patient.
- **Why 30 minutes:** it is the shorter forecast horizon. A method may not lean on more imputed
  history than the distance it is being asked to predict.
- **Rejected:** interpolating targets, which would let a model be scored against a number no
  sensor produced, and would inflate accuracy by exactly the amount the imputation is smooth.
  Also rejected: dropping every window touched by any gap, which would discard a sixth of the
  data and bias the evaluation toward the best-instrumented stretches.

### D-006 · Sub-nominal CGM intervals are counted, not failed and not ignored
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted
- **Decision:** intervals shorter than the nominal five minutes are recorded on the record as
  `short_intervals` and pinned by a test; duplicate or out-of-order timestamps still raise.
- **Evidence:** across all 24 files there are exactly **three** sub-nominal intervals, the
  shortest **179 s** (patient 540, training), and **no duplicate timestamps anywhere**. Twenty-
  three of 24 files have a minimum interval of exactly 300 s.
- **Why not fail:** rejecting a patient over one reading in 141,000 would lose eight weeks of
  data for a rounding artefact. **Why not ignore:** an undocumented cadence violation is how a
  resampling bug hides. The first draft of the loader did fail, and the count was measured
  before the rule was relaxed.

### D-007 · The reference is a CGM, not a laboratory measurement
- **Date / commit:** 2026-10-05 · `25c9809`
- **Status:** adopted
- **Decision:** every result in this repository is described as **agreement with the CGM reference**, never as accuracy against blood glucose.
- **Why it matters more than any other limitation here:** the ground truth is a Medtronic Enlite sensor. Every figure computed on OhioT1DM — this repository's and every published one — inherits that sensor's own error against a laboratory assay. A model that matched the CGM perfectly would still differ from the patient's true glucose by the sensor's error, and that error is not estimable from this dataset because the dataset contains no laboratory reference to compare against.
- **Deliberately not quoted:** no MARD figure for the Enlite appears anywhere in this repository. None has been sourced, and a number carried from memory is exactly the failure the PPG repository's D-042 was about. If one is cited later it arrives with its reference.
- **Consequence:** a reader deciding whether this could inform treatment needs the sensor's accuracy as well as the model's, and only the second is measurable here.

### D-008 · Censoring at the sensor floor biases the low end optimistically
- **Date / commit:** 2026-10-05 · `25c9809`
- **Status:** adopted (corrects an earlier statement)
- **Correction:** an earlier note in this repository said hypoglycaemia sensitivity at 54 mg/dL "cannot be evaluated" because of censoring. **That was wrong.** Both 70 and 54 mg/dL sit inside the sensor's [40, 400] reporting range and both can be evaluated.
- **What censoring actually does:** it prevents measuring error *below* 40. A true value of 30 is reported as 40, so a forecast of 40 scores as exact when it is 10 mg/dL high. Errors below 40 are unmeasurable on this dataset.
- **Measured** over all 166,533 readings: **206 (0.124%) sit exactly at the 40 mg/dL floor** and 335 (0.201%) at the 400 ceiling. The floor readings concentrate where it matters: **17.8% of all readings below 54 mg/dL are at the floor**, against 3.8% of those below 70.
- **Rule:** error-grid zones at the low end, and any metric computed on floored readings, are reported as **biased optimistic**, with the count of floored targets they rest on.

### D-009 · Per-patient figures carry their effective n
- **Date / commit:** 2026-10-05 · `25c9809`
- **Status:** adopted
- **Decision:** every per-patient result states the number of real targets it was computed on, and every pooled figure states the total.
- **Evidence:** `results/gap_analysis.csv`. 552's test split is **40.2% missing** with a **118-hour** gap; the best-instrumented splits are above 95% complete. Comparing a per-patient RMSE from 552 against one from 588 without that context compares two different quantities.
- **Precedent:** the fold counts carried in every row of the PPG repository's per-activity tables, which existed because S6 lacked three activities.

### D-010 · Covariate scope is frozen before the first fit
- **Date / commit:** 2026-10-05 · `25c9809`
- **Status:** adopted
- **Decision:** **primary models use CGM, insulin (basal, temporary basal, bolus) and meals only.** Wearable channels are a secondary, stratified analysis, evaluated only on patients who have them and reported separately, never mixed into the primary comparison.
- **Evidence:** D-004. Heart rate is present-but-empty for five 2020 patients and absent entirely for 596; acceleration exists only in the 2020 cohort at a 1-minute cadence; the 2018 band carries channels the 2020 band lacks. Only glucose, finger sticks, basal and bolus are non-empty in all 24 files.
- **The specific error this prevents:** a model that quietly used heart rate would be evaluated on six patients and then compared against published numbers computed on twelve. The comparison would be meaningless and would look like a result.
- **Also:** "missing" and "empty" channels are handled identically, as already established in D-004.

### D-011 · The decision log is guarded by a test, from day one
- **Date / commit:** 2026-10-05 · `25c9809`
- **Status:** adopted
- **Decision:** `tests/test_decision_log.py` parses this file and asserts its own rules: every entry has an id and a status, ids are unique and sequential, every status reference points at an entry that exists, an entry that supersedes or closes another leaves that predecessor's status updated, no entry is both adopted and superseded, nothing marked open is claimed resolved elsewhere, and every adopted entry cites a committed artifact or a measured number. It runs in CI.
- **Why it exists:** in the PPG repository five entries were left marked `open` after their successors landed, and the inconsistency survived a full close-out review. Prose review does not catch this class of error; parsing does.
- **Verified by reintroducing the bug:** marking a resolved predecessor `open` again makes two independent checks fail and name both entries. The guard was not written to pass on a log that already happened to be clean.
- **Rejected:** retrofitting the guard at the end of the repository, which is when the inconsistency has already had weeks to accumulate.
