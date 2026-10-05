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
- **Date / commit:** 2026-10-05 · `pending`
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
- **Date / commit:** 2026-10-05 · `pending`
- **Status:** adopted
- **Decision:** `docs/protocol.md` fixes the splits, horizons, metrics, gap rules, fitting
  rules and stop conditions before a single baseline is fitted. **SHA-256 prefix `361352fd45a09072`.**
  Any later change is a new dated section with a reason, never an edit, and results produced
  under the old protocol are re-run or struck.
- **Why now:** the PPG repository froze its evaluation design stage by stage and still had to
  correct a claim that came from a brief rather than a CSV (D-042 there). Writing the protocol
  down before the first fit removes the opportunity to choose the rule that flatters the result.
- **Rejected:** fitting a baseline "just to see" today. A baseline fitted before the protocol
  is frozen is evidence about the protocol, not about the baseline.

### D-003 · The test period starts where the published challenge says it does
- **Date / commit:** 2026-10-05 · `pending`
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
- **Date / commit:** 2026-10-05 · `pending`
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
- **Date / commit:** 2026-10-05 · `pending`
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
- **Date / commit:** 2026-10-05 · `pending`
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
