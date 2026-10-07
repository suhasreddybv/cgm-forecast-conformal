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
  rules and stop conditions before a single baseline is fitted. **SHA-256 prefix `e1f76075a0decc0c`.**
  Any later change is a new dated section with a reason, never an edit, and results produced
  under the old protocol are re-run or struck.
- **Why now:** the PPG repository froze its evaluation design stage by stage and still had to
  correct a claim that came from a brief rather than a CSV (D-042 there). Writing the protocol
  down before the first fit removes the opportunity to choose the rule that flatters the result.
- **Rejected:** fitting a baseline "just to see" today. A baseline fitted before the protocol
  is frozen is evidence about the protocol, not about the baseline.

### D-003 · The test period starts where the published challenge says it does
- **Date / commit:** 2026-10-05 · `3e6f82f`
- **Status:** adopted — implemented in `src/data/windows.py` as D-014
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
- **Correction (2026-10-05, same day):** the first version of this entry said the reference is "the thing a clinician would not treat on alone". That is true of the Enlite and **not** of CGMs generally, and the overreach is the kind a reader at a sensor company would catch immediately. Precisely: **the Enlite was an adjunctive device, not cleared for treatment decisions without a confirmatory fingerstick.** Current sensors — Dexcom G6 and G7, FreeStyle Libre 2 and 3 — carry non-adjunctive labelling, and the FDA's iCGM special controls at **21 CFR 862.1355** define what a sensor must demonstrate to be dosed from directly. OhioT1DM's reference belongs to the adjunctive generation; that is the claim, and it is narrower than the one it replaces.

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

### D-012 · Published comparisons are recorded with the inputs each method used
- **Date / commit:** 2026-10-05 · `1564f78`
- **Status:** adopted
- **Decision:** the comparison table carries an `inputs` column and states, per row, whether a published figure is a like-for-like comparison with the primary models here.
- **Why:** the primary models are frozen at CGM, insulin and meals (D-010). A published method that also consumed wearable channels, or that used CGM alone, ran a different experiment. Ranking those figures in one column against ours would compare methods and input sets at the same time, and the difference would be invisible.
- **Consequence for Tuesday's reading block:** for each published method, the inputs are recorded alongside the RMSE and MAE, and only the matching ones are described as comparable. The numbers themselves are still extracted from the sources then, not from memory.

### D-013 · Windows are anchored on the target and valued by interpolation in time
- **Date / commit:** 2026-10-06 · `6d4f969`
- **Status:** adopted
- **Decision:** one window per real CGM reading that has enough history (`src/data/windows.py`). For a target at time T and horizon h, the history is H slots at T − (h+k)·300 s, k = H−1 … 0, so the history end sits exactly h steps before the target. The value at a slot is the reading at that instant if one exists, else the linear interpolation of the two real readings bracketing it. A window is dropped if any bracketing interval exceeds 1,800 s (protocol rule 2) or a slot precedes the first reading available. Built for **h ∈ {6, 12} and H ∈ {6, 12, 24}** — 144 window sets, cached under `~/.cache/cgm-forecast-conformal/windows`, outside the tree (D-001).
- **Why interpolate in time rather than take H consecutive readings:** the cadence is not exactly 300 s. **Measured over all 24 files: 165,482 intervals are exactly 300 s, 503 are 301–360 s, 3 are shorter**, and of the 521 intervals over 360 s, 125 are not multiples of 300 s. A window of H consecutive readings silently stretches across any gap inside it, and a global 5-minute grid loses phase after every jittered interval. Anchoring on the target and interpolating handles both with one rule; the value is exact whenever a slot coincides with a reading, and a 301 s interval moves it by at most 1/301 of one step.
- **What counts as interpolated:** every window records the widest bracketing interval it rests on (`max_bracket_s`). It is reported as interpolated across a gap when that exceeds 360 s — the definition `gap_analysis.py` already uses — so one-second jitter is not reported as a gap. Protocol rule 3 ("every window records whether it was interpolated") is met with the number, not just the flag.
- **Pinned by `tests/test_windows.py`:** the closed form for an isolated gap of g missing slots — **H−1+min(g, h) targets touched; all dropped if the gap exceeds 30 minutes, all kept and flagged otherwise** — on synthetic series across both horizons and three history lengths; the 30-minute boundary is inclusive; a 301 s interval is kept unflagged; a 783 s interval is kept and flagged; targets are real readings with history ending exactly h steps before them. On the real data, for every patient × split × h × H: no target timestamp is absent from the raw stream, no history span exceeds 30 minutes, and no history slot is later than its target or outside the patient's own series. **Verified against the closed form on 559/train: 36 gaps over 30 minutes × (H−1+h) = 396 targets dropped at h=6, H=6, which is what the windower produced.**
- **Rejected:** index-space windows of H consecutive readings — the common shortcut, and the one whose blind spot D-016 measures. A global 5-minute grid anchored at the first reading, because the phase shifts after every non-multiple gap and every jittered interval.

### D-014 · The test warm-up as implemented
- **Date / commit:** 2026-10-06 · `6d4f969`
- **Status:** adopted (implements D-003)
- **Decision:** for the test split the history series is the training file followed by the test file, and candidate targets start at `eval_start_index` — the first 12 readings of a 2020 test file serve as history and are never scored (D-003); every 2018 test reading is a candidate. History never reaches past the target and never outside the patient's own two files.
- **Measured:** the training tail ends **exactly 300 s before the first test reading for all 12 patients**, so the test file continues the training file with no gap and early test targets draw full history from the training tail. **`n_dropped_no_history` is 0 for every test split** at every h and H (`results/effective_n.csv`). Both facts are pinned by `tests/test_windows.py::test_test_split_warm_up_as_the_protocol_specifies`.
- **The one case where the first candidate is not the first target:** 552's training file has a **409-minute gap ending 10 readings before the split**. At H=24 the first scored test points reach into it and are dropped — **7 at h=6, 13 at h=12; 1 at h=12, H=12; none at H=6** — under the same 30-minute rule as everywhere else. No other patient drops any candidate before its first target. The test asserts exactly this: every candidate before the first target was dropped for a long gap, never for lack of readings, and only on 552.

### D-015 · Covariates are aligned to the history slots, causally
- **Date / commit:** 2026-10-06 · `6d4f969`
- **Status:** adopted
- **Decision:** the channels the primary models may use (D-010) are carried on every window as separate `[n, H]` arrays aligned to the same slots as the glucose history: **`basal`** — U/h in effect at the slot instant, with a temporary basal overriding the schedule while active (begin inclusive, end exclusive; where two overlap, the later-begun one); **`bolus`** — units delivered in the five minutes ending at the slot, with an extended bolus spread uniformly over its delivery interval; **`carbs`** — grams logged in the five minutes ending at the slot. Nothing after the history end enters any of them: a bolus or meal between the history end and the target is unseen.
- **Measured, and why each rule exists:** of **3,733 boluses, 245 are extended** (227 `square dual`, 18 `square`), delivered over up to 3,600 s — dumping them at `ts_begin` would place an hour of insulin in one five-minute bin. Of **572 temporary basals, one pair overlaps** (591/train). **575's first basal event comes 56 minutes after its first CGM reading**, so 12 training windows per (h, H) carry NaN basal for the slots before it (57–78 slots); no other file has unknown basal, because the training file's last scheduled rate carries into the test split. NaN is kept rather than filled: a rate that was not recorded is not zero.
- **Pinned by `tests/test_windows.py`:** a step basal with a temp-basal override, bin boundaries left-open right-closed, a normal bolus landing in its bin and a 30-minute square bolus spreading 0.5 U per full bin, on synthetic channels.
- **Not done, on purpose:** no feature engineering (insulin on board, time since meal, time of day) — that is the model's business on Saturday, and building it today would be a design choice made before the baselines exist. Insulin planned inside the horizon is excluded as future information; if it is ever used it is a separate, labelled analysis.

### D-016 · The effective-n table, and the gap analysis reconciled with the windower
- **Date / commit:** 2026-10-06 · `6d4f969`
- **Status:** adopted
- **Deliverable:** `results/effective_n.csv` — per patient × split × horizon × H, the number of evaluable real targets, split into clean and interpolated-history, with what was dropped and why, as a fraction of the file's readings and of the readings the sensor should have produced. **This is the n column that accompanies every result from here on (D-009).** Aggregate counts only.
- **Result:** mean evaluable share of raw readings is **97.4% (train) / 97.3% (test) at h=6, H=6**, falling to **92.1% / 92.3% at h=12, H=24**; worst files 567/train at 84.5% and 584/test at 84.9%. **552/test keeps 88.8% of its own readings at h=12, H=24 — but that is 2,100 targets against the 3,950 its sensor should have produced (53.2%).** The cost of a long gap is paid in absolute n, not in the fraction: the readings inside 552's 118-hour gap never existed, so they are not targets to lose. The table carries both denominators so neither reading of the number can hide the other.
- **Reconciliation** (`results/window_reconciliation.csv`, pinned by `tests/test_windows.py::test_reconciles_with_the_gap_analysis`): the naive count in D-005 — "1.6% of 30-minute history windows cross a gap" — counts windows of H consecutive *readings*. Recomputed over the same candidate targets it gives **2,601 crossing windows at h=6, H=6**; the windower finds **2,578 targets touched on the history side**, within 0.9%. At H=24: 11,813 against 11,727 (0.7%); at h=12 the shortfall is 1.2–1.5%. The per-file shortfall is **never positive and at most 27 targets**, and comes from clustered gaps, which the two counts attribute to different gaps.
- **What the naive count cannot see:** a further **2,705 targets at h=6 and 4,523 at h=12 have their history *end* inside a gap** — the target lies within h steps after a gap, so in index space its "latest" reading is hours old and the window never notices. Those are real losses under the protocol and the index-space count excludes them by construction. **The true cost of gap handling at h=6, H=6 is 5,283 targets, 2.0× the naive figure**; at h=12, H=24 it is 16,181 against 11,815. The isolated-gap closed form is H−1+min(g, h) per gap, of which the min(g, h) term is this horizon-side effect.
- **Stop condition, considered:** the brief said a reconciliation off by more than a few percent means one of the two is wrong. The totals differ by a factor of two, but the component the two methods both measure agrees to within 1.5% and the remainder is located exactly, target by target, and pinned by a closed form on synthetic data. Neither implementation is wrong; the naive diagnostic measured a narrower quantity than its sentence implied. Reported here rather than stopped. **The protocol is unchanged** — the rules are the same; only the diagnostic count in D-005 was narrower than it read, and `results/gap_analysis.csv` keeps its original meaning.
- **Not run:** no persistence baseline, although it needs no fitting. Baselines land Wednesday beside the published numbers.

### D-017 · The scored set starts where the published rule says, whatever history is available
- **Date / commit:** 2026-10-07 · `pending`
- **Status:** adopted (verifies D-003 and D-014)
- **The rule, from the protocol:** 2020 cohort — the first hour of each test file "is excluded from the set of points used for evaluation" (Marling & Bunescu 2020), so the first 12 readings serve as warm-up and are never scored; 2018 cohort — "the number of test examples was equal to the number of data points in the XML testing file", so every test reading is a candidate. D-014 made history available for the earliest test targets from the training tail; that may not move the start of scoring.
- **Evidence, verified before any fit:** across all **72 test window sets** (12 patients × 2 horizons × 3 history lengths), **no scored target precedes the published start of the test period** — zero violations. Scored candidates equal the challenge counts exactly (2018: the file's reading count; 2020: the count minus 12, e.g. 540: 2,884 of 2,896). Pinned by `tests/test_windows.py::test_test_split_warm_up_as_the_protocol_specifies`, which now asserts `target_ts.min() >= cgm.ts[eval_start_index]` for every set.
- **Consequence:** `results/effective_n.csv` is unchanged, and its `evaluable_targets` column is the scored n under every baseline figure (`tests/test_baselines.py::test_scored_n_equals_the_effective_n_table`).

### D-018 · Three baselines, and why each exists
- **Date / commit:** 2026-10-07 · `pending`
- **Status:** adopted
- **Decision** (`src/eval/baselines.py`, `results/baselines.csv`, `results/baselines_per_patient.csv`): **p0 persistence** — the value at the last history slot, T − h, the same most-recent reading every later method sees, so it is the floor they are measured from; the slot value is within 1/301 of a step of the raw reading (D-013). **l1 linear extrapolation** — a least-squares line through the last k slots, k ∈ {2, 3, 6}, extended h steps; it measures how much of the task is local trend. **l2 linear autoregression** — OLS on the H slots, H ∈ {6, 12, 24}, with and without the allowed covariates (D-010, D-015). All evaluated on the test split only, scored targets only, 30 and 60 minutes.
- **Result, pooled mean of per-patient RMSE (pooled figure in brackets), n under each:** p0 **23.39 (23.55) at 30 min, n = 30,912; 38.54 (38.79) at 60 min, n = 30,579**; MAE 17.01 / 28.70, MAPE 11.4% / 19.4%. Worst patient 540 at 28.45 / 47.52, best 570 at 19.20 at 30 min. l1 is worse than p0 at every k: 29.86 / 63.35 (k=2), 28.23 / 59.13 (k=3), 28.15 / 55.86 (k=6). l2 best is **per-patient, H=24, with covariates: 19.02 (19.17) / 31.97 (32.18)** — 4.4 and 6.6 mg/dL below persistence; without covariates 19.94 / 33.78; H=12 without covariates is identical at 30 min (19.94), so history beyond an hour adds nothing to a linear model without covariates. Per cohort, 2020 is harder: p0 24.16 against 22.61 at 30 min, 40.27 against 36.81 at 60.
- **Reported, not selected:** every H and every k. No row is labelled best in the results; the README names the best configuration in prose only, and nothing was chosen on the test split.
- **Published comparison column:** present and blank in `results/baselines.csv` (`published_rmse`, `published_mae`, `published_source`, `published_inputs`, `like_for_like`), pinned blank by a test until filled from the reading table under D-012.

### D-019 · The fitting rule as applied
- **Date / commit:** 2026-10-07 · `pending`
- **Status:** adopted
- **Decision:** per-patient l2 is fitted on that patient's training windows; the population l2 is fitted **leave-one-patient-out on the other eleven patients' training windows**, the evaluated patient held out — the protocol's definition ("fitted on training patients only, with the evaluated patient held out"). The day's brief described the population fit as "on all twelve training files"; the frozen protocol takes precedence and the pooled-twelve variant was not run. Coefficients come from `np.linalg.lstsq` on training windows only; **no scaler**, because OLS is invariant to feature scaling and a scaler would be one more fitted object to keep clean. Nothing is chosen on the test split.
- **Result:** the LOPO population fit is within **0.2 mg/dL of the per-patient fit at 30 minutes** (20.10 against 19.94 at H=24 without covariates; 19.62 against 19.02 with) and **0.3–0.5 worse at 60** (34.09 against 33.78; 33.32 against 31.97). For a linear model the patients are nearly interchangeable; whether that holds for the sequence model is Saturday's question.
- **Pinned:** `tests/test_baselines.py` recovers a known linear generator exactly from training windows, checks the per-step slope of l1 across horizons, and checks the design widens to 1 + 4H columns with covariates.

### D-020 · Windows with a NaN covariate are dropped from the with-covariate rows only
- **Date / commit:** 2026-10-07 · `pending`
- **Status:** adopted
- **Decision:** a window whose basal, bolus or carbohydrate array contains NaN is dropped from any fit or score that uses covariates, and kept in the rows that do not, so the two are compared on stated n rather than on identical n. Nothing is imputed: a rate that was not recorded is not zero (D-015).
- **Measured:** the only NaN covariate in the dataset is 575's basal before its first scheduled-rate event, 56 minutes into its training file (D-015). **12 training windows per (h, H) are dropped from 575's with-covariate fits and from every LOPO fit that includes 575; zero test windows are dropped for any patient.** So every with-covariate test figure in `results/baselines.csv` is on exactly the same n as its without-covariate row, and only the training n differs, by 12. Pinned by `tests/test_baselines.py::test_nan_covariate_windows_exist_only_in_575s_training_file`.

### D-021 · Sanity checks before the baseline numbers are trusted
- **Date / commit:** 2026-10-07 · `pending`
- **Status:** adopted
- **Result — the three checks, run by `src/eval/baselines.py` after writing the tables (the script exits non-zero if any fails) and re-run on the committed CSVs by `tests/test_baselines.py`:**
  1. Persistence RMSE at 30 minutes within [10, 40] mg/dL for every patient — **pass, range 19.20–28.45**. Outside that range the windows or the target alignment would be suspect before the data.
  2. 60-minute RMSE exceeds 30-minute RMSE for every baseline and every patient — **pass, 216 comparisons, 0 violations**.
  3. Two-point extrapolation an hour out is noisier than persistence for most patients — **pass, 12 of 12**; a uniform improvement would have meant the slope was in the wrong units.
- **Not a check, but noted:** the persistence floor is 23.4 / 38.5 mg/dL. Whatever Saturday's model reports is read against that line first and against the linear AR second; a sequence model that does not clear 19.0 / 32.0 has not earned its parameters.
