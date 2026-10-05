# Data: OhioT1DM

**This dataset is under a Data Use Agreement. Nothing from it, raw or derived, may be committed to this repository or redistributed in any form.** That includes processed CSVs, cached arrays, per-patient figures or tables that expose a glucose trace, model checkpoints trained on it, and any file keyed to the original timestamps. Aggregate metrics, error grids and distribution plots are permitted; a plotted time series is not.

## Access

A DUA must be executed before the data can be obtained. It is signed both by the principal investigator and by a legal signatory for the research institution — an individual researcher cannot execute it alone.

- Dataset page, maintained by Razvan Bunescu at UNC Charlotte: https://webpages.charlotte.edu/rbunescu/ohiot1dm.html
- The original Ohio University host (`smarthealth.cs.ohio.edu`) was taken offline in May 2024; the signed DUA form is emailed to the technology coordinator at Ohio University as directed on the page above.
- Cite: Marling, C., & Bunescu, R. (2020). The OhioT1DM dataset for blood glucose level prediction: Update 2020. *CEUR Workshop Proceedings*, 2675, 71–74.

## Layout the code expects

```
data/OhioT1DM/2018/{train,test}/{559,563,570,575,588,591}-ws-{training,testing}.xml
data/OhioT1DM/2020/{train,test}/{540,544,552,567,584,596}-ws-{training,testing}.xml
```

A symlink works, or set `OHIO_T1DM_ROOT` to wherever `OhioT1DM/` lives. `data/` is gitignored except this file.

The parsed cache is written **outside the repository** — `~/.cache/cgm-forecast-conformal` by default, `CGM_CACHE` to override — so that no cached array can be added to git even with `git add -f`.

## What is in it

12 patients in two cohorts of six, eight weeks each: CGM every five minutes (Medtronic Enlite), finger sticks, basal and temporary basal insulin, boluses, meals with carbohydrate estimates, self-reported sleep, work, stressors, hypoglycaemic events, illness and exercise, and sensor-band channels.

**The cohorts wore different bands and do not share a channel set.** 2018 (Basis Peak): heart rate, GSR, skin temperature, air temperature and steps, aggregated every 5 minutes. 2020 (Empatica Embrace): GSR, skin temperature and magnitude of acceleration, aggregated every 1 minute. See [docs/protocol.md](../docs/protocol.md) and `results/gap_analysis.csv` for everything verified against the data.

**CGM values are censored to [40, 400] mg/dL** by the sensor: a reading of 40 means "40 or below".
