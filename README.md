# cgm-forecast-conformal

Glucose forecasting at 30 and 60 minutes with calibrated conformal prediction intervals, evaluated by clinical error grids rather than RMSE alone.

**Result.** _Pending. Headline number goes here once baselines and the sequence model are evaluated on held-out patients._

![Clarke and Parkes error grids, 30- and 60-minute horizons](figures/hero.png)

## Method

_To be written (5–18 Oct 2026)._ Planned: persistence and linear baselines reported first; sequence model second; split-conformal intervals with empirical vs nominal coverage; MARD, Clarke and Parkes zones, time-in-range agreement, hypoglycaemia sensitivity at 70 and 54 mg/dL. Error-grid figure is the hero image.

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# commands added as the pipeline lands
```

Data: see [data/README.md](data/README.md). No data is included in this repository.

## Limitations and failure cases

_To be written. Must include horizon-dependent coverage failures, hypoglycaemia misses and per-patient variation._

## References

- Marling, C., & Bunescu, R. (2020). The OhioT1DM dataset for blood glucose level prediction: Update 2020. *CEUR Workshop Proceedings*, 2675, 71–74.
- Clarke, W. L., et al. (1987). Evaluating clinical accuracy of systems for self-monitoring of blood glucose. *Diabetes Care*, 10(5), 622–628.
- Parkes, J. L., et al. (2000). A new consensus error grid to evaluate the clinical significance of inaccuracies in the measurement of blood glucose. *Diabetes Care*, 23(8), 1143–1148.

## Licence

MIT
