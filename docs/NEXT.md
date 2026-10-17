# Next — one question, January 2027

*Written 18 October 2026 at the close of Week 5. The repository is complete as published; this is the single item marked for return, after the application cycle.*

**The question.** Does training the forecaster with an asymmetric or quantile loss move hypoglycaemia onset sensitivity at two false alarms per day above the **27%** the point forecast achieves 30 minutes ahead (D-065)?

**Why this and nothing else.** Every post-hoc interval in Week 5 is centred on the frozen point forecast, and every one of them misses low readings in the same direction because the forecast sits above them (D-061). A signed score cannot reach below a symmetric one (D-062); conditioning on the rate of change widens where residuals are large, not where the low readings are (D-063). No residual-only score can move the interval's centre; only training can (D-069).

**What it unlocks.** A model trained with a pinball loss at a lower quantile, or an asymmetric loss that penalises over-prediction below 70 mg/dL, produces a lower quantile directly — which makes **conformalised quantile regression** available, the method excluded in Week 5 because the model is a point forecaster (D-046).

**What stays fixed.** The evaluation protocol (`docs/protocol.md`, hash `adf2584bff0271d3`), both evaluation variants (primary and conforming, D-027), the validation split (D-025), the clinical scorer (D-022–024), the onset-detection definition (D-047), the claim rule (D-064) and the existing operating curve as the test bed. The new model is scored on that curve, under both variants, at two false alarms per day, against the sweep — the same rule, the same bootstrap over six patients.

**What is not the target.** Pooled RMSE. A model that lowers RMSE and lowers onset sensitivity has repeated Week 4's result (D-039) and is reported as such.

Everything else in `docs/decisions.md` stays as logged.
