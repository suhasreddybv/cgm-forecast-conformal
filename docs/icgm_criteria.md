# iCGM special controls — 21 CFR 862.1355(b)

Transcribed 15 Oct 2026 by Suhas from the regulation text as reproduced at https://fda.innolitics.com/productcode/SAF (paragraph (b) of §862.1355). Verify against eCFR before citing in the memo: https://www.ecfr.gov/current/title-21/chapter-I/subchapter-H/part-862/subpart-B/section-862.1355

`src/eval/icgm.py` parses the tables in §2 and §3 of this file at run time; no threshold is typed into code (D-059).

## 1. What the criteria are, and the one fact about them that matters most

The accuracy requirements in (b)(1)(v) apply to an iCGM compared against **parallel blood samples measured on an FDA-accepted laboratory method traceable to a higher-order reference** ((b)(1)(ii)). Every criterion is stated **binned by the iCGM's own reading**, not by the reference value — the question the regulation asks is "when the device says X, how far can the truth be?", because X is what the user acts on. Every percentage must be reported as a **lower one-sided 95% confidence bound**, not a point estimate.

## 2. Accuracy criteria — (b)(1)(v), adult population

| Para | iCGM reading range | Tolerance | Lower one-sided 95% CI bound must exceed |
|---|---|---|---|
| (A) | < 70 mg/dL | within ±15 mg/dL of reference | 85% |
| (B) | 70–180 mg/dL | within ±15% of reference | 70% |
| (C) | > 180 mg/dL | within ±15% of reference | 80% |
| (D) | < 70 mg/dL | within ±40 mg/dL | 98% |
| (E) | 70–180 mg/dL | within ±40% | 99% |
| (F) | > 180 mg/dL | within ±40% | 99% |
| (G) | whole measuring range | within ±20% | 87% |

**Hard exclusions — zero tolerance, not a percentage:**

| Para | Rule |
|---|---|
| (H) | When the iCGM reads < 70, no reference value may be > 180 |
| (I) | When the iCGM reads > 180, no reference value may be < 70 |

## 3. Rate-of-change criteria — (b)(1)(v)(J)–(K)

| Para | Rule |
|---|---|
| (J) | No more than 1% of iCGM measurements may indicate a positive rate of change > +1 mg/dL/min when the true rate (from reference) is < −2 mg/dL/min |
| (K) | No more than 1% may indicate a negative rate of change < −1 mg/dL/min when the true rate is > +2 mg/dL/min |

These are the "trend arrow pointing the wrong way" rules.

## 4. Other controls relevant to this repo

- **(b)(1)(iv)** — consistent performance throughout the sensor wear period (performance drift).
- **(b)(1)(vii)** — no clinically significant gaps in data availability that would prevent a connected device from achieving its intended use. *This is the regulatory form of the gap analysis and the post-outage recovery result.*
- **(b)(4)** — acceptable performance in the presence of interfering substances.
- **(b)(7) labeling** must report accuracy in these ranges: **< 54, 54 to < 70, 70–180, > 180 to 250, > 250 mg/dL** ((7)(i)); positive and negative rate-of-change accuracy ((7)(ii)); frequency and duration of data gaps ((7)(iii)); **true, false, missed and correct alert rates** ((7)(iv)); observed sensor life ((7)(v)).

The labeling bins in (7)(i) are **exactly the true-glucose bins used in Tuesday's coverage breakdown.** The alert-rate quartet in (7)(iv) is **exactly the onset-detection table** Saturday produces: true alert rate = sensitivity, missed alert rate = 1 − sensitivity, false alert rate = false alarms, correct alert rate = PPV.

## 5. The caveat — why everything computed here is an analogue

These criteria govern a **sensor against a laboratory reference**. This repository compares a **forecast against a CGM** (Medtronic Enlite, adjunctive-labelled). The CGM is itself off from blood glucose by its own error. So:
- A forecast that *passes* an analogue has shown agreement with a sensor, not accuracy against blood; the real criterion is stricter by the sensor's own error.
- A forecast that *fails* an analogue fails something easier than the regulation asks.
- Nothing computed here is a claim that any device meets, or could meet, §862.1355. The word "analogue" appears on every table and figure.

## 6. Forecast analogues — definitions for Saturday (Claude Code implements from this section, not from memory)

Substitute **"forecast" for "iCGM measurement"** and **"CGM reading at the target time" for "reference blood glucose"**, and compute each criterion three times: on the **point forecast**, on the **interval lower bound**, and on the **interval upper bound** (a bound is an estimate too; whether the lower bound meets (A) and (H) is the hypoglycaemia question in regulatory form).

- **(A)–(G):** bin by the estimate's value, compute the percentage within tolerance of the CGM target, and report the **lower one-sided 95% bound** (Wilson). Note beside it that the targets are serially correlated, so the bound is optimistic; report the plain percentage too.
- **(H), (I):** count violations. Any violation is reported as a count and as a rate per 1,000 estimates — the regulation's requirement is zero.
- **(J), (K):** implied forecast rate = (estimate − origin reading) / horizon minutes; true rate = (target − origin reading) / horizon minutes. Count the wrong-direction cases as defined, as a percentage of all estimates; the requirement is ≤ 1%.
- **(7)(i) bins** are the breakdown bins already in the repo; reuse them.
- **(7)(iv) alert rates** come from the onset-detection table at each nominal level; label the four columns with the regulation's words.
- **(b)(1)(vii)** maps to the gap analysis and the post-outage recovery figure; cite both.

Expected and worth stating before computing: the model almost never forecasts below 70 (242 and 28 forecasts in 54–70), so (A) and (D) on the **point forecast** will rest on a tiny n and are not reportable under the n ≥ 100 rule; on the **lower bound** they will have a real n. That contrast is itself the finding.

Every output file and figure carries the word "analogue" in its name and title.
