# Step 6: AHP weight sensitivity and literature cross-check (Leh 2010)

Figure: `data/validation/ahp_sensitivity_2010.png` (on Saharsh's PC, inside the fanrisk folder).

## What was tested
The base case is the fixed AHP map (`ROUTE_SURFACE_BLEND = 0`). Each test scores it two ways: against the traced 2010 footprint inside the survey window, and against HEC-RAS R04 wet cells (depth ≥ 0.3 m).

- **Monte Carlo on the judgements** (after Feizizadeh & Blaschke 2014). Each Saaty pairwise entry of both matrices is moved at random by ±1 or ±2 steps on the 1/9…9 ladder. The fan matrix is 6×6 and the catchment matrix is 7×7. Weights are re-solved, and samples with CR > 0.10 are rejected. The ±1 run kept 340 samples and the ±2 run kept 83 (most ±2 draws fail CR).
- **One-at-a-time (OAT)** (after Chen, Yu & Khan 2010). Each fan criterion weight is multiplied by 0.5 to 1.5 and the other weights are renormalised to sum to 1.

Baseline: CR fan 0.043, catchment 0.020.

| Matrix | Weights |
|---|---|
| Fan | HAND⁻¹ 0.306, upstream area 0.295, fan slope 0.212, local relief 0.100, event rainfall 0.087 |
| Catchment | slope 0.372, upstream 0.307, relief 0.249, ruggedness 0.073 |

## Results

### Monte Carlo AUC (5 / 50 / 95 %, minimum)

| Test | vs 2010 footprint | vs HEC-RAS R04 |
|---|---|---|
| Baseline | 0.849 | 0.812 |
| ±1 step (n=340) | 0.847 / 0.849 / 0.850, min 0.845 | 0.803 / 0.812 / 0.817, min 0.796 |
| ±2 steps (n=83) | 0.841 / 0.849 / 0.850, min 0.824 | 0.783 / 0.811 / 0.819, min 0.729 |

On cells the AHP map covers (91 % of the window), the AUC against the footprint is about 0.93 (0.928 to 0.936) and against HEC-RAS 0.78 to 0.82.

Fan weights under ±1 step (5–95 %):

| Criterion | Range |
|---|---|
| HAND⁻¹ | 0.23–0.35 |
| slope | 0.16–0.26 |
| upstream | 0.25–0.36 |
| relief | 0.08–0.14 |
| rain | 0.07–0.13 |

### One-at-a-time
- Against the footprint, every change stays between 0.926 and 0.936.
- Against HEC-RAS, HAND is the most influential criterion. Its weight at ×0.5 gives 0.760 and at ×0.7 gives 0.800.
- Fan slope is next: at ×1.5 it gives 0.794.
- Upstream area and relief are nearly flat.

### Stability of the "high" class
Of the cells rated high in the baseline:
- under ±1 step, 59 % stay high in at least 90 % of samples, and only 5 % of mapped cells flip class between 10 % and 90 % of the time;
- under ±2 steps, 39 % stay high and 9 % of mapped cells flip.

### Reading
The validation result does not depend on the exact judgements. Every accepted sample is above the 0.7 "acceptable" level, nearly all are in the 0.8–0.9 "excellent" band, and all are above the Ladakh benchmark of 0.766. The core of the high zone is stable; the edges are what move. HAND should be described as the criterion that carries the hydraulic agreement.

## Literature cross-check

| Our choice / result | Literature | Agreement | Source status |
|---|---|---|---|
| AUC 0.85 (0.93 on covered cells) | Hosmer, Lemeshow & Sturdivant (2013), p.177: 0.7–0.8 acceptable, 0.8–0.9 excellent, ≥0.9 outstanding | "Excellent" | verified |
| AUC vs a regional benchmark | Negi et al. 2020, Ladakh flash-flood susceptibility AUC 0.766 | Ours is higher, though our test is local and a different setup | from project papers |
| MC on Saaty judgements, CR ≤ 0.1 filter | Feizizadeh & Blaschke 2014, IJGIS 28(3):610–638 | Same method | verified |
| OAT weight scaling with renormalisation | Chen, Yu & Khan 2010, Env. Modelling & Software 25:1582–1591 | Same method | verified |
| CR ≤ 0.10 | Saaty (1980) | Baseline 0.043 / 0.020 | standard |
| Hazard classes 0.75 / 1.25 / 2.5, DF 0.5 | Defra/EA FD2321/TR1; FD2320 (DF 0.5 for d ≤ 0.25 m, 1.0 above) | Breaks match. DF 0.5 is the low end for the deep channel, so DF 1.0 would raise classes in d > 0.25 m cells | verified |
| Manning's n = 0.06 uniform | Chow (1959): mountain streams with cobbles and boulders, about 0.04–0.07 | Within range | from memory; check against the book before citing |
| Melton ratio 0.3 / 0.6 for flood / debris-flood / debris-flow | Wilford et al. 2004, Landslides 1:61–66 | Thresholds used as published | citation verified; thresholds not re-read |
| Footprint widths 104 / 546 / 650 m | Bhatt et al. 2011 damage mapping: about 130 / 400 / 650 m | Close; the middle reach is wider in ours | measured from their figure |
| Peak flows | Thayyen et al. 2013 | Hydrology within 2 % with published inputs (Assignment 3) | project papers |

## Caveats to write into the report
- The AUC against the footprint is flat partly because 9 % of the window has no AHP coverage, which is scored as 0. The covered-cell AUC (0.93) shows the ranking itself.
- AHP rates the first 1.2 km of the Sabu channel "moderate" where HEC-RAS says "extreme".
- DF = 0.5 everywhere understates the hazard class in deep cells (FD2320 uses 1.0 above 0.25 m).
