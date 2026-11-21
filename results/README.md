# Results

Outputs of our own runs, kept so you can see what fanrisk produces before you install it. Everything here was made by the commands in the main README; re-running them writes the same files into `runs/`.

## Leh 2010 (`leh_2010/`)

| File | What it shows |
| --- | --- |
| `viewer.html` | Interactive web map: AHP risk, HEC-RAS depth and hazard, the 2010 footprint, 1,958 buildings hit, critical facilities. Open in a browser. |
| `validation_summary.json` | Scores against the 2010 footprint: AHP risk AUC 0.85; HEC-RAS (2 m burn) AUC 0.77, POD 0.51, FAR 0.34, CSI 0.40 |
| `risk/ahp_C2_risk_abs.tif` (+ `.qml` style) | AHP risk map of the Sabu complex, 0-1; opens in QGIS with its colours |
| `risk/master_table.csv` | 29 morphometric and hydrological values for every feeder catchment |
| `hecras_R04/` | HEC-RAS run R04 (1,070 m3/s observed 2010 hydrograph, 2 m channel burn, n 0.06): max depth, max speed, hazard class, readback summary (Sabu ST1 depth 5.04 m vs 5.2 m observed) |
| `exposure/` | Buildings, roads and facilities hit (`exposure_report.txt` is the short version) |
| `sensitivity/` | Channel-burn and roughness experiments (`burn_summary.csv`, `case_summary.csv`); AHP weight sensitivity and the literature cross-check |
| `run_log.txt` | The full log of the run, with every setting that differs from the defaults |

## Dharali 2025 (`dharali_2025/`)

| File | What it shows |
| --- | --- |
| `event_check_automatic.txt` | Fully automatic run: no fan found (narrow trunk valley) |
| `event_check_manual_apex.txt` | With one manual apex: catchment 18.3 km2 (published 17), Melton 0.82 (0.80), debris flow |
| `viewer.html`, `risk/` | AHP risk map of the Kheer Gad fan |
| `run_log_*.txt` | Run logs |

## Kedarnath 2013 (`kedarnath_2013/`)

| File | What it shows |
| --- | --- |
| `event_check_automatic.txt` | Fully automatic run: the town lands on a 3 km2 side-stream fan; the trunk apex is missed |
| `event_check_manual_apex.txt` | With one manual apex on the Mandakini: catchment 18.6 km2, Melton 0.72, debris flow |
| `hecras_R04/`, `hecras_R05/` | HEC-RAS runs with the tool's own storm (109 m3/s clear water; 199 m3/s bulked): max depth, max speed, hazard class, `readback.log`, and `compare_map.png` / `hotspots.kml` comparing them with the AHP map |
| `hecras_batch_summary.csv` | Run times, volume error, wet area, AHP-vs-HEC-RAS AUC per run |
| `hecras_model_area.geojson` | The 2D flow area the tool chose |
| `viewer.html`, `risk/` | Web map and AHP risk map of the Mandakini fan complex |
| `run_log_*.txt` | Run logs |

## Figures (`figures/`)

| Figure | What it shows |
| --- | --- |
| `pipeline.png` | The pipeline, DEM to web map |
| `validation_2010_R04.png` | First scoring against the 2010 footprint (before the routing fix): HEC-RAS AUC 0.71, AHP 0.23 |
| `ahp_routing_fix_before_after.png` | Routing the flood on the real DEM instead of the channel-filled fan surface: AHP AUC 0.19 to 0.85 |
| `ahp_sensitivity_2010.png` | Monte Carlo on the AHP judgements and one-at-a-time weight changes |
| `channel_burn_sections.png`, `channel_burn_results.png` | Burning the channel 0, 2 or 4 m into the terrain; 2 m fits 2010 best |
| `exposure_map_leh2010.png` | Buildings and facilities hit, AHP vs HEC-RAS |
| `kedarnath_2013_R04.png` | Kedarnath: AHP risk, HEC-RAS hazard rating and their agreement (run R04) |

Background imagery in the figures and web maps: (c) Esri, Maxar, Earthstar Geographics.
