# -*- coding: utf-8 -*-
"""Step 7 summary: one table comparing HEC-RAS run R04 on the plain DEM and with the channel burned in (run from the fanrisk folder).
    python tools/burn_summary.py
Reads runs/leh_2010_burn<D>/hecras/runs/R04/readback/readback_summary.json and runs/leh_2010_burn<D>/validation/validation_summary.json,
prints the table and writes runs/burn_summary.csv."""
import csv
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


rows = []
for d in sorted(glob.glob(os.path.join(HERE, "runs", "leh_2010_burn*")), key=lambda p: float(re.findall(r"burn([\d.]+)", p)[0])):
    depth = re.findall(r"burn([\d.]+)", d)[0]
    rb = load(os.path.join(d, "hecras", "runs", "R04", "readback", "readback_summary.json"))
    va = load(os.path.join(d, "validation", "validation_summary.json"))
    r = dict(burn_m=depth)
    if rb:
        m = rb["model"]
        r.update(wet_km2=round(m["wet_area_km2"], 3), max_depth_m=round(m["max_depth_m"], 2), max_speed_ms=round(m["max_speed_ms"], 2),
                 extreme_km2=round(m["hazard_class_area_km2"]["extreme"], 3))
        for p in rb.get("probes", []):
            r.update(ST1_depth_m="%.2f (obs %.1f)" % (p["depth_model_m"], p["depth_obs_m"]), ST1_speed_ms="%.2f (obs %.1f)" % (p["speed_model_ms"], p["speed_obs_ms"]))
    if va:
        for name, v in va["maps"].items():
            if v["kind"] == "depth":
                r.update(AUC=round(v["AUC"], 3), POD=round(v["POD"], 3), FAR=round(v["FAR"], 3), CSI=round(v["CSI"], 3),
                         points_hit_pct=round(v.get("points_hit_pct", float("nan")), 1), holdout_hit_pct=round(v.get("holdout_points_hit_pct", float("nan")), 1))
    rows.append(r)

if not rows:
    raise SystemExit("no runs/leh_2010_burn* folders yet: run run_burn_test.bat first")
keys = []
for r in rows:
    keys += [k for k in r if k not in keys]
print("Step 7 - HEC-RAS R04 (1070 m3/s) on plain vs channel-burned terrain, scored against the 2010 footprint")
print("  ".join("%-16s" % k for k in keys))
for r in rows:
    print("  ".join("%-16s" % r.get(k, "-") for k in keys))
with open(os.path.join(HERE, "runs", "burn_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, keys)
    w.writeheader()
    w.writerows(rows)
print("written runs/burn_summary.csv")
