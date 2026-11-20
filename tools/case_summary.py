# -*- coding: utf-8 -*-
"""One table for several HEC-RAS R04 cases (any runs/leh_2010_* folders), e.g.
    python tools/case_summary.py runs\\leh_2010_n004 runs\\leh_2010_burn2 runs\\leh_2010_n008
Reads each folder's hecras/runs/R04/readback/readback_summary.json, validation/validation_summary.json and the site's hecras model
(Manning n, channel burn), prints the table and writes runs/case_summary.csv."""
import csv
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def model_facts(d):
    out = {}
    g = glob.glob(os.path.join(d, "hecras", "model", "*.g0?"))
    if g:
        m = re.search(r"Storage Area Mannings=\s*([\d.]+)", open(g[0], errors="ignore").read())
        if m:
            out["mannings_n"] = float(m.group(1))
    log = os.path.join(d, "report", "run_log.txt")
    if os.path.exists(log):
        m = re.search(r"RAS_BURN_DEPTH_M\s*=\s*([\d.]+)", open(log, errors="ignore").read())
        if m:
            out["burn_m"] = float(m.group(1))
    return out


rows = []
for d in sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, "runs", "leh_2010_*"))):
    d = d if os.path.isabs(d) else os.path.join(HERE, d)
    r = dict(case=os.path.basename(d))
    r.update(model_facts(d))
    rb = load(os.path.join(d, "hecras", "runs", "R04", "readback", "readback_summary.json"))
    va = load(os.path.join(d, "validation", "validation_summary.json"))
    if rb:
        m = rb["model"]
        r.update(wet_km2=round(m["wet_area_km2"], 3), max_depth_m=round(m["max_depth_m"], 2), max_speed_ms=round(m["max_speed_ms"], 2),
                 extreme_km2=round(m["hazard_class_area_km2"]["extreme"], 3))
        for p in rb.get("probes", []):
            r.update(ST1_depth_m=round(p["depth_model_m"], 2), ST1_speed_ms=round(p["speed_model_ms"], 2))
    if va:
        for name, v in va["maps"].items():
            if v["kind"] == "depth":
                r.update(AUC=round(v["AUC"], 3), POD=round(v["POD"], 3), FAR=round(v["FAR"], 3), CSI=round(v["CSI"], 3),
                         points_hit_pct=round(v.get("points_hit_pct", float("nan")), 1))
    rows.append(r)
keys = []
for r in rows:
    keys += [k for k in r if k not in keys]
print("HEC-RAS R04 (1070 m3/s) cases vs the 2010 footprint; Sabu ST1 observed depth 5.2 m, speed 8.3 m/s (Thayyen et al. 2013)")
print("  ".join("%-14s" % k for k in keys))
for r in rows:
    print("  ".join("%-14s" % r.get(k, "-") for k in keys))
with open(os.path.join(HERE, "runs", "case_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, keys)
    w.writeheader()
    w.writerows(rows)
print("written runs/case_summary.csv")
