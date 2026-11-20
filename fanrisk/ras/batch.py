# -*- coding: utf-8 -*-
"""Run the HEC-RAS 2D scenarios of a handoff automatically (ras-commander), export the Max maps and score them.

For every run in the handoff's scenario_table.csv (flow scale x rheology):
    1. copy the HEC-RAS project into <runs_dir>/<run> (the original project is never modified),
    2. write that run's inflow hydrograph into the unsteady-flow file and make the simulation long enough,
    3. run the plan with ras-commander (HEC-RAS stays closed),
    4. check the run (completed?, error lines, volume-accounting error),
    5. export Depth (Max), Velocity (Max) and Arrival Time as GeoTIFF into <runs_dir>/<run>/maps,
    6. score them (fanrisk.ras.readback) against the site's validation case and its AHP map,
and write <runs_dir>/batch_summary.csv with one row per run.

Flow scales and 'bulked clear water' are automatic (bulked flow = clear-water flow / (1 - Cv)). 'non-Newtonian (Bingham)' needs a plan made
once in the HEC-RAS GUI: put its number in [hecras.rheology_plans] of site.toml; until then those runs are skipped with a message.
This is ras_batch.py as tested on HEC-RAS 7.0.2; the settings now come from site.toml (configure()).
"""
import csv
import datetime as _dt
import json
import math
import os
import re
import shutil
import sys
import time
from pathlib import Path

# =============================================================================
# SETTINGS (filled by configure() from the site file)
# =============================================================================
PROJECT_PRJ = r""                    # the HEC-RAS project (.prj) - never modified
RAS_VERSION = "7.0.2"                # the installed HEC-RAS version
HANDOFF_DIR = r""                    # the handoff folder written by 'fanrisk run'
RUNS_DIR = r""                       # "" = <project folder>_runs (next to the project)
BASE_PLAN = "01"                     # the clear-water plan
RHEOLOGY_PLANS = {"non-Newtonian (Bingham)": None}
INFLOW_BC_NAME = ""                  # name of the inflow boundary line; "" = the only flow-hydrograph boundary in the unsteady file
RUN_ONLY = []
TAIL_MIN = 30.0                      # keep simulating this long after the inflow hydrograph ends (the flood must leave the model)
HYDRO_INTERVAL_MIN = 1               # interval of the written hydrograph
NUM_CORES = 4                        # cores per run
MAX_RUNTIME_S = 3 * 3600             # stop a run that takes longer than this (None = no limit)
EXPORT_ARRIVAL = True                # also export Arrival Time (falls back to depth + velocity if the export refuses it)
ARRIVAL_DEPTH_M = 0.3                # depth that counts as 'arrived' (converted to feet for an English-unit project)
KEEP_RUN_RESULTS = True              # False = delete each run's copied project after its maps are exported (keeps maps + scores)
SMOKE_MINUTES = 10                   # length of the 'smoke' test run
SITE = None                          # the Site (scoring settings, AHP map, observations)

MON = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
M3S_TO_CFS = 35.31466672


def configure(site, runs=None):
    """Take every setting from a Site (config.load_site)."""
    global PROJECT_PRJ, RAS_VERSION, HANDOFF_DIR, RUNS_DIR, BASE_PLAN, RHEOLOGY_PLANS, INFLOW_BC_NAME, RUN_ONLY, TAIL_MIN, HYDRO_INTERVAL_MIN, NUM_CORES
    global MAX_RUNTIME_S, EXPORT_ARRIVAL, ARRIVAL_DEPTH_M, KEEP_RUN_RESULTS, SMOKE_MINUTES, SITE
    h = site.hecras
    SITE = site
    PROJECT_PRJ = site.file("project")
    if not PROJECT_PRJ:
        raise SystemExit("[hecras] project is missing in site.toml: the .prj of the HEC-RAS model (build it once, or make it with 'fanrisk hecras fit').")
    if not os.path.exists(PROJECT_PRJ):
        raise SystemExit("[hecras] project not found: %s (make it with 'fanrisk hecras fit')." % PROJECT_PRJ)
    HANDOFF_DIR = site.handoff_dir()
    if not HANDOFF_DIR:
        raise SystemExit("No handoff found in %s: run 'fanrisk run' with [stages] hecras_handoff = true first." % os.path.join(site.out_dir, "hecras"))
    RUNS_DIR = site.file("runs_dir") or os.path.join(site.out_dir, "hecras", "runs")
    RAS_VERSION, BASE_PLAN = str(h["version"]), str(h["base_plan"]).zfill(2)
    RHEOLOGY_PLANS = {"non-Newtonian (Bingham)": None}
    RHEOLOGY_PLANS.update({k: (str(v).zfill(2) if v else None) for k, v in (h.get("rheology_plans") or {}).items()})
    INFLOW_BC_NAME, RUN_ONLY = str(h["inflow_bc_name"]), list(runs if runs else h["runs"])
    TAIL_MIN, HYDRO_INTERVAL_MIN, NUM_CORES = float(h["tail_min"]), int(h["hydrograph_interval_min"]), int(h["cores"])
    MAX_RUNTIME_S = float(h["max_runtime_h"]) * 3600 if h["max_runtime_h"] else None
    EXPORT_ARRIVAL, ARRIVAL_DEPTH_M = bool(h["export_arrival"]), float(h["arrival_depth_m"])
    KEEP_RUN_RESULTS, SMOKE_MINUTES = bool(h["keep_run_results"]), int(h["smoke_minutes"])


# =============================================================================
# plain-text HEC-RAS files (no ras-commander needed; tested without HEC-RAS)
# =============================================================================
def read_text(p):
    with open(p, "r", encoding="latin-1", newline="") as fh:
        return fh.read()


def write_text(p, s):
    with open(p, "w", encoding="latin-1", newline="") as fh:
        fh.write(s)


def eol_of(s):
    return "\r\n" if "\r\n" in s else "\n"


def project_files(prj):
    """Project name, folder, units, plan list and plan -> (flow file, geometry file) from the .prj and .p## files."""
    prj = Path(prj)
    s = read_text(prj)
    units = "SI" if re.search(r"^\s*SI Units", s, re.M) else ("English" if re.search(r"^\s*English Units", s, re.M) else "unknown")
    plans = re.findall(r"^Plan File=p(\d+)", s, re.M)
    out = dict(name=prj.stem, folder=prj.parent, units=units, plans={})
    for p in plans:
        pf = prj.parent / ("%s.p%s" % (prj.stem, p))
        if pf.exists():
            out["plans"][p] = plan_info(pf)
    return out


def plan_info(pf):
    s = read_text(pf)
    g = lambda key: (re.search(r"^%s=(.*)$" % re.escape(key), s, re.M).group(1).strip() if re.search(r"^%s=" % re.escape(key), s, re.M) else "")
    return dict(path=Path(pf), title=g("Plan Title"), short_id=g("Short Identifier"), flow=g("Flow File"), geom=g("Geom File"),
                sim=g("Simulation Date"), comp_interval=g("Computation Interval"), output_interval=g("Output Interval"),
                mapping_interval=g("Mapping Interval"))


DEFAULT_DAY = "01JAN2000"                 # used when the plan has no simulation date (HEC-RAS then shows 01JAN1000 in its maps)


def parse_ras_dt(day, hhmm):
    """'01JAN2010' / '1jan2010' + '0130' / '01:30' / '1:30' / '' -> datetime. A blank date means DEFAULT_DAY."""
    day = (day or "").strip() or DEFAULT_DAY
    m = re.match(r"^(\d{1,2})([A-Za-z]{3})(\d{2,4})$", day)
    if not m:
        raise ValueError("cannot read the date '%s'" % day)
    yr = int(m.group(3))
    yr = yr + 2000 if len(m.group(3)) == 2 else yr
    d = _dt.datetime(yr, MON.index(m.group(2).upper()) + 1, int(m.group(1)))
    t = (hhmm or "").strip()
    if not t:
        return d
    if ":" in t:
        h, mi = (t.split(":") + ["0"])[:2]
    else:
        t = t.zfill(4)
        h, mi = t[:-2], t[-2:]
    return d + _dt.timedelta(hours=int(h or 0), minutes=int(mi or 0))


def fmt_ras_dt(t):
    return "%02d%s%04d" % (t.day, MON[t.month - 1], t.year), "%02d%02d" % (t.hour, t.minute)


def sim_window(sim):
    """'01JAN1000,0000,01JAN1000,0130' -> (start, end) datetimes."""
    a = [x.strip() for x in sim.split(",")]
    a += [""] * (4 - len(a))
    try:
        s, e = parse_ras_dt(a[0], a[1]), parse_ras_dt(a[2] or a[0], a[3])
    except Exception as ex:
        raise ValueError("cannot read 'Simulation Date=%s' (%s)" % (sim, ex))
    if e <= s:
        raise ValueError("'Simulation Date=%s' ends before it starts - set the start and end in the plan (Unsteady Flow Analysis window)" % sim)
    return s, e


def set_sim_end(plan_path, new_end):
    s = read_text(plan_path)
    m = re.search(r"^Simulation Date=(.*)$", s, re.M)
    a = [x.strip() for x in m.group(1).split(",")]
    a += [""] * (4 - len(a))
    if not a[0]:                                     # only a blank start date is written out; an existing start is kept exactly as HEC-RAS wrote it
        a[0], a[1] = fmt_ras_dt(parse_ras_dt(a[0], a[1]))
    d, h = fmt_ras_dt(new_end)
    a[2], a[3] = d, h
    write_text(plan_path, s[:m.start(1)] + ",".join(a) + s[m.end(1):])


def interval_minutes(txt):
    m = re.match(r"\s*(\d+)\s*(SEC|MIN|HOUR|DAY|WEEK)", txt.upper())
    if not m:
        raise ValueError("unknown interval '%s'" % txt)
    return int(m.group(1)) * {"SEC": 1 / 60.0, "MIN": 1.0, "HOUR": 60.0, "DAY": 1440.0, "WEEK": 10080.0}[m.group(2)]


def flow_bc_blocks(u_text):
    """Every flow-hydrograph boundary of an unsteady file: dict(name, location, start/end line, interval line, header line, value lines, values, use_dss, fixed_start)."""
    L = u_text.splitlines()
    starts = [i for i, l in enumerate(L) if l.startswith("Boundary Location=")]
    out = []
    for k, i in enumerate(starts):
        j = starts[k + 1] if k + 1 < len(starts) else len(L)
        hdr = next((x for x in range(i, j) if L[x].startswith("Flow Hydrograph=")), None)
        if hdr is None:
            continue
        n = int(re.search(r"=\s*(\d+)", L[hdr]).group(1))
        nlines = int(math.ceil(n / 10.0))
        vals = []
        for x in range(hdr + 1, hdr + 1 + nlines):
            row = L[x]
            vals += [float(row[c:c + 8]) for c in range(0, len(row), 8) if row[c:c + 8].strip()]
        loc = [f.strip() for f in L[i].split("=", 1)[1].split(",")]
        name = next((f for f in reversed(loc) if f), "")
        intv = next((x for x in range(i, hdr) if L[x].startswith("Interval=")), None)
        g = lambda key: next((L[x].split("=", 1)[1].strip() for x in range(i, j) if L[x].startswith(key + "=")), "")
        out.append(dict(name=name, location=L[i], area=loc[5] if len(loc) > 5 else "", start=i, end=j, interval_line=intv, header=hdr,
                        nlines=nlines, n=n, values=vals, interval=L[intv].split("=", 1)[1].strip() if intv is not None else "",
                        use_dss=g("Use DSS").lower() == "true", fixed_start=g("Use Fixed Start Time").lower() == "true"))
    return out


def pick_bc(blocks, name):
    if name:
        b = [x for x in blocks if x["name"].lower() == name.lower()]
        if len(b) != 1:
            raise RuntimeError("inflow boundary '%s' not found exactly once; flow-hydrograph boundaries: %s" % (name, [x["name"] for x in blocks]))
        return b[0]
    if len(blocks) != 1:
        raise RuntimeError("the unsteady file has %d flow-hydrograph boundaries (%s): set INFLOW_BC_NAME" % (len(blocks), [x["name"] for x in blocks]))
    return blocks[0]


def fmt8(v):
    for f in ("%8.2f", "%8.1f", "%8.0f"):
        s = f % v
        if len(s) <= 8:
            return s
    raise ValueError("value %g does not fit an 8-character HEC-RAS field" % v)


def replace_hydrograph(u_text, bc_name, values, interval_txt):
    """New unsteady-file text with the flow hydrograph of one boundary replaced (values at interval_txt)."""
    eol = eol_of(u_text)
    L = u_text.splitlines()
    b = pick_bc(flow_bc_blocks(u_text), bc_name)
    rows = ["".join(fmt8(v) for v in values[k:k + 10]) for k in range(0, len(values), 10)]
    new = L[:b["header"]] + ["Flow Hydrograph= %d " % len(values)] + rows + L[b["header"] + 1 + b["nlines"]:]
    if b["interval_line"] is not None:
        new[b["interval_line"]] = "Interval=" + interval_txt
    else:
        new.insert(b["header"], "Interval=" + interval_txt)
    shift = len(new) - len(L)                                # the hydrograph must follow the simulation clock, not a fixed date
    for x in range(b["start"], b["end"] + shift):
        if new[x].startswith("Use Fixed Start Time="):
            new[x] = "Use Fixed Start Time=False"
    return eol.join(new) + (eol if u_text.endswith(("\n", "\r")) else "")


def load_hydro_csv(p):
    t, q = [], []
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            t.append(float(r["time_min"]))
            q.append(float(r["Q_m3s"]))
    return t, q


def resample(t, q, dt, t_end):
    """Values every dt minutes from 0 to t_end (linear; 0 after the hydrograph ends)."""
    n = int(math.floor(t_end / dt + 1e-9)) + 1
    out = []
    for k in range(n):
        x = k * dt
        if x <= t[0]:
            out.append(q[0] if x == t[0] else 0.0)
        elif x >= t[-1]:
            out.append(q[-1] if x == t[-1] else 0.0)
        else:
            i = next(i for i in range(1, len(t)) if t[i] >= x)
            out.append(q[i - 1] + (q[i] - q[i - 1]) * (x - t[i - 1]) / (t[i] - t[i - 1]))
    return out


def read_scenarios(handoff):
    p = Path(handoff) / "scenario_table.csv"
    if not p.exists():
        raise SystemExit("No scenario_table.csv in %s - the handoff needs a hydrograph (set RAS_OBSERVED_HYDROGRAPH or event rainfall)." % handoff)
    with open(p, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("flow_scale", "peak_m3s", "volume_m3", "Cv", "bulking_factor", "n_channel"):
            r[k] = float(r[k]) if r.get(k) not in (None, "", "nan") else float("nan")
    return rows


def plan_for(run, base_plan, rheology_plans):
    """(plan number, hydrograph multiplier, note) for one scenario row; plan None = skip."""
    rh = run["rheology"]
    if rh == "clear water":
        return base_plan, 1.0, ""
    if rh == "bulked clear water":
        cv = run["Cv"]
        return base_plan, 1.0 / (1.0 - cv), "bulked: clear-water flow / (1 - Cv), Cv = %.2f" % cv
    p = rheology_plans.get(rh)
    if not p:
        return None, 1.0, "skipped: no HEC-RAS plan for '%s' yet (make it in the GUI and add it to RHEOLOGY_PLANS)" % rh
    return p, 1.0, "rheology from plan %s (set Cv %.2f there)" % (p, run["Cv"])


def copy_project(src_folder, dst_folder):
    """Copy the project without old plan results (they are recomputed)."""
    if os.path.exists(dst_folder):
        shutil.rmtree(dst_folder)
    ign = shutil.ignore_patterns("*.p[0-9][0-9].hdf", "*.p[0-9][0-9].tmp.hdf", "*.bak", "*.single-plan.orig", "*_runs", "~*")
    shutil.copytree(src_folder, dst_folder, ignore=ign)


def prepare_run(prj, plan_no, hydro_csv, mult, units, bc_name, tail_min, dt_min, log):
    """Write the hydrograph into the plan's unsteady file (in a COPY) and make the simulation long enough. Returns facts for the summary."""
    pinfo = plan_info(Path(prj).parent / ("%s.p%s" % (Path(prj).stem, plan_no)))
    u_path = Path(prj).parent / ("%s.%s" % (Path(prj).stem, pinfo["flow"]))
    t, q = load_hydro_csv(hydro_csv)
    q = [v * mult for v in q]
    conv = M3S_TO_CFS if units == "English" else 1.0
    start, end = sim_window(pinfo["sim"])
    need = max(t) + tail_min
    have = (end - start).total_seconds() / 60.0
    note = []
    if have < need:
        set_sim_end(pinfo["path"], start + _dt.timedelta(minutes=need))
        note.append("simulation lengthened from %.0f to %.0f min" % (have, need))
        have = need
    vals = [v * conv for v in resample(t, q, dt_min, have)]
    u_old = read_text(u_path)
    bc = pick_bc(flow_bc_blocks(u_old), bc_name)
    if bc["use_dss"]:
        raise RuntimeError("the inflow boundary '%s' reads from DSS (Use DSS=True): switch it to a table in the GUI" % bc["name"])
    if bc["fixed_start"]:
        note.append("inflow hydrograph switched from 'Fixed Start Time' to 'Use Simulation Time' in this copy")
    u_new = replace_hydrograph(u_old, bc["name"], vals, "%dMIN" % dt_min)
    write_text(u_path, u_new)
    back = pick_bc(flow_bc_blocks(read_text(u_path)), bc["name"])      # read back what HEC-RAS will read
    if back["fixed_start"]:
        raise RuntimeError("could not switch off 'Fixed Start Time' for the inflow boundary in %s - nothing was run" % u_path.name)
    if back["n"] != len(vals) or max(abs(a - b) for a, b in zip(back["values"], vals)) > 0.06 or back["interval"].upper() != "%dMIN" % dt_min:
        raise RuntimeError("the hydrograph written to %s does not read back correctly - nothing was run" % u_path.name)
    peak_in = max(q) * conv
    if max(vals) < 0.97 * peak_in:
        note.append("WARNING: the %d-min interval cuts the peak from %.0f to %.0f" % (dt_min, peak_in, max(vals)))
    return dict(unsteady_file=u_path.name, bc=bc["name"], values=len(vals), peak_written=max(vals), units=units, sim_minutes=have,
                volume_written_m3=sum((a + b) / 2.0 for a, b in zip(vals[:-1], vals[1:])) * dt_min * 60.0 / conv, notes="; ".join(note))


# =============================================================================
# ras-commander part (Windows + HEC-RAS)
# =============================================================================
def rc():
    try:
        import ras_commander as R
    except ImportError as e:
        raise SystemExit("ras-commander is not importable in this Python (%s). Use the Python where you ran test_ras.py." % e)
    return R


def call_flex(fn, *a, **kw):
    """Call fn, dropping keyword arguments this ras-commander version does not know (and saying so)."""
    while True:
        try:
            return fn(*a, **kw)
        except TypeError as e:
            m = re.search(r"unexpected keyword argument '(\w+)'", str(e))
            if not m or m.group(1) not in kw:
                raise
            print("   (this ras-commander version has no '%s' option - continuing without it)" % m.group(1))
            kw.pop(m.group(1))


def open_project(prj):
    R = rc()
    obj = R.RasPrj()
    call_flex(R.init_ras_project, str(prj), RAS_VERSION, ras_object=obj)
    return obj


def run_messages(R, hdf):
    msgs = ""
    try:
        msgs = R.HdfResultsPlan.get_compute_messages(str(hdf)) or ""
    except Exception as e:
        msgs = ""
        print("   (could not read compute messages: %s)" % e)
    bad = [l.strip() for l in msgs.splitlines()
           if re.search(r"\b(ERROR|FAILED|UNSTABLE|ABORT)", l.upper()) and not re.search(r"VOLUME ACCOUNTING ERROR|VOLUME ERROR|ERROR AS PERCENT", l.upper())]
    m = re.search(r"Volume Accounting Error as percentage:\s*([-\d.Ee+]+)", msgs, re.I) or re.search(r"percentage[^:\n]*:\s*([-\d.Ee+]+)", msgs, re.I)
    vol = float(m.group(1)) if m else float("nan")
    done = bool(re.search(r"Complete Process|Finished Unsteady Flow Simulation", msgs, re.I))
    return dict(error_lines=bad[:10], volume_error_pct=vol, finished_message=done, messages=msgs)


def compute_one(prj_copy, plan_no, log):
    R = rc()
    obj = open_project(prj_copy)
    t0 = time.time()
    res = call_flex(R.RasCmdr.compute_plan, plan_no, ras_object=obj, num_cores=NUM_CORES, clear_geompre=False, force_rerun=True,
                    verify=True, max_runtime=MAX_RUNTIME_S)
    dt = time.time() - t0
    ok = bool(res)
    verified = getattr(res, "completion_verified", None)
    hdf = Path(prj_copy).parent / ("%s.p%s.hdf" % (Path(prj_copy).stem, plan_no))
    rm = run_messages(R, hdf) if hdf.exists() else dict(error_lines=["no result HDF written"], volume_error_pct=float("nan"), finished_message=False, messages="")
    details = getattr(res, "execution_details", {}) or {}
    if details.get("runtime_timed_out"):
        rm["error_lines"].insert(0, "stopped after MAX_RUNTIME_S")
    return dict(success=ok and hdf.exists(), completion_verified=verified, minutes=dt / 60.0, hdf=str(hdf), **{k: v for k, v in rm.items() if k != "messages"}), rm["messages"]


def export_maps(prj_copy, plan_no, out_dir, units, log):
    hdf_only = os.environ.get("FANRISK_MAPS", "").lower() == "hdf"
    if not hdf_only:
        R = rc()
        obj = open_project(prj_copy)
        try:
            call_flex(R.RasMap.ensure_results_plan_layer, plan_no, ras_object=obj)
        except Exception as e:
            print("   (results layer not registered automatically: %s)" % e)
    os.makedirs(out_dir, exist_ok=True)
    arr_depth = ARRIVAL_DEPTH_M * (3.28084 if units == "English" else 1.0)
    tries = [dict(depth=True, velocity=True, arrival_time=True, arrival_depth=arr_depth)] if EXPORT_ARRIVAL else []
    tries.append(dict(depth=True, velocity=True))
    last = None
    if hdf_only:                                          # skip RAS Mapper (e.g. its helper is blocked by Windows security)
        tries, last = [], "skipped (FANRISK_MAPS=hdf)"
    for flags in tries:
        try:
            summ = call_flex(R.RasMap.store_all_maps, plan_no, ras_object=obj, mode="selected", profile="Max", output_path=str(out_dir),
                             timeout=1800, **flags)
            pl = (summ.get("plans") or {}).get(plan_no) or next(iter((summ.get("plans") or {}).values()), {})
            if summ.get("success") and pl.get("success"):
                return dict(maps_ok=True, map_files=pl.get("files", []), arrival="arrival_time" in flags)
            last = pl.get("error") or "store_all_maps reported failure"
        except Exception as e:
            last = "%s: %s" % (type(e).__name__, e)
        if "arrival_time" in flags:
            print("   arrival-time export refused (%s); retrying with depth + velocity only" % last)
        if "WinError 225" in str(last):                     # Windows security blocked RAS Mapper's helper: retrying cannot help
            break
    print("   RAS Mapper export failed (%s)\n   -> making the Max maps directly from the results HDF instead" % last)
    try:
        from . import hdfmaps
        files = hdfmaps.maps_from_hdf(*_hdf_inputs(prj_copy, plan_no), str(out_dir), print)
        return dict(maps_ok=True, map_files=files, arrival=False, via="results HDF")
    except Exception as e:
        return dict(maps_ok=False, map_error="%s | HDF fallback: %s: %s" % (last, type(e).__name__, e))


def _hdf_inputs(prj_copy, plan_no):
    prj_copy = Path(prj_copy)
    hdf = prj_copy.with_name("%s.p%s.hdf" % (prj_copy.stem, str(plan_no).zfill(2)))
    if not hdf.exists():
        raise RuntimeError("no results file %s" % hdf)
    tdir = prj_copy.parent / "Terrain"
    tifs = sorted(tdir.glob("terrain_handoff.tif")) or sorted(f for f in tdir.glob("*.tif") if not f.name.lower().startswith("terrain."))
    if not tifs:
        raise RuntimeError("no terrain GeoTIFF in %s" % tdir)
    return str(hdf), str(tifs[0])


def score(run_dir, log):
    """Score the exported maps of one run in this process (no GDAL needed). The boundary lines come from the run's own geometry file."""
    from . import geometry, readback
    lines = []

    def keep(msg):
        lines.append(str(msg))
    try:
        bc = None
        proj = Path(run_dir) / "project"
        if proj.exists():
            try:
                bc = list(geometry.read_bc_lines(geometry.geometry_file(str(proj))).values()) or None
            except Exception:
                bc = None
        ahp = os.path.join(SITE.out_dir, "risk", "ahp_C%d_risk_abs.tif" % json.load(open(os.path.join(HANDOFF_DIR, "handoff_summary.json")))["complex"])
        readback.score_maps(HANDOFF_DIR, str(Path(run_dir) / "maps"), str(Path(run_dir) / "readback"), readback.settings_from_site(SITE), ahp_risk_file=ahp,
                            footprint=SITE.file("flood_footprint"), points=SITE.file("flood_points"), allow_snapshot=bool(SITE.score["allow_snapshot"]),
                            bc_lines=bc, log=keep)
        ok, err = True, ""
    except (Exception, SystemExit) as e:
        ok, err = False, "%s: %s" % (type(e).__name__, e)
        lines.append(err)
    with open(Path(run_dir) / "readback.log", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return dict(score_ok=True) if ok else dict(score_ok=False, score_error=err)


def runs_dir():
    return Path(RUNS_DIR) if RUNS_DIR else Path(str(Path(PROJECT_PRJ).parent) + "_runs")


def write_summary(rows_by_run):
    rows = []
    for rid, r in sorted(rows_by_run.items()):
        row = dict(r)
        rbj = runs_dir() / rid / "readback" / "readback_summary.json"
        if rbj.exists():
            j = json.load(open(rbj))
            row.update(wet_area_km2=j["model"]["wet_area_km2"], max_depth_m=j["model"]["max_depth_m"], max_speed_ms=j["model"]["max_speed_ms"],
                       wet_cells_at_edge=j["model"].get("wet_cells_at_model_edge"), maximum_maps=j.get("files_used", {}).get("maximum_maps"))
            for p in j.get("probes", []):
                row.update(probe_depth_m=p["depth_model_m"], probe_depth_obs_m=p["depth_obs_m"], probe_speed_ms=p["speed_model_ms"], probe_speed_obs_ms=p["speed_obs_ms"])
                break
            if "footprint" in j:
                row.update(footprint_CSI=j["footprint"]["CSI"], footprint_POD=j["footprint"]["POD"], footprint_FAR=j["footprint"]["FAR"])
            if "points" in j:
                row.update(points_hit_rate=j["points"]["hit_rate"])
            if "ahp" in j:
                row.update(ahp_auc_wet=j["ahp"]["roc_auc_risk_vs_wet"], ahp_rho_hazard=j["ahp"]["spearman_risk_vs_hazard_rating"])
        rows.append(row)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    p = runs_dir() / "batch_summary.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and v != v) else v) for k, v in r.items()})
    return p


def state_path():
    return runs_dir() / "batch_state.json"


def load_state():
    p = state_path()
    return json.load(open(p)) if p.exists() else {}


def save_state(st):
    os.makedirs(runs_dir(), exist_ok=True)
    with open(state_path(), "w") as fh:
        json.dump(st, fh, indent=2, default=str)


# =============================================================================
# commands
# =============================================================================
def cmd_check():
    ok = True

    def say(flag, text):
        nonlocal ok
        ok &= flag is not False
        print("  [%s] %s" % ({True: " OK ", False: "FAIL", None: "WARN"}[flag], text))
    print("1. Python and ras-commander")
    say(True, "Python %s at %s" % (sys.version.split()[0], sys.executable))
    try:
        R = rc()
        say(True, "ras-commander %s" % getattr(R, "__version__", "?"))
        for name in ("init_ras_project", "RasPrj", "RasCmdr", "RasMap", "HdfResultsPlan"):
            say(hasattr(R, name), "ras_commander.%s %s" % (name, "found" if hasattr(R, name) else "MISSING"))
    except SystemExit as e:
        say(False, str(e))
        return
    print("2. Your project (read-only)")
    prj = Path(PROJECT_PRJ)
    if not prj.exists():
        say(False, "PROJECT_PRJ does not exist: %s" % prj)
        return
    pf = project_files(prj)
    say(True if pf["units"] == "SI" else None, "units: %s%s" % (pf["units"], "" if pf["units"] == "SI" else " - hydrographs are converted to cfs, but the scoring reads depths as METRES: use an SI project"))
    say(bool(pf["plans"]), "plans: %s" % ", ".join("%s (%s, flow %s, geometry %s)" % (k, v["short_id"] or v["title"], v["flow"], v["geom"]) for k, v in pf["plans"].items()))
    for k, v in pf["plans"].items():                       # every plan's window, to help choose BASE_PLAN
        hdf = prj.parent / ("%s.p%s.hdf" % (prj.stem, k))
        try:
            s_, e_ = sim_window(v["sim"])
            win = "%s -> %s (%.0f min)" % (s_.strftime("%d%b%Y %H:%M"), e_.strftime("%d%b%Y %H:%M"), (e_ - s_).total_seconds() / 60)
        except Exception:
            win = "NO VALID SIMULATION DATES ('%s')" % v["sim"]
        print("        plan %s  '%s' / %s: %s | results file: %s" % (k, v["title"], v["short_id"], win, "yes (%.0f MB)" % (hdf.stat().st_size / 1e6) if hdf.exists() else "no"))
    plans_needed = [BASE_PLAN] + [p for p in RHEOLOGY_PLANS.values() if p]
    for p in plans_needed:
        if p not in pf["plans"]:
            say(False, "plan %s is not in the project" % p)
            continue
        pi = pf["plans"][p]
        try:
            s, e = sim_window(pi["sim"])
            say(True, "plan %s: simulation %s -> %s (%.0f min), computation interval %s, mapping interval %s" % (p, s, e, (e - s).total_seconds() / 60, pi["comp_interval"], pi["mapping_interval"] or "?"))
        except Exception as ex:
            say(False, "plan %s: %s - send this line to fix it" % (p, ex))
        u = prj.parent / ("%s.%s" % (prj.stem, pi["flow"]))
        if not u.exists():
            say(False, "plan %s: unsteady file %s is missing (is it an unsteady plan?)" % (p, u.name))
            continue
        blocks = flow_bc_blocks(read_text(u))
        try:
            b = pick_bc(blocks, INFLOW_BC_NAME)
            say(not b["use_dss"], "plan %s: inflow boundary '%s' on 2D area '%s', %d values at %s%s" % (p, b["name"], b["area"], b["n"], b["interval"], ", reads DSS - switch to a table" if b["use_dss"] else ""))
            if b["fixed_start"]:
                say(None, "plan %s: inflow hydrograph uses 'Fixed Start Time' - the runs switch it to 'Use Simulation Time' in their copies; do the same in the GUI for your own runs" % p)
        except RuntimeError as ex:
            say(False, "plan %s: %s" % (p, ex))
        if "Normal Depth" not in read_text(u) and "Friction Slope=" not in read_text(u):
            say(None, "plan %s: no normal-depth (friction slope) outflow boundary found in %s" % (p, u.name))
        hdf = prj.parent / ("%s.p%s.hdf" % (prj.stem, p))
        if hdf.exists():
            rm = run_messages(R, hdf)
            say(None if rm["error_lines"] else True, "plan %s has old results: volume error %s %%, error lines: %s" % (p, "%.4f" % rm["volume_error_pct"] if rm["volume_error_pct"] == rm["volume_error_pct"] else "n/a", rm["error_lines"][:3] or "none"))
    rasmap = prj.with_suffix(".rasmap")
    if rasmap.exists():
        rs = read_text(rasmap)
        m = re.search(r'<RASProjectionFilename Filename="([^"]+)"', rs)
        pth = (prj.parent / m.group(1).replace("\\", "/").lstrip("./")) if m else None
        say(True if (pth and pth.exists()) else None, "projection: %s" % (str(pth) if pth and pth.exists() else
            "NOT SET - RAS Mapper > Project > Set Projection (use the terrain's EPSG:32643); without it the exported maps may lack georeferencing"))
        say(True if "<Terrains" in rs and "Type=\"TerrainLayer\"" in rs else None, "terrain layer registered in RAS Mapper")
    print("3. ras-commander can open the project")
    try:
        obj = open_project(prj)
        say(True, "init_ras_project OK (%d plans in plan_df)" % len(obj.plan_df))
    except Exception as e:
        say(False, "init_ras_project failed: %s" % e)
    print("4. Handoff package")
    hd = Path(HANDOFF_DIR) if HANDOFF_DIR else None
    if not hd or not hd.exists():
        say(False, "HANDOFF_DIR not set / missing")
    else:
        for f in ("scenario_table.csv", "handoff_summary.json", "flow_area_polygon.geojson", "hydrograph_Qx1.00.csv"):
            say((hd / f).exists(), f)
        if (hd / "scenario_table.csv").exists():
            sc = read_scenarios(hd)
            for r in sc:
                p, mult, note = plan_for(r, BASE_PLAN, RHEOLOGY_PLANS)
                print("        %s  Q x %.2f  %-24s -> %s %s" % (r["run"], r["flow_scale"], r["rheology"], ("plan " + p) if p else "SKIPPED", note))
    print("5. Hydrograph write test (temporary copy, deleted afterwards)")
    if hd and (hd / "hydrograph_Qx1.00.csv").exists() and BASE_PLAN in pf["plans"]:
        tmp = Path(str(runs_dir())) / "_check_copy"
        try:
            copy_project(prj.parent, tmp)
            info = prepare_run(tmp / prj.name, BASE_PLAN, hd / "hydrograph_Qx1.00.csv", 1.0, pf["units"], INFLOW_BC_NAME, TAIL_MIN, HYDRO_INTERVAL_MIN, print)
            say(True, "wrote %d values (peak %.0f) into %s and read them back; simulation %.0f min %s" % (info["values"], info["peak_written"], info["unsteady_file"], info["sim_minutes"], info["notes"]))
        except Exception as e:
            say(False, "write test failed: %s" % e)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    print("6. Scoring")
    from .. import gio
    say(True, "maps are read with the '%s' backend; observations: %s" % (gio.backend(), (SITE.case or {}).get("name") or "none (no validation case)"))
    print("\nNOTE: keep the project CLOSED in HEC-RAS while runs execute (HEC-RAS 7.0.2 single-plan guard).")
    print("SETUP %s" % ("OK" if ok else "HAS PROBLEMS - fix the FAIL lines first"))


def cmd_smoke():
    pf = project_files(PROJECT_PRJ)
    d = runs_dir() / "_smoke"
    print("copying the project to %s ..." % d)
    copy_project(Path(PROJECT_PRJ).parent, d)
    prj = d / Path(PROJECT_PRJ).name
    pi = plan_info(d / ("%s.p%s" % (prj.stem, BASE_PLAN)))
    s, e = sim_window(pi["sim"])
    set_sim_end(pi["path"], s + _dt.timedelta(minutes=SMOKE_MINUTES))
    print("plan file now says: Simulation Date=%s" % plan_info(pi["path"])["sim"])
    print("running plan %s for %d simulated minutes ..." % (BASE_PLAN, SMOKE_MINUTES))
    res, msgs = compute_one(prj, BASE_PLAN, print)
    with open(d / "compute_messages.txt", "w") as fh:
        fh.write(msgs)
    print("success %s | completion verified %s | %.1f min wall time | volume error %s %% | error lines %s" % (
        res["success"], res["completion_verified"], res["minutes"], res["volume_error_pct"], res["error_lines"] or "none"))
    print("SMOKE TEST %s" % ("PASSED - HEC-RAS runs from Python" if res["success"] else "FAILED - see the messages above / in HEC-RAS"))
    if not res["success"]:
        print(SMOKE_HELP)
        raise SystemExit(1)


def _has_depth_map(d):
    d = Path(d)
    return d.is_dir() and any("depth" in f.name.lower() and f.suffix.lower() in (".tif", ".tiff", ".vrt") for f in d.iterdir())


SMOKE_HELP = """Usual causes, in order:
  1. HEC-RAS terms not accepted: open HEC-RAS (the installed version, e.g. 7.0.2) from the Start menu once, click 'I Agree', close it.
  2. psutil missing: pip install psutil
  3. the project is open in the HEC-RAS window: close HEC-RAS before running."""


def cmd_run(do_compute=True, do_maps=True, do_score=True):
    if not HANDOFF_DIR:
        raise SystemExit("Set HANDOFF_DIR first.")
    pf = project_files(PROJECT_PRJ)
    st = load_state()
    rows = read_scenarios(HANDOFF_DIR)
    if RUN_ONLY:
        rows = [r for r in rows if r["run"] in RUN_ONLY]
    for r in rows:
        rid = r["run"]
        rd = runs_dir() / rid
        S = st.setdefault(rid, {})
        S.update(run=rid, flow_scale=r["flow_scale"], rheology=r["rheology"], Cv=r["Cv"])
        plan, mult, note = plan_for(r, BASE_PLAN, RHEOLOGY_PLANS)
        print("\n=== %s: Q x %.2f, %s %s" % (rid, r["flow_scale"], r["rheology"], ("-> plan %s" % plan) if plan else ""))
        if note:
            print("   " + note)
        if not plan:
            S.update(status="skipped", note=note)
            save_state(st)
            continue
        prj = rd / "project" / Path(PROJECT_PRJ).name
        try:
            if do_compute:
                copy_project(Path(PROJECT_PRJ).parent, rd / "project")
                info = prepare_run(prj, plan, Path(HANDOFF_DIR) / r["hydrograph_file"], mult, pf["units"], INFLOW_BC_NAME, TAIL_MIN, HYDRO_INTERVAL_MIN, print)
                print("   hydrograph: %d values, peak %.0f, %.2f million m3; simulation %.0f min %s" % (info["values"], info["peak_written"], info["volume_written_m3"] / 1e6, info["sim_minutes"], info["notes"]))
                res, msgs = compute_one(prj, plan, print)
                with open(rd / "compute_messages.txt", "w") as fh:
                    fh.write(msgs)
                S.update(plan=plan, hydro_note=note, peak_written=info["peak_written"], sim_minutes=info["sim_minutes"], compute_minutes=round(res["minutes"], 2),
                         success=res["success"], completion_verified=res["completion_verified"], volume_error_pct=res["volume_error_pct"],
                         error_lines=" | ".join(res["error_lines"]), status="computed" if res["success"] else "FAILED")
                save_state(st)
                print("   computed: %s, verified %s, %.1f min, volume error %s %%%s" % (res["success"], res["completion_verified"], res["minutes"],
                      res["volume_error_pct"], (" | errors: " + "; ".join(res["error_lines"][:3])) if res["error_lines"] else ""))
                if not res["success"]:
                    continue
            if do_maps:
                if not prj.exists():
                    print("   no computed project for %s - run it first" % rid)
                    continue
                m = export_maps(prj, plan, rd / "maps", pf["units"], print)
                S.update(maps_ok=m["maps_ok"], map_error=m.get("map_error", ""), arrival_exported=m.get("arrival", False))
                save_state(st)
                print("   maps: %s" % ("%d files -> %s" % (len(m["map_files"]), rd / "maps") if m["maps_ok"] else "FAILED - " + str(m["map_error"])))
                if m["maps_ok"] and not KEEP_RUN_RESULTS:
                    shutil.rmtree(rd / "project", ignore_errors=True)
            if do_score:
                if not (rd / "maps").exists():
                    print("   no maps for %s - export them first" % rid)
                    continue
                sc = score(rd, print)
                S.update(**sc)
                save_state(st)
                print("   score: %s" % ("OK -> %s" % (rd / "readback") if sc["score_ok"] else "FAILED - " + sc["score_error"]))
        except Exception as e:
            S.update(status="FAILED", error="%s: %s" % (type(e).__name__, e))
            save_state(st)
            print("   FAILED: %s: %s" % (type(e).__name__, e))
    p = write_summary(st)
    print("\nsummary: %s" % p)
    if do_compute or do_maps:                             # a batch that produced no maps must not look like a success to scripts
        bad = [r["run"] for r in rows if plan_for(r, BASE_PLAN, RHEOLOGY_PLANS)[0] and not _has_depth_map(runs_dir() / r["run"] / "maps")]
        if bad:
            print("NO MAPS for %s: the run failed (see compute_messages.txt in its folder)." % ", ".join(bad))
            print(SMOKE_HELP)
            raise SystemExit(1)


def main(argv):
    cmd = (argv[1] if len(argv) > 1 else "check").lower()
    if cmd == "check":
        cmd_check()
    elif cmd == "smoke":
        cmd_smoke()
    elif cmd == "run":
        cmd_run(True, True, True)
    elif cmd == "maps":
        cmd_run(False, True, False)
    elif cmd == "score":
        cmd_run(False, False, True)
    elif cmd == "summary":
        print(write_summary(load_state()))
    else:
        raise SystemExit("usage: fanrisk hecras check | smoke | run | maps | score | summary")
