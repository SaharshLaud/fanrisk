# -*- coding: utf-8 -*-
"""Site configuration: one site.toml -> everything a run needs.

Four kinds of numbers are kept apart on purpose:
    method defaults   fanrisk/defaults.py          how fans are found and hazard is scored (not tied to a place)
    region profile    fanrisk/profiles/*.toml      starting assumptions for a type of terrain (runoff, roughness)
    site              your site.toml               the DEM, the storm, optional layers, what to simulate
    validation case   fanrisk/cases/*.toml         observations of one documented flood, used only to TEST results
Precedence: defaults < region profile < site.toml sections < [advanced].
"""
import copy
import difflib
import os

from .defaults import ENGINE_DEFAULTS

try:
    import tomllib as _toml
except ImportError:                                      # Python < 3.11
    try:
        import tomli as _toml
    except ImportError:
        _toml = None

HERE = os.path.dirname(os.path.abspath(__file__))


class ConfigError(Exception):
    pass


# site.toml key -> engine setting
ENGINE_KEYS = {
    "storm": dict(rainfall_mm="RAINFALL_DESIGN_MM", duration_min="STORM_DURATION_MIN", areal_reduction_exponent="STORM_ARF_EXPONENT",
                  core_km2="STORM_CORE_KM2", active_fraction="RAIN_ACTIVE_FRAC", rainfall_in_hazard="RAIN_IN_HAZARD"),
    "runoff": dict(method="RUNOFF_METHOD", coefficient="RUNOFF_COEFF", curve_number="CURVE_NUMBER", antecedent_moisture="CN_AMC",
                   initial_abstraction_ratio="CN_IA_RATIO"),
    "targets": dict(lonlat="AHP_TARGET_LONLAT", min_fan_km2="AHP_MIN_FAN_KM2", max_targets="AHP_MAX_TARGETS", apex_mode="APEX_MODE",
                    manual_apexes="MANUAL_APEXES", extra_feeders="EXTRA_FEEDERS"),
    "ahp": dict(class_breaks="RISK_CLASS_BREAKS", scenario_maps="SCENARIO_MAPS"),
    "hecras": dict(feeder="RAS_FEEDER", reach_km="RAS_REACH_KM", half_width_m="RAS_HALF_WIDTH_M", start_m="RAS_START_M", q_scales="RAS_Q_SCALES",
                   cell_size_m="RAS_CELL_SIZE_M", burn_depth_m="RAS_BURN_DEPTH_M"),
}
# site.toml keys that are file paths (relative to the site file)
FILE_KEYS = {
    "site": ("dem",),
    "storm": ("rainfall_raster",),
    "runoff": ("curve_number_raster", "landcover_raster", "soil_group_raster"),
    "layers": ("historical_rainfall", "bare_fraction", "built_density", "glacial_lakes", "reference_fans"),
    "validation": ("flood_footprint", "flood_points", "survey_window"),
    "hecras": ("centerline", "fine_dem", "project", "runs_dir"),
}
OTHER_KEYS = {
    "site": ("name", "output_dir", "region_profile", "validation_case"),
    "stages": ("risk", "hecras_handoff"),
    "storm": (), "layers": (), "validation": (),
    "runoff": ("default_soil_group", "cn_table"),
    "targets": ("mode", "ids"),
    "ahp": ("catchment_matrix", "fan_matrix"),
    "hecras": ("complex", "inflow", "version", "base_plan", "mannings_n", "rheology_plans", "inflow_bc_name", "runs", "tail_min", "hydrograph_interval_min", "cores",
               "max_runtime_h", "export_arrival", "arrival_depth_m", "smoke_minutes", "keep_run_results"),
    "score": ("wet_depth_m", "point_buffer_m", "probe_radius_m", "debris_factor", "hazard_breaks", "tolerance", "allow_snapshot"),
}
HECRAS_DEFAULTS = dict(complex=0, inflow="pipeline", mannings_n=None, version="7.0.2", base_plan="01", rheology_plans={}, inflow_bc_name="", runs=[], tail_min=30.0,
                       hydrograph_interval_min=1, cores=4, max_runtime_h=3.0, export_arrival=True, arrival_depth_m=0.3, smoke_minutes=10,
                       keep_run_results=True)
SCORE_DEFAULTS = dict(wet_depth_m=0.3, point_buffer_m=30.0, probe_radius_m=60.0, debris_factor=0.5, hazard_breaks=[0.75, 1.25, 2.5], tolerance=0.35,
                      allow_snapshot=False)


def load_toml(path):
    if _toml is None:
        raise ConfigError("Reading .toml files needs Python 3.11+ (or: pip install tomli).")
    if not os.path.exists(path):
        raise ConfigError("file not found: %s" % path)
    with open(path, "rb") as fh:
        try:
            return _toml.load(fh)
        except Exception as e:
            hint = ""
            with open(path, encoding="utf-8", errors="ignore") as th:
                if any('"' in l and "\\" in l.split("#")[0] for l in th):
                    hint = " Hint: write Windows paths in single quotes ('C:\\data\\dem.tif') or with forward slashes (\"C:/data/dem.tif\")."
            raise ConfigError("%s is not valid TOML: %s.%s" % (path, e, hint))


def available(kind):
    d = os.path.join(HERE, kind)
    return sorted(f[:-5] for f in os.listdir(d) if f.endswith(".toml"))


def _find(kind, name, root):
    """A profile / case by name (shipped with fanrisk) or by path (absolute, or relative to the site file)."""
    if not name:
        return None
    for p in (os.path.join(HERE, kind, name + ".toml"), name if os.path.isabs(name) else os.path.join(root, name)):
        if os.path.isfile(p):
            return p
    raise ConfigError("%s '%s' not found. Shipped with fanrisk: %s. Or give the path of your own .toml file." % (kind[:-1], name, ", ".join(available(kind))))


def _suggest(word, options):
    m = difflib.get_close_matches(word, list(options), n=1)
    return " (did you mean '%s'?)" % m[0] if m else ""


class Site(object):
    """Everything one run needs. P = engine settings; files = absolute paths ('' = not given)."""

    def __init__(self):
        self.path = self.root = self.name = self.out_dir = ""
        self.P, self.files, self.stages = {}, {}, {}
        self.profile, self.case = {}, None
        self.hecras, self.score = {}, {}
        self.target_mode, self.target_ids = "all", []
        self.cn_table, self.default_soil_group = None, None
        self.catchment_matrix = self.fan_matrix = None
        self.sources = {}                                  # engine setting -> where its value came from
        self.notes = []

    def file(self, key):
        return self.files.get(key, "")

    def handoff_dir(self, cid=None, feeder=None):
        base = os.path.join(self.out_dir, "hecras")
        if cid is None:
            if os.path.isdir(base):
                d = sorted(x for x in os.listdir(base) if x.startswith("handoff_") and os.path.isfile(os.path.join(base, x, "handoff_summary.json")))
                if d:
                    return os.path.join(base, d[0])
            return ""
        return os.path.join(base, "handoff_C%d_%s" % (cid, feeder))


def load_site(path):
    path = os.path.abspath(path)
    if os.path.isdir(path):
        path = os.path.join(path, "site.toml")
    doc = load_toml(path)
    s = Site()
    s.path, s.root = path, os.path.dirname(path)
    known = set(ENGINE_KEYS) | set(FILE_KEYS) | set(OTHER_KEYS) | {"advanced"}
    for sec in doc:
        if sec not in known:
            raise ConfigError("%s: unknown section [%s]%s. Sections: %s" % (os.path.basename(path), sec, _suggest(sec, known), ", ".join(sorted(known))))
        if not isinstance(doc[sec], dict):
            raise ConfigError("%s: '%s' must be a section, written as [%s]" % (os.path.basename(path), sec, sec))
    site = doc.get("site", {})
    s.name = str(site.get("name") or os.path.basename(s.root))
    out = site.get("output_dir") or "output"
    s.out_dir = os.path.normpath(out if os.path.isabs(out) else os.path.join(s.root, out))

    P = copy.deepcopy(ENGINE_DEFAULTS)
    s.sources = {}

    def put(key, val, src):
        if key not in ENGINE_DEFAULTS:
            raise ConfigError("%s: '%s' is not an engine setting%s. See fanrisk/defaults.py for the list." % (src, key, _suggest(key, ENGINE_DEFAULTS)))
        P[key] = tuple(val) if isinstance(ENGINE_DEFAULTS[key], tuple) and isinstance(val, list) else val
        s.sources[key] = src

    prof_path = _find("profiles", site.get("region_profile") or "generic", s.root)
    prof = load_toml(prof_path)
    s.profile = dict(prof.get("profile", {}), path=prof_path)
    for k, v in prof.get("settings", {}).items():
        put(k, v, "region profile '%s'" % s.profile.get("name", os.path.basename(prof_path)))

    for sec, body in doc.items():
        if sec == "advanced":
            continue
        ok = set(ENGINE_KEYS.get(sec, {})) | set(FILE_KEYS.get(sec, ())) | set(OTHER_KEYS.get(sec, ()))
        for k, v in body.items():
            if k not in ok:
                raise ConfigError("%s: unknown key '%s' in [%s]%s. Keys: %s" % (os.path.basename(path), k, sec, _suggest(k, ok), ", ".join(sorted(ok))))
            if k in ENGINE_KEYS.get(sec, {}):
                if v == "" or v is None:
                    continue
                if k == "lonlat" and not v:
                    continue
                if k in ("manual_apexes", "extra_feeders"):
                    v = [tuple(row) for row in v]
                put(ENGINE_KEYS[sec][k], v, "site.toml [%s] %s" % (sec, k))
            elif k in FILE_KEYS.get(sec, ()):
                if v:
                    s.files[k] = v if os.path.isabs(str(v)) else os.path.normpath(os.path.join(s.root, str(v)))
    for k, v in doc.get("advanced", {}).items():
        put(k, v, "site.toml [advanced]")

    st = doc.get("stages", {})
    s.stages = dict(risk=bool(st.get("risk", True)), hecras_handoff=bool(st.get("hecras_handoff", False)))
    P["RUN_AHP"], P["RUN_RAS_HANDOFF"] = s.stages["risk"], s.stages["hecras_handoff"]

    tg = doc.get("targets", {})
    s.target_mode, s.target_ids = str(tg.get("mode", "all")).lower(), [int(i) for i in tg.get("ids", [])]
    if s.target_mode not in ("all", "top", "ids"):
        raise ConfigError("[targets] mode must be 'all', 'top' or 'ids'")
    if s.target_mode == "ids" and not s.target_ids:
        raise ConfigError("[targets] mode = 'ids' needs ids = [fan complex numbers] (run once with mode = 'all' to see them)")
    P["AHP_TARGET"] = s.target_ids if s.target_mode == "ids" else s.target_mode

    ro = doc.get("runoff", {})
    s.cn_table, s.default_soil_group = ro.get("cn_table"), ro.get("default_soil_group")
    ah = doc.get("ahp", {})
    s.catchment_matrix, s.fan_matrix = ah.get("catchment_matrix"), ah.get("fan_matrix")

    s.hecras = dict(HECRAS_DEFAULTS)
    s.hecras.update({k: v for k, v in doc.get("hecras", {}).items() if k in OTHER_KEYS["hecras"]})
    s.score = dict(SCORE_DEFAULTS)
    s.score.update({k: v for k, v in doc.get("score", {}).items() if k in OTHER_KEYS["score"]})
    cx = int(s.hecras.get("complex") or 0)
    P["RAS_COMPLEX"] = cx if cx > 0 else None

    case_path = _find("cases", site.get("validation_case") or "", s.root)
    if case_path:
        c = load_toml(case_path)
        s.case = dict(c.get("case", {}), path=case_path, observed_hydrograph=c.get("observed_hydrograph"), probes=c.get("probes", []),
                      widths=c.get("widths", []))
        for k, v in (c.get("files") or {}).items():
            if v and k in FILE_KEYS["validation"] and k not in s.files:
                s.files[k] = v if os.path.isabs(v) else os.path.normpath(os.path.join(s.root, v))
        if "tolerance" in s.case and "tolerance" not in doc.get("score", {}):
            s.score["tolerance"] = float(s.case["tolerance"])
    inflow = str(s.hecras.get("inflow", "pipeline")).lower()
    if inflow not in ("pipeline", "observed"):
        raise ConfigError("[hecras] inflow must be 'pipeline' (the storm in [storm]) or 'observed' (the hydrograph of the validation case)")
    if inflow == "observed":
        if not (s.case and s.case.get("observed_hydrograph")):
            raise ConfigError("[hecras] inflow = 'observed' needs a validation case with an [observed_hydrograph] ([site] validation_case)")
        P["RAS_OBSERVED_HYDROGRAPH"] = dict(s.case["observed_hydrograph"])
        s.sources["RAS_OBSERVED_HYDROGRAPH"] = "validation case '%s'" % os.path.basename(case_path)
    s.P = P
    _validate(s, doc)
    return s


def _validate(s, doc):
    P, f = s.P, s.files
    if not f.get("dem"):
        raise ConfigError("[site] dem is missing: give the path of your DEM (GeoTIFF)")
    for k, p in f.items():
        if k in ("runs_dir", "project"):                 # made later ('fanrisk hecras fit'); the HEC-RAS commands check it themselves
            continue
        if not os.path.exists(p):
            raise ConfigError("file not found for '%s': %s" % (k, p))
    rain = f.get("rainfall_raster") or P.get("RAINFALL_DESIGN_MM")
    if rain:
        if "STORM_DURATION_MIN" not in s.sources:
            raise ConfigError("[storm] duration_min is missing: a rainfall depth needs the storm duration (minutes) to become a peak flow")
        if float(P["STORM_DURATION_MIN"]) <= 0:
            raise ConfigError("[storm] duration_min must be > 0")
    m = str(P["RUNOFF_METHOD"]).lower()
    if m not in ("coefficient", "scs_cn"):
        raise ConfigError("[runoff] method must be 'coefficient' or 'scs_cn'")
    if m == "scs_cn" and rain:
        has = P.get("CURVE_NUMBER") is not None or f.get("curve_number_raster") or (f.get("landcover_raster") and f.get("soil_group_raster"))
        if not has:
            raise ConfigError("runoff method 'scs_cn' (from %s) needs a curve number: [runoff] curve_number = <number>, or curve_number_raster, or "
                              "landcover_raster + soil_group_raster" % s.sources.get("RUNOFF_METHOD", "the settings"))
    if bool(f.get("landcover_raster")) != bool(f.get("soil_group_raster")):
        raise ConfigError("[runoff] landcover_raster and soil_group_raster are used together: give both or neither")
    if m == "coefficient" and not (0.0 < float(P["RUNOFF_COEFF"]) <= 1.0):
        raise ConfigError("[runoff] coefficient must be between 0 and 1")
    if P.get("CURVE_NUMBER") is not None and not (1 <= float(P["CURVE_NUMBER"]) <= 100):
        raise ConfigError("[runoff] curve_number must be between 1 and 100")
    b = list(P["RISK_CLASS_BREAKS"])
    if len(b) != 4 or any(b[i] >= b[i + 1] for i in range(3)) or b[0] <= 0 or b[-1] >= 1:
        raise ConfigError("[ahp] class_breaks must be 4 rising numbers between 0 and 1")
    for nm, mat, n in (("catchment_matrix", s.catchment_matrix, 7), ("fan_matrix", s.fan_matrix, 6)):
        if mat is None:
            continue
        if len(mat) != n or any(len(r) != n for r in mat):
            raise ConfigError("[ahp] %s must be %d x %d" % (nm, n, n))
        for i in range(n):
            for j in range(n):
                if mat[i][j] <= 0 or abs(mat[i][j] * mat[j][i] - 1.0) > 0.02:
                    raise ConfigError("[ahp] %s is not reciprocal at row %d, column %d (a[i][j] x a[j][i] must be 1)" % (nm, i + 1, j + 1))
    if s.stages["hecras_handoff"] and not s.stages["risk"]:
        raise ConfigError("[stages] hecras_handoff needs risk = true")
    if s.stages["hecras_handoff"] and not rain and "RAS_OBSERVED_HYDROGRAPH" not in s.sources:
        s.notes.append("no rainfall and no observed hydrograph: the HEC-RAS handoff will have geometry but no inflow hydrograph")
    if m == "coefficient" and float(P["RUNOFF_COEFF"]) >= 1.0 and rain and s.sources.get("RUNOFF_COEFF", "").startswith("region profile 'generic'"):
        s.notes.append("runoff: all rain runs off (generic profile, worst case). Describe the ground in [runoff] or pick a region profile.")
    if not rain and s.stages["risk"]:
        s.notes.append("no rainfall given: the maps show susceptibility with the storm trigger at its worst case (1.0), and no peak flows")


def settings_report(s):
    """Lines listing every engine setting that differs from the method defaults, and where it came from."""
    rows = []
    for k in sorted(s.sources):
        if k in ENGINE_DEFAULTS and s.P[k] != ENGINE_DEFAULTS[k]:
            rows.append("  %-26s = %-22s <- %s" % (k, s.P[k], s.sources[k]))
        elif k not in ENGINE_DEFAULTS:
            rows.append("  %-26s = %-22s <- %s" % (k, s.P.get(k), s.sources[k]))
    return rows
