# -*- coding: utf-8 -*-
"""DEM -> fans and catchments -> AHP risk maps -> HEC-RAS 2D handoff, driven by one Site (see config.py). Files in, files out: no QGIS needed.

Output folder:
    fans/      fanrisk.gpkg (all vector layers; GeoJSON files without geopandas) and the terrain-analysis rasters
    risk/      ahp_C<id>_*.tif (+ .qml colour styles), master_table.csv, ahp_summary.json, ahp_validation.csv
    hecras/    handoff_C<id>_<feeder>/ (terrain, flow area, boundary lines, hydrographs, scenario table, README)
    report/    run_log.txt, settings_used.json, run_summary.json
"""
import csv
import json
import math
import os
import re
import time

import numpy as np

from . import __version__, engine as E, gio, hydrology
from .config import settings_report

NODATA = gio.NODATA


class Logger(object):
    def __init__(self, path=None, echo=True):
        self.t0, self.lines, self.path, self.echo = time.time(), [], path, echo
        if path:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w", encoding="utf-8").close()

    def __call__(self, msg):
        line = "[%6.1fs] %s" % (time.time() - self.t0, msg)
        self.lines.append(line)
        if self.echo:
            try:
                print(line, flush=True)
            except UnicodeEncodeError:
                print(line.encode("ascii", "replace").decode("ascii"), flush=True)
        if self.path:
            with open(self.path, "a", encoding="utf-8") as fh:
                fh.write(line + "\n")


def _num(v):
    return None if v is None or v != v else float(v)


def write_qml(path_tif, breaks):
    """Colour style for a 0-1 risk raster (QGIS picks up <name>.qml next to the GeoTIFF automatically)."""
    cols = ["#1a9641", "#a6d96a", "#ffffbf", "#fdae61", "#d7191c"]
    labels = ["very low", "low", "moderate", "high", "very high"]
    ups = list(breaks) + ["inf"]
    items = "\n".join('     <item alpha="255" value="%s" label="%s" color="%s"/>' % (u, l, c) for u, l, c in zip(ups, labels, cols))
    txt = ("<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n<qgis version=\"3.34.0\" styleCategories=\"Symbology\">\n <pipe>\n"
           "  <rasterrenderer type=\"singlebandpseudocolor\" band=\"1\" opacity=\"0.8\" alphaBand=\"-1\" classificationMin=\"0\" classificationMax=\"1\">\n"
           "   <rasterTransparency/>\n   <rastershader>\n    <colorrampshader colorRampType=\"DISCRETE\" classificationMode=\"2\" clip=\"0\" "
           "minimumValue=\"0\" maximumValue=\"1\" labelPrecision=\"2\">\n" + items + "\n    </colorrampshader>\n   </rastershader>\n  </rasterrenderer>\n"
           "  <brightnesscontrast brightness=\"0\" contrast=\"0\" gamma=\"1\"/>\n  <rasterresampler maxOversampling=\"2\"/>\n </pipe>\n <blendMode>0</blendMode>\n</qgis>\n")
    with open(os.path.splitext(path_tif)[0] + ".qml", "w") as fh:
        fh.write(txt)


# =============================================================================
# inputs
# =============================================================================
def load_inputs(site, log):
    """DEM and every optional layer, on the DEM grid."""
    P = site.P
    arr, gt, crs, nd = gio.load_dem(site.file("dem"), P["REPROJECT_GEOGRAPHIC"], log)
    dem, note = E.clean_dem(arr, nd, P["NODATA_VALUES"], P["ZERO_IS_NODATA"])
    log("DEM: %s  (%s, %d x %d cells of %.1f m)" % (os.path.basename(site.file("dem")), crs.authid() or "no EPSG code", dem.shape[1], dem.shape[0], abs(gt[1])))
    if note:
        log("  DEM cleaning: " + note)
    if dem.size > 6e6:
        log("WARNING: %.1f million cells - the hydrology will be slow and memory hungry. Consider clipping the DEM to the study area." % (dem.size / 1e6))
    I = dict(dem=dem, gt=gt, crs=crs)

    def to_dem(lon, lat):
        x, y = gio.transform_xy([lon], [lat], gio.WGS84, crs)
        return float(x[0]), float(y[0])
    I["to_dem"] = to_dem
    I["manual"] = []
    for row in P["MANUAL_APEXES"]:
        fid, system, lon, lat = row
        x, y = to_dem(lon, lat)
        I["manual"].append(dict(name=str(fid), system=str(system), x=x, y=y))
    I["reference"] = None
    if site.file("reference_fans"):
        feats = gio.read_vector(site.file("reference_fans"), crs)
        I["reference"] = [dict(name="ref%d" % (i + 1), mask=gio.rasterize([f["geometry"]], dem.shape, gt)) for i, f in enumerate(feats)]
        log("Reference fans: %d polygon(s)" % len(I["reference"]))
    return I


def load_risk_layers(site, I, out_dir, log):
    P, dem, gt, crs = site.P, I["dem"], I["gt"], I["crs"]
    L = dict(rain=None, hist=None, bare=None, built=None, lakes=[], rain_design=False, cn=None)
    for key, fk in (("rain", "rainfall_raster"), ("hist", "historical_rainfall"), ("bare", "bare_fraction"), ("built", "built_density")):
        if site.file(fk):
            L[key] = gio.raster_on_grid(site.file(fk), dem.shape, gt, crs)
            log("  input layer %-8s <- %s" % (key, os.path.basename(site.file(fk))))
    if L["rain"] is None and P["RAINFALL_DESIGN_MM"]:
        L["rain"] = np.where(np.isfinite(dem), float(P["RAINFALL_DESIGN_MM"]), np.nan)
        L["rain_design"] = True
        log("  design storm: %.0f mm in %.0f min everywhere" % (float(P["RAINFALL_DESIGN_MM"]), float(P["STORM_DURATION_MIN"])))
    if site.file("glacial_lakes"):
        for f in gio.read_vector(site.file("glacial_lakes"), crs):
            g = f["geometry"]
            v = {str(k).lower(): val for k, val in f["properties"].items()}.get("volume_m3")
            L["lakes"].append(dict(area_m2=gio.geom_area(g), volume_m3=float(v) if isinstance(v, (int, float)) and v > 0 else None,
                                   cells=np.flatnonzero(gio.rasterize([g], dem.shape, gt).ravel())))
        log("  glacial lakes: %d polygon(s)" % len(L["lakes"]))
    if str(P["RUNOFF_METHOD"]).lower() == "scs_cn" and L["rain"] is not None:
        if site.file("curve_number_raster"):
            L["cn"] = gio.raster_on_grid(site.file("curve_number_raster"), dem.shape, gt, crs, "nearest")
            src = os.path.basename(site.file("curve_number_raster"))
        elif site.file("landcover_raster"):
            lc = gio.raster_on_grid(site.file("landcover_raster"), dem.shape, gt, crs, "nearest")
            sg = gio.raster_on_grid(site.file("soil_group_raster"), dem.shape, gt, crs, "nearest")
            L["cn"] = hydrology.cn_from_landcover(lc, sg, site.cn_table, site.default_soil_group)
            src = "land cover x soil group"
        if L["cn"] is not None:
            valid = np.isfinite(dem)
            log("  runoff: SCS curve number from %s (%s); antecedent moisture %s" % (src, hydrology.describe_cn(L["cn"], valid), P["CN_AMC"]))
            if P.get("CURVE_NUMBER") is not None:        # fill gaps of the raster with the single number
                L["cn"] = np.where(np.isfinite(L["cn"]), L["cn"], float(P["CURVE_NUMBER"]))
            gio.write_raster(os.path.join(out_dir, "curve_number.tif"), np.where(valid, L["cn"], np.nan), gt, crs, "float32", NODATA)
        else:
            log("  runoff: SCS curve number %s everywhere; antecedent moisture %s" % (P["CURVE_NUMBER"], P["CN_AMC"]))
    elif L["rain"] is not None:
        log("  runoff: coefficient %.2f (share of the rain that runs off)" % float(P["RUNOFF_COEFF"]))
    return L


def load_validation(site, I, log):
    dem, gt, crs = I["dem"], I["gt"], I["crs"]
    if not (site.file("flood_footprint") or site.file("flood_points")):
        return None
    presence = np.zeros(dem.shape, dtype=bool)
    for fk, kind in (("flood_footprint", "footprint"), ("flood_points", "points")):
        if site.file(fk):
            feats = gio.read_vector(site.file(fk), crs)
            presence |= gio.rasterize([f["geometry"] for f in feats], dem.shape, gt)
            log("  validation %-9s <- %s: %d feature(s)" % (kind, os.path.basename(site.file(fk)), len(feats)))
    log("  validation: %d flooded cell(s) on the DEM grid (%.2f km2)" % (presence.sum(), presence.sum() * abs(gt[1] * gt[5]) / 1e6))
    return dict(presence=presence)


# =============================================================================
# outputs
# =============================================================================
def fan_outputs(res, P, I, out_dir, log):
    """Terrain-analysis rasters and the fan / catchment vector layers. Returns the layer list (written later, with the AHP layers)."""
    dem, gt, crs = I["dem"], I["gt"], I["crs"]
    grid, fans, cplx, cands = res["grid"], res["fans"], res["complexes"], res["candidates"]
    paper = P["APEX_METHOD"].lower() in ("paper", "both")
    rasters = [("flowacc_km2", res["acc_km2"], "float32", NODATA), ("streams", res["streams"].astype(np.uint8), "uint8", 0),
               ("fan_cones", res["cone"], "float32", NODATA), ("fan_id", res["fan_id"], "int32", 0), ("complex_id", res["complex_id"], "int32", 0)]
    if P["WRITE_DEBUG_RASTERS"]:
        rasters += [("fan_surface", res["zsurf"], "float32", NODATA), ("slope_deg", res["slope"], "float32", NODATA)]
        if paper:
            rasters += [("slope_smoothed", res["slope_s"], "float32", NODATA), ("steep_zone", res["steep"].astype(np.uint8), "uint8", 0),
                        ("feeder_channels", res["feeder"].astype(np.uint8), "uint8", 0), ("curvature", res["curv"], "float32", NODATA)]
    for name, arr, dt, nd in rasters:
        gio.write_raster(os.path.join(out_dir, name + ".tif"), arr, gt, crs, dt, nd)
    x0, y0, dx, dy = grid.x0, grid.y0, grid.dx, grid.dy

    def poly(mask, r0, c0):
        return gio.polygonize(mask, (x0 + c0 * dx, dx, 0.0, y0 - r0 * dy, 0.0, -dy))

    def smoothed(g):
        return gio.smooth_geometry(g) if (g is not None and P["SMOOTH_OUTLINES"]) else g

    def cells_poly(cells):
        cr, cc = cells // grid.W, cells % grid.W
        a0, b0 = cr.min(), cc.min()
        m = np.zeros((cr.max() + 1 - a0, cc.max() + 1 - b0), dtype=bool)
        m[cr - a0, cc - b0] = True
        return poly(m, a0, b0)

    pt = lambda x, y: {"type": "Point", "coordinates": [float(x), float(y)]}
    rows = dict(cplx=[], fan=[], ccatch=[], catch=[], drain=[], prof=[], end=[], apex=[])
    for cx in cplx:
        r0, r1, c0, c1 = cx["win"]
        g = smoothed(poly(cx["mask"], r0, c0))
        if g is not None:
            rows["cplx"].append((g, [cx["cid"], len(cx["members"]), cx["fan_ids"], cx["area_km2"], cx["catch_area_km2"], cx["catch_relief_m"], cx["melton"],
                                     int(cx["catch_complete"]), cx["flags"], _num(cx["iou_ref"])]))
        gc = cells_poly(cx["catch_cells"])
        if gc is not None:
            rows["ccatch"].append((gc, [cx["cid"], len(cx["members"]), cx["catch_area_km2"], cx["catch_relief_m"], cx["melton"], int(cx["catch_complete"])]))
    for fan in fans:
        cand, ra = fan["cand"], fan["ra"]
        r0, r1, c0, c1 = fan["win"]
        g = smoothed(poly(fan["mask"], r0, c0))
        if g is not None:
            rows["fan"].append((g, [fan["fid"], fan["complex_id"], cand["name"], cand["system"], cand["source"], cand["method"], fan["area_km2"],
                                    fan["length_m"], fan["span_deg"], fan["cone_slope_deg"], fan["mean_slope"], ra["z_apex"], fan["toe_elev_min"],
                                    int(ra["n_trunc"]), fan["catch_area_km2"], fan["catch_relief_m"], fan["melton"], int(fan["catch_complete"]),
                                    fan["drain_len_km"], fan["flags"], _num(fan["iou_ref"])]))
        gc = cells_poly(fan["catch_cells"])
        if gc is not None:
            rows["catch"].append((gc, [fan["fid"], fan["complex_id"], cand["name"], fan["catch_area_km2"], fan["catch_relief_m"], fan["melton"],
                                       int(fan["catch_complete"])]))
        x1, y1, x2, y2 = fan["drain_segments"]
        if x1.size:
            gd = {"type": "MultiLineString", "coordinates": [[[round(float(a), 2), round(float(b), 2)], [round(float(c), 2), round(float(d), 2)]]
                                                             for a, b, c, d in zip(x1, y1, x2, y2)]}
            rows["drain"].append((gd, [fan["fid"], cand["name"], fan["drain_len_km"]]))
        for k in range(ra["D"].size):
            if ra["D"][k] <= 0:
                continue
            th = ra["ang"][k]
            ex, ey = cand["x"] + ra["D"][k] * math.cos(th), cand["y"] + ra["D"][k] * math.sin(th)
            brg = (90.0 - math.degrees(th)) % 360.0
            rows["prof"].append(({"type": "LineString", "coordinates": [[float(cand["x"]), float(cand["y"])], [float(ex), float(ey)]]},
                                 [fan["fid"], brg, float(ra["D"][k]), float(ra["z_apex"] - ra["Zend"][k]), ra["reason"][k]]))
            rows["end"].append((pt(ex, ey), [fan["fid"], brg, float(ra["D"][k]), float(ra["Zend"][k]), ra["reason"][k]]))
    for cand in cands:
        if "x" not in cand:
            cand["x"], cand["y"] = grid.rc_to_xy(cand["row"], cand["col"])
            cand["z"] = float(dem[cand["row"], cand["col"]])
            cand["acc_km2"] = float(res["ctx"].acc[cand["cell"]] * res["ctx"].cell_km2)
    lon = lat = np.full(len(cands), np.nan)
    if cands and (crs.epsg or crs.wkt):
        try:
            lon, lat = gio.transform_xy([c["x"] for c in cands], [c["y"] for c in cands], crs, gio.WGS84)
        except Exception as e:
            log("  (apex longitude / latitude not written: %s)" % e)
    for k, cand in enumerate(cands):
        rows["apex"].append((pt(cand["x"], cand["y"]), [int(cand.get("fan_id", 0)), cand["name"], cand["system"], cand["source"], cand["method"],
                                                         cand["status"], _num(lon[k]), _num(lat[k]), cand["z"], cand["acc_km2"], _num(cand["mean_curv"]),
                                                         _num(cand["feeder_len_m"])]))
    I_, D, S = "integer", "double", "string"
    return [
        dict(name="fan_complexes", rows=rows["cplx"], fields=[("complex_id", I_), ("n_fans", I_), ("fan_ids", S), ("area_km2", D), ("catch_km2", D),
                                                              ("relief_m", D), ("melton", D), ("catch_full", I_), ("flags", S), ("iou_ref", D)]),
        dict(name="fans", rows=rows["fan"], fields=[("fan_id", I_), ("complex_id", I_), ("name", S), ("system", S), ("source", S), ("method", S),
                                                    ("area_km2", D), ("length_m", D), ("span_deg", D), ("cone_slope", D), ("mean_slope", D),
                                                    ("apex_elev", D), ("toe_elev", D), ("n_trunc", I_), ("catch_km2", D), ("relief_m", D), ("melton", D),
                                                    ("catch_full", I_), ("drain_km", D), ("flags", S), ("iou_ref", D)]),
        dict(name="apexes", rows=rows["apex"], fields=[("fan_id", I_), ("name", S), ("system", S), ("source", S), ("method", S), ("status", S),
                                                       ("lon", D), ("lat", D), ("elev_m", D), ("acc_km2", D), ("curv_mean", D), ("feeder_m", D)]),
        dict(name="fan_end_points", rows=rows["end"], fields=[("fan_id", I_), ("bearing", D), ("dist_m", D), ("elev_m", D), ("end_reason", S)]),
        dict(name="radial_profiles", rows=rows["prof"], fields=[("fan_id", I_), ("bearing", D), ("length_m", D), ("relief_m", D), ("end_reason", S)]),
        dict(name="source_drainage", rows=rows["drain"], fields=[("fan_id", I_), ("name", S), ("length_km", D)]),
        dict(name="complex_catchments", rows=rows["ccatch"], fields=[("complex_id", I_), ("n_fans", I_), ("area_km2", D), ("relief_m", D), ("melton", D),
                                                                     ("complete", I_)]),
        dict(name="catchments", rows=rows["catch"], fields=[("fan_id", I_), ("complex_id", I_), ("name", S), ("area_km2", D), ("relief_m", D),
                                                            ("melton", D), ("complete", I_)]),
    ]


def risk_outputs(ahp, P, I, out_dir, log):
    """AHP rasters (+ styles), tables, summary. Returns the feeder / inflow-line layers."""
    dem, gt, crs = I["dem"], I["gt"], I["crs"]
    wgt = lambda r0, c0: (gt[0] + c0 * gt[1], gt[1], 0.0, gt[3] + r0 * gt[5], 0.0, gt[5])
    per_feeder = P["SCENARIO_MAPS"] is True or (P["SCENARIO_MAPS"] == "auto" and len(ahp["targets"]) <= 2)
    sc = {"region": "abs", "fan": "rel"}                       # abs = absolute physical scales, rel = relative inside the fan
    for t in ahp["targets"]:
        cid = t["cid"]
        r0, r1, c0, c1 = t["win"]
        items = [("risk_abs", t["R_region"]), ("risk_rel", t["R_fan"]), ("hazard_abs", t["H_region"]), ("hazard_rel", t["H_fan"]), ("hand_m", t["hand"])]
        if t.get("flow_km2") is not None:
            items.append(("flow_km2", np.where(t["mask"], t["flow_km2"], np.nan)))
        items.append(("corridor", np.where(t["mask"], t["corridor"].astype(np.float32), np.nan)))
        if np.nanmin(t["E"]) < 1.0:
            items.append(("exposure", t["E"]))
        for (kind, scope), arr in sorted(t["worst"].items()):
            items.append(("worst_%s_%s" % (kind, sc[scope]), arr))
        if per_feeder:
            items += [("scen_%s_%s" % (s["feeder"], s["kind"]), s["arr"]) for s in t["scen"] if s["scope"] == "region"]
        for name, arr in items:
            path = os.path.join(out_dir, "ahp_C%d_%s.tif" % (cid, name))
            gio.write_raster(path, np.asarray(arr, dtype=np.float32), wgt(r0, c0), crs, "float32", NODATA)
            if name not in ("hand_m", "flow_km2", "corridor"):
                write_qml(path, P["RISK_CLASS_BREAKS"])
        R0, R1, C0, C1 = t["clip_win"]
        gio.write_raster(os.path.join(out_dir, "ahp_C%d_dem_clip.tif" % cid), dem[R0:R1, C0:C1].astype(np.float32), wgt(R0, C0), crs, "float32", NODATA)
    rf, ri = [], []
    for t in ahp["targets"]:
        cmap = {f["name"]: f for f in t["feeders"]}
        for f in t["feeders"]:
            m, h = f["metrics"], f["hydro"] or {}
            rf.append(({"type": "Point", "coordinates": [float(f["x"]), float(f["y"])]},
                       [int(t["cid"]), f["name"], f["source"], int(f["fan_id"]), m["area_km2"], m["relief_m"], m["mean_slope_deg"], m["melton"],
                        m["drainage_density"], m["main_channel_km"], _num(m["kirpich_tc_min"]), f["C"], _num(f["rain_mm"]), _num(h.get("curve_number")),
                        _num(h.get("peak_m3s")), _num(h.get("runoff_mm")), _num(h.get("volume_m3")), _num(h.get("tp_min")), f["T"], int(f["lakes"]),
                        f["lake_volume_m3"], f["G"]]))
        for nm, x1, y1, x2, y2 in t["inflow"]:
            hq = cmap[nm]["hydro"] or {}
            ri.append(({"type": "LineString", "coordinates": [[float(x1), float(y1)], [float(x2), float(y2)]]},
                       [int(t["cid"]), nm, cmap[nm]["C"], _num(hq.get("peak_m3s")), _num(hq.get("volume_m3"))]))
    I_, D, S = "integer", "double", "string"
    layers = [dict(name="feeders", rows=rf, fields=[("complex_id", I_), ("name", S), ("source", S), ("fan_id", I_), ("catch_km2", D), ("relief_m", D),
                                                    ("slope_deg", D), ("melton", D), ("drain_dens", D), ("main_km", D), ("tc_min", D), ("C_source", D),
                                                    ("rain_mm", D), ("curve_no", D), ("peak_q_m3s", D), ("runoff_mm", D), ("volume_m3", D), ("tp_min", D),
                                                    ("T_trigger", D), ("lakes", I_), ("lake_m3", D), ("G_glof", D)]),
              dict(name="inflow_lines", rows=ri, fields=[("complex_id", I_), ("feeder", S), ("C_source", D), ("peak_q_m3s", D), ("volume_m3", D)])]
    keys = []
    for row in ahp["master"]:
        keys += [k for k in row if k not in keys]
    with open(os.path.join(out_dir, "master_table.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for row in ahp["master"]:
            w.writerow({k: (None if (isinstance(v, float) and v != v) else v) for k, v in row.items()})
    vrows = []
    for t in ahp["targets"]:
        for v in t.get("validation", []):
            r = dict(complex_id=t["cid"], map=v["map"], flooded_cells=v["flooded_cells"], covered_pct=v["covered_pct"],
                     flooded_share_of_fan_pct=v["flooded_share_pct"], auc_roc=v["auc_roc"], auc_success_rate=v["auc_success"],
                     flooded_cells_sit_in_top_pct=v["mean_rank_pct"])
            for k, c in enumerate(v["classes"]):
                r["class%d_area_pct" % (k + 1)], r["class%d_flooded_pct" % (k + 1)] = c["area_pct"], c["flooded_pct"]
            vrows.append(r)
    if vrows:
        with open(os.path.join(out_dir, "ahp_validation.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(vrows[0].keys()))
            w.writeheader()
            for r in vrows:
                w.writerow({k: (None if (isinstance(v, float) and v != v) else v) for k, v in r.items()})
    wc, ic = ahp["weights"]["catchment"]
    wf_, icf = ahp["weights"]["fan"]
    summary = dict(weights_catchment=wc, weights_fan=wf_, consistency_ratio=dict(catchment=ic["cr"], fan=icf["cr"]),
                   dropped_criteria=dict(catchment=ic["dropped"], fan=icf["dropped"]),
                   targets=[dict(complex_id=t["cid"], class_area_km2=dict(absolute=t["classes"]["region"], relative=t["classes"]["fan"]),
                                 hand_sensitivity=t["sens"], validation=t.get("validation", []),
                                 feeders=[dict(name=f["name"], catchment_km2=f["metrics"]["area_km2"], C_source=f["C"], T_trigger=f["T"], G_glof=f["G"],
                                               event=f["hydro"]) for f in t["feeders"]]) for t in ahp["targets"]],
                   class_breaks=list(P["RISK_CLASS_BREAKS"]),
                   note="Scores are relative susceptibility (0-1), not calibrated probabilities or inundation depths. "
                        "'abs' maps use fixed physical scales (comparable between fans); 'rel' maps are normalised inside each fan.")
    with open(os.path.join(out_dir, "ahp_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    log("Risk files -> %s: risk_abs / risk_rel / hazard / flow_km2 / worst-case%s GeoTIFFs (+ .qml styles), master_table.csv, ahp_summary.json%s" % (
        out_dir, " / per-feeder scenario" if per_feeder else "", ", ahp_validation.csv" if vrows else ""))
    return layers


def handoff_outputs(res, ahp, site, I, log):
    """The HEC-RAS 2D package for one flood path. Returns its folder."""
    P, crs = site.P, I["crs"]
    line = None
    if site.file("centerline"):
        line = gio.first_line(gio.read_vector(site.file("centerline"), crs))
        if line is None:
            raise RuntimeError("no usable line in %s" % site.file("centerline"))
    pkg = E.ras_build(res, ahp, P, line_xy=line, log=log)
    folder = site.handoff_dir(pkg["cid"], pkg["feeder"]["name"])
    os.makedirs(folder, exist_ok=True)
    for fn, tx in E.ras_texts(pkg, P, crs.epsg).items():
        with open(os.path.join(folder, fn), "w") as fh:
            fh.write(tx)
    gio.write_raster(os.path.join(folder, "terrain.tif"), pkg["terrain"], pkg["gt_window"], crs, "float32", NODATA)
    burned_tif = os.path.join(folder, "terrain_burned.tif")
    if pkg["terrain_burned"] is not None:
        gio.write_raster(burned_tif, pkg["terrain_burned"], pkg["gt_window"], crs, "float32", NODATA)
        cut = pkg["terrain"] - pkg["terrain_burned"]
        cut = cut[np.isfinite(cut) & (cut > 0.01)]
        log("channel burned into terrain_burned.tif: %d cells, %.0f m wide strip lowered to %.1f m below the bed profile (median cut %.1f m)" % (
            cut.size, P["RAS_CHANNEL_WIDTH_M"], P["RAS_BURN_DEPTH_M"], float(np.median(cut)) if cut.size else 0.0))
    elif os.path.exists(burned_tif):                      # a burned terrain from an earlier run would be picked up by 'hecras fit'
        try:
            os.remove(burned_tif)
        except OSError:
            os.replace(burned_tif, burned_tif + ".old")
        log("burn_depth_m = 0: removed the earlier terrain_burned.tif so 'hecras fit' uses terrain.tif")
    gio.write_raster(os.path.join(folder, "landcover.tif"), pkg["landcover"], pkg["gt_window"], crs, "uint8", 0)
    r = pkg["reach"]
    if site.file("fine_dem"):
        try:
            _a, g_f, c_f, _n = gio.read_raster_raw(site.file("fine_dem"))
            res_f = abs(g_f[1]) if (c_f.same(crs) or not c_f.geographic and not c_f.epsg) else float(P["RAS_CELL_SIZE_M"])
            x0, x1, y0, y1 = r["ring"][:, 0].min() - 100, r["ring"][:, 0].max() + 100, r["ring"][:, 1].min() - 100, r["ring"][:, 1].max() + 100
            shp = (int(math.ceil((y1 - y0) / res_f)), int(math.ceil((x1 - x0) / res_f)))
            fine = gio.raster_on_grid(site.file("fine_dem"), shp, (x0, res_f, 0.0, y1, 0.0, -res_f), crs)
            gio.write_raster(os.path.join(folder, "terrain_fine.tif"), fine, (x0, res_f, 0.0, y1, 0.0, -res_f), crs, "float32", NODATA)
            log("  terrain_fine.tif written from %s at %.1f m" % (os.path.basename(site.file("fine_dem")), res_f))
        except Exception as e:
            log("WARNING: the fine DEM could not be clipped (%s: %s); terrain.tif uses the pipeline DEM." % (type(e).__name__, e))
    ring = [[float(x), float(y)] for x, y in r["ring"]]
    shapes = [("flow_area_polygon.shp", {"type": "Polygon", "coordinates": [ring]}, "2D flow area"),
              ("centerline.shp", {"type": "LineString", "coordinates": [[float(a), float(b)] for a, b in zip(r["x"], r["y"])]}, "centreline"),
              ("bc_inflow.shp", {"type": "LineString", "coordinates": [[float(a), float(b)] for a, b in r["inflow"]]}, "inflow"),
              ("bc_outflow.shp", {"type": "LineString", "coordinates": [[float(a), float(b)] for a, b in r["outflow"]]}, "outflow")]
    try:
        if not all(gio.write_shapefile(os.path.join(folder, fn), g, crs, lab) for fn, g, lab in shapes):
            log("  (shapefiles not written - geopandas is not installed; the .geojson copies hold the same lines)")
    except Exception as e:
        log("WARNING: shapefiles could not be written (%s: %s); use the .geojson copies." % (type(e).__name__, e))
    hy = pkg["hydro"]
    log("HEC-RAS handoff -> %s" % folder)
    log("  C%d %s: reach %.2f km%s, 2D flow area %.2f km2 (about %.0f cells at %.0f m), bed %.0f -> %.0f m" % (
        pkg["cid"], pkg["feeder"]["name"], r["length_m"] / 1000, " (your line)" if line else " (DEM flow path)", r["area_km2"],
        pkg["mesh"][10.0] * (10.0 / P["RAS_CELL_SIZE_M"]) ** 2, P["RAS_CELL_SIZE_M"], r["z"][0], r["z"][-1]))
    if hy:
        log("  inflow: %s; peak %.0f m3/s at %.0f min, volume %.2f million m3. Flow type %s (Melton %.2f). %d scenario runs listed." % (
            hy["source"], hy["peak"], hy["tp_min"], hy["volume_m3"] / 1e6, pkg["flow_type"], pkg["melton"], len(pkg["scenarios"])))
    for n in pkg["notes"]:
        log("  NOTE: " + n)
    return folder


# =============================================================================
# the run
# =============================================================================
def run(site, echo=True):
    """Run the stages switched on in the site file. Returns a summary dict (also written to report/run_summary.json)."""
    P = site.P
    out = site.out_dir
    d = {k: os.path.join(out, k) for k in ("fans", "risk", "hecras", "report")}
    for k in ("fans", "report"):
        os.makedirs(d[k], exist_ok=True)
    log = Logger(os.path.join(d["report"], "run_log.txt"), echo)
    log("fanrisk %s - site '%s' - file backend: %s" % (__version__, site.name, gio.backend()))
    log("Region profile: %s - %s" % (site.profile.get("name"), site.profile.get("description", "")))
    if site.case:
        log("Validation case: %s (used for testing only)" % site.case.get("name"))
    rep = settings_report(site)
    if rep:
        log("Settings that differ from the method defaults:")
        for line in rep:
            log(line)
    for n in site.notes:
        log("NOTE: " + n)
    with open(os.path.join(d["report"], "settings_used.json"), "w") as fh:
        json.dump(dict(fanrisk=__version__, site=site.name, site_file=site.path, region_profile=site.profile, validation_case=(site.case or {}).get("name"),
                       files=site.files, stages=site.stages, settings=P, setting_sources=site.sources), fh, indent=2, default=str)

    I = load_inputs(site, log)
    res = E.run_pipeline(I["dem"], I["gt"], P, manual=I["manual"], reference=I["reference"], log=log)
    fans, cplx = res["fans"], res["complexes"]
    summary = dict(site=site.name, fanrisk=__version__, n_fans=len(fans), n_complexes=len(cplx), output_dir=out, risk=None, handoff=None, errors=[],
                   complexes=[dict(id=cx["cid"], n_fans=len(cx["members"]), area_km2=cx["area_km2"], catchment_km2=cx["catch_area_km2"],
                                   melton=cx["melton"], iou_ref=cx["iou_ref"]) for cx in cplx])
    ahp, ahp_layers = None, []
    if P["RUN_AHP"] and cplx:
        os.makedirs(d["risk"], exist_ok=True)
        old = (E.CATCH_MATRIX, E.FAN_MATRIX)
        try:
            if site.catchment_matrix is not None:
                E.CATCH_MATRIX = np.array(site.catchment_matrix, dtype=float)
                log("AHP: catchment pairwise matrix from site.toml")
            if site.fan_matrix is not None:
                E.FAN_MATRIX = np.array(site.fan_matrix, dtype=float)
                log("AHP: fan pairwise matrix from site.toml")
            L = load_risk_layers(site, I, d["risk"], log)
            txy = I["to_dem"](*P["AHP_TARGET_LONLAT"]) if P["AHP_TARGET_LONLAT"] else None
            extra = []
            if P["EXTRA_FEEDERS"]:
                ex = []
                for nm, lon, lat in P["EXTRA_FEEDERS"]:
                    x, y = I["to_dem"](lon, lat)
                    ex.append(dict(name=str(nm), system="", x=x, y=y))
                ctx = res["ctx"]
                extra = [dict(name=c["name"], cell=c["cell"]) for c in E.snap_manual(ex, ctx.grid, ctx.streams_flat, ctx.acc, P, log) if c["status"] == "candidate"]
            ahp = E.run_ahp(res, P, L, txy, extra, load_validation(site, I, log), log=log)
            for nm, (w_, info) in (("catchment", ahp["weights"]["catchment"]), ("fan", ahp["weights"]["fan"])) if ahp else ():
                if info["cr"] > 0.10:
                    log("WARNING: the %s pairwise matrix is inconsistent (CR %.3f > 0.10): revise the judgements." % (nm, info["cr"]))
        except Exception as e:                           # the fan results must survive a failure of the risk stage
            log("ERROR: the risk stage failed (%s: %s) - the fan results are unaffected." % (type(e).__name__, e))
            summary["errors"].append("risk: %s: %s" % (type(e).__name__, e))
            ahp = None
        finally:
            E.CATCH_MATRIX, E.FAN_MATRIX = old

    log("Writing fan outputs to %s" % d["fans"])
    layers = fan_outputs(res, P, I, d["fans"], log)
    if ahp is not None and ahp["targets"]:
        try:
            ahp_layers = risk_outputs(ahp, P, I, d["risk"], log)
            summary["risk"] = dict(targets=[dict(complex_id=t["cid"], class_area_km2=t["classes"]["region"],
                                                 feeders=[dict(name=f["name"], catchment_km2=f["metrics"]["area_km2"], C_source=f["C"], T_trigger=f["T"],
                                                               rain_mm=f["rain_mm"], event=f["hydro"]) for f in t["feeders"]]) for t in ahp["targets"]])
        except Exception as e:
            log("ERROR: could not write the risk outputs (%s: %s)" % (type(e).__name__, e))
            summary["errors"].append("risk outputs: %s: %s" % (type(e).__name__, e))
        if P["RUN_RAS_HANDOFF"]:
            try:
                summary["handoff"] = handoff_outputs(res, ahp, site, I, log)
            except Exception as e:                       # never lose the fan / risk results for the handoff
                log("ERROR: the HEC-RAS handoff failed (%s: %s) - everything else is unaffected." % (type(e).__name__, e))
                summary["errors"].append("handoff: %s: %s" % (type(e).__name__, e))
    paths = gio.write_layers(d["fans"], "fanrisk.gpkg", layers + ahp_layers, I["crs"], log)
    log("Vector layers -> %s" % (paths[0] if len(paths) == 1 else os.path.dirname(paths[0]) if paths else "(none)"))

    log("RESULTS: %d fans in %d fan complexes" % (len(fans), len(cplx)))
    log("  FAN COMPLEXES (merged flood plains / bajadas) and ALL catchments feeding them")
    log("  %-4s %-6s %10s %11s %7s  %s" % ("id", "fans", "area km2", "catch km2", "Melton", "flags"))
    for cx in cplx[:30]:
        log("  %-4d %-6d %10.2f %11.1f %7.2f  %s" % (cx["cid"], len(cx["members"]), cx["area_km2"], cx["catch_area_km2"], cx["melton"],
                                                    cx["flags"] + ("" if cx["iou_ref"] is None else "  IoU(ref) %.2f" % cx["iou_ref"])))
    for c in res["candidates"]:
        if c["source"] == "manual" and c["status"] != "accepted":
            log("  MANUAL APEX %s was not accepted: %s" % (c["name"], c["status"]))
    if ahp is not None and ahp["targets"]:
        log("  FEEDERS (C = source susceptibility 0-1, T = storm trigger 0-1)")
        log("  %-9s %10s %6s %9s %11s %9s %6s" % ("feeder", "catch km2", "C", "rain mm", "peak m3/s", "tp min", "T"))
        for t in ahp["targets"]:
            for f in t["feeders"]:
                h = f["hydro"] or {}
                fm = lambda v, s_: (s_ % v) if v is not None else "-"
                log("  C%d-%-6s %10.2f %6.2f %9s %11s %9s %6.2f" % (t["cid"], f["name"], f["metrics"]["area_km2"], f["C"], fm(f["rain_mm"], "%.0f"),
                                                                    fm(h.get("peak_m3s"), "%.0f"), fm(h.get("tp_min"), "%.0f"), f["T"]))
    if summary["errors"]:
        log("FINISHED WITH ERRORS: " + " | ".join(summary["errors"]))
    else:
        log("Done. Open %s in QGIS: the risk GeoTIFFs load with their colours." % out)
    with open(os.path.join(d["report"], "run_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    summary["_res"], summary["_ahp"] = res, ahp
    return summary
