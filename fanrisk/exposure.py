# -*- coding: utf-8 -*-
"""Who is hit: OpenStreetMap buildings, roads and critical facilities overlaid on the flood maps.

    fanrisk exposure <site.toml> --osm data/map.osm [--buildings <ms_buildings_*.csv.gz>...] [--run R04 | --maps <folder>] [--out <folder>]

Buildings: Microsoft Global ML Building Footprints when given (or found as ms_buildings_*.csv.gz next to the DEM; get them with
tools/get_ms_buildings.py) - OSM maps only a fraction of the houses in Ladakh. OSM always supplies roads and critical facilities.

Layers it can use (each is optional, the ones that exist are used):
    AHP risk     runs/<site>/risk/ahp_C*_risk_abs.tif (all complexes; 'high' = RISK_CLASS_BREAKS[2] and above)
    HEC-RAS      the run's readback/hazard_class.tif (0 dry, 1 low, 2 moderate, 3 significant, 4 extreme) and Depth (Max)
    observed     [validation] flood_footprint (the 2010 footprint) - only inside its survey window
The OSM file is an .osm XML export (openstreetmap.org 'Export', the API map call, or Overpass with 'out body' and recursion).
Writes into <out>/: exposure_summary.json, buildings.csv, critical_facilities.csv, roads.csv, exposure_buildings.geojson,
exposure_report.txt and exposure_map.png. NOTE: OSM shows today's buildings, not those standing in 2010.
"""
import csv
import glob
import json
import os
import xml.etree.ElementTree as ET

import numpy as np

from . import gio

HAZ = ("dry", "low", "moderate", "significant", "extreme")
CRITICAL = {
    "amenity": {"school", "college", "university", "kindergarten", "hospital", "clinic", "doctors", "pharmacy", "police", "fire_station",
                "place_of_worship", "community_centre", "townhall", "marketplace", "bus_station", "fuel", "shelter", "social_facility"},
    "healthcare": None,                                   # any value
    "building": {"school", "hospital", "college", "university", "temple", "monastery", "mosque", "church", "government", "public"},
    "emergency": None,
    "office": {"government"},
}


def parse_osm(path):
    """-> nodes {id: (lon, lat)}, ways [dict(id, tags, nodes)], tagged nodes [dict(id, tags, lon, lat)]."""
    nodes, ways, pnodes = {}, [], []
    for _ev, el in ET.iterparse(path, events=("end",)):
        if el.tag == "node":
            lon, lat = float(el.get("lon")), float(el.get("lat"))
            nodes[el.get("id")] = (lon, lat)
            tags = {t.get("k"): t.get("v") for t in el.findall("tag")}
            if tags:
                pnodes.append(dict(id="n" + el.get("id"), tags=tags, lon=lon, lat=lat))
            el.clear()
        elif el.tag == "way":
            ways.append(dict(id="w" + el.get("id"), tags={t.get("k"): t.get("v") for t in el.findall("tag")}, nodes=[n.get("ref") for n in el.findall("nd")]))
            el.clear()
    return nodes, ways, pnodes


def read_ms_buildings(paths, bbox, log=print):
    """Microsoft footprints (gzipped GeoJSON lines) inside bbox (west, south, east, north) -> list of lon/lat rings."""
    import gzip
    w, so, e, n = bbox
    rings = []
    for p in paths:
        k = 0
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                g = json.loads(line).get("geometry") or {}
                cs = g.get("coordinates") or []
                ring = cs[0] if g.get("type") == "Polygon" and cs else (cs[0][0] if g.get("type") == "MultiPolygon" and cs else None)
                if not ring:
                    continue
                lon0, lat0 = ring[0][0], ring[0][1]
                if w <= lon0 <= e and so <= lat0 <= n:
                    rings.append(ring)
                    k += 1
        log("  %s: %d buildings in the area" % (os.path.basename(p), k))
    return rings


def _critical_kind(tags):
    for k, vals in CRITICAL.items():
        v = tags.get(k)
        if v and (vals is None or v in vals):
            return "%s=%s" % (k, v)
    return ""


def _area_centroid(x, y):
    a = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    if abs(a) < 1e-9:
        return 0.0, float(np.mean(x)), float(np.mean(y))
    cx = np.sum((x[:-1] + x[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1])) / (6 * a)
    cy = np.sum((y[:-1] + y[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1])) / (6 * a)
    return abs(a), float(cx), float(cy)


class Layer(object):
    """A raster sampled by map coordinates (nearest cell), NaN outside."""

    def __init__(self, arr, gt):
        self.a, self.gt = arr, gt

    def at(self, x, y):
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        c = np.floor((x - self.gt[0]) / self.gt[1]).astype(int)
        r = np.floor((y - self.gt[3]) / self.gt[5]).astype(int)
        ok = (r >= 0) & (r < self.a.shape[0]) & (c >= 0) & (c < self.a.shape[1])
        out = np.full(x.shape, np.nan)
        out[ok] = self.a[r[ok], c[ok]]
        return out


def _ahp_layer(risk_dir, shape, gt0):
    """All complexes' risk maps on the DEM grid (each complex map covers only its own window), highest value where they overlap."""
    fs = sorted(glob.glob(os.path.join(risk_dir, "ahp_C*_risk_abs.tif")))
    if not fs:
        return None, []
    best = np.full(shape, np.nan)
    from .ras import readback as rb
    X, Y = rb.grid_xy(gt0, shape)
    for f in fs:
        a, gt, _ = gio.read_raster(f)
        s = a if (a.shape == tuple(shape) and np.allclose(gt, gt0)) else rb.sample_nearest(a, gt, X, Y)
        best = np.where(np.isfinite(s) & ~(s <= best), s, best)
    return Layer(best, gt0), [os.path.basename(f) for f in fs]


def run(site, osm_path, maps_dir="", run_name="R04", out_dir="", buildings=None, log=print):
    from .ras import readback as rb
    if not os.path.exists(osm_path):
        raise SystemExit("OSM file not found: %s" % osm_path)
    out_dir = out_dir or os.path.join(site.out_dir, "exposure")
    os.makedirs(out_dir, exist_ok=True)
    _dem, gt_dem, crs, _nd = gio.load_dem(site.file("dem"), site.P["REPROJECT_GEOGRAPHIC"], lambda m: None)
    wgs = gio.Crs(epsg=4326, geographic=True)
    log("reading %s ..." % osm_path)
    nodes, ways, pnodes = parse_osm(osm_path)
    log("  %d nodes, %d ways, %d tagged nodes" % (len(nodes), len(ways), len(pnodes)))

    # ---- layers
    breaks = site.P["RISK_CLASS_BREAKS"]
    high = float(breaks[2])
    ahp, ahp_files = _ahp_layer(os.path.join(site.out_dir, "risk"), _dem.shape, gt_dem)
    runs = site.file("runs_dir") or os.path.join(site.out_dir, "hecras", "runs")
    run_dir = os.path.dirname(os.path.abspath(maps_dir)) if maps_dir else os.path.join(runs, run_name)
    hz = dep = None
    hp = os.path.join(run_dir, "readback", "hazard_class.tif")
    if os.path.exists(hp):
        a, g, _ = gio.read_raster(hp)
        hz = Layer(a, g)
    md = maps_dir or os.path.join(run_dir, "maps")
    if os.path.isdir(md):
        p, _k = rb.find_tif(md, "depth", "", False)
        if p:
            a, g, _ = gio.read_raster(p)
            dep = Layer(a, g)
    fp_mask = win_mask = None
    if site.file("flood_footprint"):
        fp_mask = Layer(gio.rasterize([f["geometry"] for f in gio.read_vector(site.file("flood_footprint"), crs)], _dem.shape, gt_dem).astype(float), gt_dem)
        if site.file("survey_window"):
            win_mask = Layer(gio.rasterize([f["geometry"] for f in gio.read_vector(site.file("survey_window"), crs)], _dem.shape, gt_dem).astype(float), gt_dem)
    log("  layers: AHP %s | HEC-RAS hazard %s | depth %s | 2010 footprint %s" % (
        "%d complexes" % len(ahp_files) if ahp else "none", "yes (%s)" % run_dir if hz else "none", "yes" if dep else "none", "yes" if fp_mask else "none"))
    # the area all layers describe: AHP-mapped or HEC-RAS domain
    def classify(x, y):
        r = {}
        v = ahp.at(x, y) if ahp else np.full(np.size(x), np.nan)
        r["ahp_risk"] = v
        r["ahp_high"] = np.nan_to_num(v) >= high
        h = hz.at(x, y) if hz else np.full(np.size(x), np.nan)
        r["hec_class"] = h
        d = dep.at(x, y) if dep else np.full(np.size(x), np.nan)
        r["hec_depth_m"] = d
        f = fp_mask.at(x, y) if fp_mask else np.full(np.size(x), np.nan)
        w = win_mask.at(x, y) if win_mask else np.full(np.size(x), np.nan)
        r["in_2010_footprint"] = np.nan_to_num(f) > 0.5
        r["in_survey_window"] = np.nan_to_num(w) > 0.5
        return r

    # ---- buildings and critical facilities
    B, C = [], []
    for w in ways:
        t = w["tags"]
        pts = [nodes[n] for n in w["nodes"] if n in nodes]
        if len(pts) < 3:
            continue
        crit = _critical_kind(t)
        if "building" not in t and not crit:
            continue
        lon, lat = np.array(pts).T
        x, y = gio.transform_xy(lon, lat, wgs, crs)
        if x[0] != x[-1] or y[0] != y[-1]:
            x, y = np.r_[x, x[0]], np.r_[y, y[0]]
        area, cx, cy = _area_centroid(x, y)
        rec = dict(id=w["id"], kind=crit or "building", name=t.get("name", ""), building=t.get("building", ""), area_m2=round(area, 1), x=cx, y=cy)
        (C if crit else B).append(rec)
        if crit and "building" in t:
            B.append(dict(rec, kind="building"))
    n_osm_buildings = len(B)
    if buildings is None:
        buildings = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(site.file("dem"))), "ms_buildings_*.csv.gz")))
    source = "OpenStreetMap"
    if buildings:
        H, W = _dem.shape
        lon, lat = gio.transform_xy([gt_dem[0], gt_dem[0] + W * gt_dem[1]], [gt_dem[3] + H * gt_dem[5], gt_dem[3]], crs, wgs)
        rings = read_ms_buildings(buildings, (min(lon), min(lat), max(lon), max(lat)), log)
        if rings:
            B = []
            for k, ring in enumerate(rings):
                lo, la = np.array(ring, float)[:, :2].T
                x, y = gio.transform_xy(lo, la, wgs, crs)
                if x[0] != x[-1] or y[0] != y[-1]:
                    x, y = np.r_[x, x[0]], np.r_[y, y[0]]
                area, cx, cy = _area_centroid(x, y)
                B.append(dict(id="ms%d" % k, kind="building", name="", building="", area_m2=round(area, 1), x=cx, y=cy))
            source = "Microsoft Global ML Building Footprints (%d) - OSM has %d" % (len(B), n_osm_buildings)
    log("  buildings from: %s" % source)
    for n in pnodes:
        crit = _critical_kind(n["tags"])
        if crit:
            x, y = gio.transform_xy([n["lon"]], [n["lat"]], wgs, crs)
            C.append(dict(id=n["id"], kind=crit, name=n["tags"].get("name", ""), building="", area_m2=0.0, x=float(x[0]), y=float(y[0])))
    for L in (B, C):
        if L:
            k = classify(np.array([r["x"] for r in L]), np.array([r["y"] for r in L]))
            for i, r in enumerate(L):
                for key, arr in k.items():
                    v = arr[i]
                    r[key] = (bool(v) if arr.dtype == bool else (None if v != v else round(float(v), 3)))

    # ---- roads: length per class, sampled every 10 m
    R = []
    seg_rows = []
    for w in ways:
        t = w["tags"]
        if "highway" not in t:
            continue
        pts = [nodes[n] for n in w["nodes"] if n in nodes]
        if len(pts) < 2:
            continue
        lon, lat = np.array(pts).T
        x, y = gio.transform_xy(lon, lat, wgs, crs)
        s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        if s[-1] <= 0:
            continue
        ss = np.arange(5.0, s[-1], 10.0) if s[-1] > 10 else np.array([s[-1] / 2])
        xs, ys = np.interp(ss, s, x), np.interp(ss, s, y)
        k = classify(xs, ys)
        seg = s[-1] / len(ss)
        hc = np.nan_to_num(k["hec_class"], nan=-1)
        R.append(dict(id=w["id"], highway=t["highway"], name=t.get("name", ""), bridge=t.get("bridge", "") == "yes", length_m=round(float(s[-1]), 1),
                      m_ahp_high=round(float(k["ahp_high"].sum() * seg), 1), m_hec_wet=round(float((hc >= 1).sum() * seg), 1),
                      m_hec_significant_plus=round(float((hc >= 3).sum() * seg), 1), m_hec_extreme=round(float((hc >= 4).sum() * seg), 1),
                      m_2010_footprint=round(float(k["in_2010_footprint"].sum() * seg), 1)))

    # ---- summary
    def count(L, key, test):
        return int(sum(1 for r in L if r.get(key) is not None and test(r[key])))

    def cls_counts(L):
        return {n: count(L, "hec_class", lambda v, k=k: int(round(v)) == k) for k, n in enumerate(HAZ) if k > 0}

    S = dict(osm_file=os.path.abspath(osm_path), building_source=source, run_dir=run_dir if hz else "", ahp_maps=ahp_files, ahp_high_threshold=high,
             buildings=dict(total=len(B), ahp_high=count(B, "ahp_high", bool), hec_wet=count(B, "hec_class", lambda v: v >= 1),
                            hec_by_class=cls_counts(B), in_2010_footprint=count(B, "in_2010_footprint", bool),
                            in_survey_window=count(B, "in_survey_window", bool)),
             critical=dict(total=len(C), ahp_high=count(C, "ahp_high", bool), hec_wet=count(C, "hec_class", lambda v: v >= 1),
                           in_2010_footprint=count(C, "in_2010_footprint", bool)),
             roads_km=dict(total=round(sum(r["length_m"] for r in R) / 1000, 2), ahp_high=round(sum(r["m_ahp_high"] for r in R) / 1000, 2),
                           hec_wet=round(sum(r["m_hec_wet"] for r in R) / 1000, 2), hec_significant_plus=round(sum(r["m_hec_significant_plus"] for r in R) / 1000, 2),
                           hec_extreme=round(sum(r["m_hec_extreme"] for r in R) / 1000, 2), in_2010_footprint=round(sum(r["m_2010_footprint"] for r in R) / 1000, 2)),
             bridges_hit=[dict(name=r["name"] or r["id"], highway=r["highway"]) for r in R if r["bridge"] and (r["m_hec_wet"] > 0 or r["m_ahp_high"] > 0)])
    # building-level check against the 2010 footprint, inside the surveyed window only
    inw = [r for r in B if r.get("in_survey_window")]
    if inw and fp_mask is not None:
        obs = np.array([bool(r["in_2010_footprint"]) for r in inw])
        mod = float(breaks[1])
        for key, test in (("AHP high or above", lambda r: bool(r.get("ahp_high"))),
                          ("AHP moderate or above", lambda r: (r.get("ahp_risk") or 0) >= mod),
                          ("HEC-RAS depth >= 0.3 m", lambda r: (r.get("hec_class") or 0) >= 1),
                          ("HEC-RAS any water", lambda r: (r.get("hec_depth_m") or 0) > 0.01)):
            pred = np.array([test(r) for r in inw])
            tp, fn, fp_ = int((pred & obs).sum()), int((~pred & obs).sum()), int((pred & ~obs).sum())
            S.setdefault("building_check_2010", {})[key] = dict(buildings_in_window=len(inw), flooded_2010=int(obs.sum()), hits=tp, misses=fn, false_alarms=fp_,
                                                                 POD=round(tp / max(1, tp + fn), 3), FAR=round(fp_ / max(1, tp + fp_), 3))

    # ---- files
    with open(os.path.join(out_dir, "exposure_summary.json"), "w") as fh:
        json.dump(S, fh, indent=2)
    cols = ["id", "kind", "name", "building", "area_m2", "x", "y", "ahp_risk", "ahp_high", "hec_class", "hec_depth_m", "in_2010_footprint", "in_survey_window"]
    for fn, L in (("buildings.csv", B), ("critical_facilities.csv", C)):
        with open(os.path.join(out_dir, fn), "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, cols, extrasaction="ignore")
            wr.writeheader()
            wr.writerows(L)
    with open(os.path.join(out_dir, "roads.csv"), "w", newline="", encoding="utf-8") as fh:
        if R:
            wr = csv.DictWriter(fh, list(R[0].keys()))
            wr.writeheader()
            wr.writerows(R)
    feats = [dict(type="Feature", geometry=dict(type="Point", coordinates=[r["x"], r["y"]]), properties={k: r.get(k) for k in cols if k not in ("x", "y")})
             for r in B + [dict(c, kind=c["kind"]) for c in C] if r.get("ahp_high") or (r.get("hec_class") or 0) >= 1 or r.get("in_2010_footprint")]
    doc = dict(type="FeatureCollection", name="exposed", features=feats)
    if crs.epsg:
        doc["crs"] = dict(type="name", properties=dict(name="urn:ogc:def:crs:EPSG::%d" % crs.epsg))
    with open(os.path.join(out_dir, "exposure_buildings.geojson"), "w") as fh:
        json.dump(doc, fh)
    lines = report_lines(S, C)
    with open(os.path.join(out_dir, "exposure_report.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    for l in lines:
        log(l)
    try:
        figure(os.path.join(out_dir, "exposure_map.png"), B, C, ahp, hz, fp_mask, high, site.name, source)
    except Exception as e:                                # the numbers matter more than the picture
        log("  (no figure: %s)" % e)
    log("Written to %s" % out_dir)
    return S


def report_lines(S, C):
    b, c, r = S["buildings"], S["critical"], S["roads_km"]
    L = ["EXPOSURE (today's buildings and roads)",
         "  buildings mapped: %d  [%s]" % (b["total"], S["building_source"]),
         "    in AHP 'high' or above:          %d" % b["ahp_high"],
         "    wet in the HEC-RAS run:          %d   (by hazard class: %s)" % (b["hec_wet"], ", ".join("%s %d" % kv for kv in b["hec_by_class"].items())),
         "    inside the 2010 footprint:       %d   (survey window holds %d)" % (b["in_2010_footprint"], b["in_survey_window"]),
         "  critical facilities: %d; AHP high %d, HEC-RAS wet %d, in the 2010 footprint %d" % (c["total"], c["ahp_high"], c["hec_wet"], c["in_2010_footprint"]),
         "  roads: %.1f km mapped; AHP high %.2f km, HEC-RAS wet %.2f km (significant+ %.2f, extreme %.2f), in the 2010 footprint %.2f km" % (
             r["total"], r["ahp_high"], r["hec_wet"], r["hec_significant_plus"], r["hec_extreme"], r["in_2010_footprint"])]
    if S["bridges_hit"]:
        L.append("  bridges on a flooded / high-risk road: " + ", ".join("%s (%s)" % (x["name"], x["highway"]) for x in S["bridges_hit"]))
    for key, v in S.get("building_check_2010", {}).items():
        L.append("  building check vs 2010, %-23s: %4d of %d flooded buildings flagged (POD %.2f), %4d false alarms (FAR %.2f)" % (
            key, v["hits"], v["flooded_2010"], v["POD"], v["false_alarms"], v["FAR"]))
    hit = [x for x in C if x.get("ahp_high") or (x.get("hec_class") or 0) >= 1 or x.get("in_2010_footprint")]
    if hit:
        L.append("  critical facilities at risk:")
        for x in sorted(hit, key=lambda x: -(x.get("hec_class") or 0)):
            L.append("    %-28s %-26s HEC-RAS %-11s AHP %s%s" % ((x["name"] or x["id"])[:28], x["kind"][:26],
                                                         HAZ[int(x["hec_class"])] if x.get("hec_class") is not None else "-",
                                                         "high" if x.get("ahp_high") else "below high", " | in 2010 footprint" if x.get("in_2010_footprint") else ""))
    L.append("  NOTE: the building layers show what stands today, not in 2010. Roads and facilities: OpenStreetMap (completeness in Ladakh is uneven).")
    return L


def figure(path, B, C, ahp, hz, fp_mask, high, title, source=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    xs = np.array([r["x"] for r in B + C]) if (B or C) else np.array([0.0])
    ys = np.array([r["y"] for r in B + C]) if (B or C) else np.array([0.0])
    focus = [L for L in (hz, fp_mask) if L is not None]
    if focus:                                             # zoom on the modelled reach and the observed flood
        fx, fy = [], []
        for L in focus:
            r, c = np.where(np.nan_to_num(L.a) > 0) if L is fp_mask else np.where(np.isfinite(L.a))
            fx += [L.gt[0] + c.min() * L.gt[1], L.gt[0] + (c.max() + 1) * L.gt[1]]
            fy += [L.gt[3] + (r.max() + 1) * L.gt[5], L.gt[3] + r.min() * L.gt[5]]
        xs, ys = np.array(fx), np.array(fy)
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    for k, (lay, name) in enumerate(((ahp, "AHP risk"), (hz, "HEC-RAS hazard class"))):
        a = ax[k]
        if lay is not None:
            g = lay.gt
            ext = (g[0], g[0] + lay.a.shape[1] * g[1], g[3] + lay.a.shape[0] * g[5], g[3])
            if k == 0:
                a.imshow(np.where(lay.a >= high, 1.0, np.where(np.isfinite(lay.a), 0.0, np.nan)), extent=ext, cmap=ListedColormap(["#eeeeee", "#f4a582"]), vmin=0, vmax=1,
                         interpolation="nearest")
            else:
                a.imshow(np.where(lay.a >= 1, lay.a, np.nan), extent=ext, cmap=ListedColormap(["#d1e5f0", "#92c5de", "#f4a582", "#b2182b"]), vmin=0.5, vmax=4.5,
                         interpolation="nearest")
        if fp_mask is not None:
            g = fp_mask.gt
            ext = (g[0], g[0] + fp_mask.a.shape[1] * g[1], g[3] + fp_mask.a.shape[0] * g[5], g[3])
            a.contour(np.flipud(fp_mask.a), levels=[0.5], colors="k", linewidths=0.8, extent=ext)
        key = "ahp_high" if k == 0 else "hec_class"
        hitb = [r for r in B if (r.get(key) if k == 0 else (r.get(key) or 0) >= 1)]
        a.scatter([r["x"] for r in B], [r["y"] for r in B], s=2, c="#9a9a9a", lw=0)
        a.scatter([r["x"] for r in hitb], [r["y"] for r in hitb], s=5, c="#1a1a1a", lw=0, label="building hit (%d)" % len(hitb))
        hc = [r for r in C if (r.get(key) if k == 0 else (r.get(key) or 0) >= 1)]
        a.scatter([r["x"] for r in hc], [r["y"] for r in hc], s=40, marker="^", c="#2166ac", edgecolors="white", lw=0.6, label="critical facility hit (%d)" % len(hc))
        a.set_title("%s - buildings and facilities hit (black line: 2010 footprint)" % name, loc="left", fontsize=10)
        a.set_xlim(xs.min() - 700, xs.max() + 700)
        a.set_ylim(ys.min() - 700, ys.max() + 700)
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
        a.legend(loc="lower right", fontsize=8, frameon=True)
    fig.suptitle("Exposure - %s. Buildings: %s; facilities: OpenStreetMap" % (title, source.split(" (")[0] or "OpenStreetMap"), x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
