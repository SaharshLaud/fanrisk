# -*- coding: utf-8 -*-
"""Check a run against a documented event when there is no flood footprint, only a place and a few published numbers.

    python tools/check_event.py <site.toml> --lon 78.7972 --lat 31.0407 [--catchment-km2 17] [--melton 0.8] [--deposit-ha 20]

Reports: did the tool find a fan at the place; which complex / feeder drains onto it; its catchment area and Melton ratio against
the published values; the flow type the tool assigns (Melton breaks); the AHP risk class at the place and around it; how much of
that fan the AHP marks high or above, against the published deposit area; and where the complex ranks among all complexes.
Writes <output_dir>/event_check/event_check.json and event_check.txt.
"""
import argparse
import csv
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from fanrisk import config, gio  # noqa: E402

CLASSES = ("very low", "low", "moderate", "high", "very high")


def ring_contains(ring, x, y):
    inside = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i][:2]
        x2, y2 = ring[(i + 1) % n][:2]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside


def polys(g):
    if g["type"] == "Polygon":
        return [g["coordinates"]]
    if g["type"] == "MultiPolygon":
        return g["coordinates"]
    return []


def dist_to_poly(g, x, y):
    best = float("inf")
    for p in polys(g):
        if ring_contains(p[0], x, y) and not any(ring_contains(h, x, y) for h in p[1:]):
            return 0.0
        for ring in p:
            a = np.asarray(ring, float)[:, :2]
            for (x1, y1), (x2, y2) in zip(a[:-1], a[1:]):
                dx, dy = x2 - x1, y2 - y1
                t = 0.0 if dx == dy == 0 else max(0.0, min(1.0, ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)))
                best = min(best, math.hypot(x - x1 - t * dx, y - y1 - t * dy))
    return best


def read_layer(site, name):
    """A vector layer of the run: from fans/fanrisk.gpkg (geopandas installed) or fans/vectors/<name>.geojson. [] if absent."""
    gpkg = os.path.join(site.out_dir, "fans", "fanrisk.gpkg")
    gj = os.path.join(site.out_dir, "fans", "vectors", name + ".geojson")
    if os.path.exists(gpkg):
        try:
            import geopandas as gpd
            from shapely.geometry import mapping
            gdf = gpd.read_file(gpkg, layer=name)
            cols = [c for c in gdf.columns if c != gdf.geometry.name]
            return [dict(geometry=json.loads(json.dumps(mapping(r[gdf.geometry.name]))), properties={c: r[c] for c in cols})
                    for _, r in gdf.iterrows() if r[gdf.geometry.name] is not None]
        except Exception:
            pass
    return gio.read_vector(gj) if os.path.exists(gj) else []


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("site")
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--catchment-km2", type=float)
    ap.add_argument("--melton", type=float)
    ap.add_argument("--deposit-ha", type=float)
    ap.add_argument("--radius-m", type=float, default=250.0, help="neighbourhood around the place for the risk summary")
    a = ap.parse_args(argv)
    site = config.load_site(a.site)
    out_dir = os.path.join(site.out_dir, "event_check")
    os.makedirs(out_dir, exist_ok=True)
    fans = read_layer(site, "fans")
    feeders = read_layer(site, "feeders")
    _d, gt_dem, crs, _n = gio.load_dem(site.file("dem"), site.P["REPROJECT_GEOGRAPHIC"], lambda m: None)
    x, y = gio.transform_xy([a.lon], [a.lat], gio.Crs(epsg=4326, geographic=True), crs)
    x, y = float(x[0]), float(y[0])
    L = []
    say = L.append
    R = dict(place=dict(lon=a.lon, lat=a.lat, x=x, y=y), published=dict(catchment_km2=a.catchment_km2, melton=a.melton, deposit_ha=a.deposit_ha))
    say("EVENT CHECK - %s" % site.name)
    say("  place: %.4f N, %.4f E  (map x %.0f, y %.0f)" % (a.lat, a.lon, x, y))
    say("  fans found on the DEM: %d" % len(fans))
    if not fans:
        say("  NO FAN FOUND - the check fails at the first step")
        R["fan_found"] = False
    else:
        d = [(dist_to_poly(f["geometry"], x, y), f) for f in fans]
        dist, fan = min(d, key=lambda t: t[0])
        p = fan["properties"]
        R.update(fan_found=dist == 0.0, distance_to_fan_m=dist, fan=p)
        say("  nearest fan: fan %s in complex C%s, %s (%.0f m away)" % (p.get("fan_id"), p.get("complex_id"), "the place is ON it" if dist == 0 else "the place is NOT on it", dist))
        say("    fan area %.2f km2, mean slope %.1f deg, apex %.0f m" % (p.get("area_km2") or 0, p.get("mean_slope") or 0, p.get("apex_elev") or 0))
        fd = [f["properties"] for f in feeders if f["properties"].get("fan_id") == p.get("fan_id")]
        fd = max(fd, key=lambda q: q.get("catch_km2") or 0) if fd else None
        cat = (fd or p).get("catch_km2") or p.get("catch_km2")
        mel = (fd or p).get("melton") or p.get("melton")
        lo, hi = site.P["RAS_MELTON_BREAKS"]
        ftype = "clear-water flood" if mel < lo else ("debris flood" if mel < hi else "debris flow")
        R.update(feeder=fd, catchment_km2=cat, melton=mel, flow_type=ftype)
        say("  catchment draining onto it: %.1f km2%s" % (cat, "   (published %.1f km2: %+.0f%%)" % (a.catchment_km2, 100 * (cat - a.catchment_km2) / a.catchment_km2) if a.catchment_km2 else ""))
        say("  Melton ratio: %.2f%s" % (mel, "   (published %.2f)" % a.melton if a.melton else ""))
        say("  flow type from the Melton breaks %.1f / %.1f: %s" % (lo, hi, ftype.upper()))
        cf = [f["properties"] for f in feeders if f["properties"].get("complex_id") == p.get("complex_id")]
        if len(cf) > 1:
            big = max(cf, key=lambda q: q.get("catch_km2") or 0)
            R["complex_largest_feeder"] = big
            say("  (complex C%s has %d feeders; the largest, %s, drains %.1f km2 with Melton %.2f)" % (
                p.get("complex_id"), len(cf), big.get("name"), big.get("catch_km2") or 0, big.get("melton") or 0))
    # AHP risk at and around the place, and on the fan
    fs = sorted(glob.glob(os.path.join(site.out_dir, "risk", "ahp_C*_risk_abs.tif")))
    br = list(site.P["RISK_CLASS_BREAKS"])
    best_val, best_name, fan_high_ha = None, None, None
    for f in fs:
        arr, g, _ = gio.read_raster(f)
        c = int(math.floor((x - g[0]) / g[1]))
        r = int(math.floor((y - g[3]) / g[5]))
        rad = int(math.ceil(a.radius_m / abs(g[1])))
        if not (0 <= r < arr.shape[0] and 0 <= c < arr.shape[1]):
            continue
        v = arr[r, c]
        win = arr[max(0, r - rad):r + rad + 1, max(0, c - rad):c + rad + 1]
        if np.isfinite(v) and (best_val is None or v > best_val):
            best_val, best_name = float(v), os.path.basename(f)
            vals = win[np.isfinite(win)]
            R["risk_at_place"] = dict(map=best_name, value=best_val, cls=CLASSES[int(np.digitize(best_val, br))],
                                      neighbourhood_max=float(vals.max()) if vals.size else None,
                                      neighbourhood_share_high=float((vals >= br[2]).mean()) if vals.size else None)
            if R.get("fan"):
                H, W = arr.shape
                X = g[0] + (np.arange(W) + 0.5) * g[1]
                Y = g[3] + (np.arange(H) + 0.5) * g[5]
                m = gio.rasterize([fan["geometry"]], arr.shape, g)
                hi_cells = m & (np.nan_to_num(arr, nan=-1) >= br[2])
                fan_high_ha = float(hi_cells.sum() * abs(g[1] * g[5]) / 1e4)
                R["fan_high_ha"] = fan_high_ha
    if best_val is None:
        say("  AHP risk at the place: NOT MAPPED (no complex map covers it)")
    else:
        q = R["risk_at_place"]
        say("  AHP risk at the place: %.2f = %s  (%s)" % (q["value"], q["cls"].upper(), q["map"]))
        say("    within %.0f m: highest %.2f, %.0f%% of cells high or above" % (a.radius_m, q["neighbourhood_max"], 100 * q["neighbourhood_share_high"]))
        if fan_high_ha is not None:
            say("  area of that fan marked high or above: %.1f ha%s" % (fan_high_ha, "   (published deposit about %.0f ha)" % a.deposit_ha if a.deposit_ha else ""))
    # rank of the complex
    mt = os.path.join(site.out_dir, "risk", "master_table.csv")
    if os.path.exists(mt) and R.get("fan"):
        rows = list(csv.DictReader(open(mt)))
        cid = str(R["fan"].get("complex_id"))
        sc = sorted({(r["complex_id"], max(float(q["C_source"]) for q in rows if q["complex_id"] == r["complex_id"])) for r in rows}, key=lambda t: -t[1])
        rank = [i for i, (k, _v) in enumerate(sc, 1) if k == cid]
        if rank:
            R["complex_rank"] = dict(rank=rank[0], of=len(sc), source_score=dict(sc)[cid])
            say("  complex C%s ranks %d of %d complexes by catchment source score (C = %.2f)" % (cid, rank[0], len(sc), dict(sc)[cid]))
    with open(os.path.join(out_dir, "event_check.json"), "w") as fh:
        json.dump(R, fh, indent=2, default=float)
    with open(os.path.join(out_dir, "event_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("Written to %s" % out_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
