# -*- coding: utf-8 -*-
"""Score maps against an observed flood: a footprint polygon, damage points, and the window that was actually surveyed.

    fanrisk validate <site.toml> [--maps <HEC-RAS Max maps folder>]

Cells inside the footprint are positives, cells inside the survey window but outside the footprint are negatives; nothing outside the
window is scored (nobody looked there). Scores:
    AUC        how well a continuous map ranks flooded above dry cells (0.5 = random, 1 = perfect)
    POD / FAR / CSI   for a yes/no map (HEC-RAS depth >= wet_depth_m; AHP class >= high)
    points     share of damage points the map marks (all points and the held-out third separately)
An AHP map that does not cover part of the window scores 0 there ("not at risk") - leaving those cells out would flatter it.
"""
import json
import os

import numpy as np

from . import gio


def auc(score, pos, neg):
    """Mann-Whitney AUC with average ranks for ties."""
    s1, s0 = score[pos], score[neg]
    n1, n0 = s1.size, s0.size
    if n1 == 0 or n0 == 0:
        return float("nan")
    v = np.concatenate([s1, s0])
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(v.size)
    sv = v[order]
    i = 0
    while i < sv.size:                                    # average ranks of ties
        j = i
        while j + 1 < sv.size and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def categorical(yes, pos, neg):
    tp, fn, fp = int((yes & pos).sum()), int((~yes & pos).sum()), int((yes & neg).sum())
    return dict(POD=tp / max(1, tp + fn), FAR=fp / max(1, tp + fp), CSI=tp / max(1, tp + fn + fp))


def _ring_mask(features, shape, gt):
    return gio.rasterize([f["geometry"] for f in features], shape, gt)


def score(grid_shape, gt, crs, footprint, window, maps, points=None, domain=None, high=0.55, wet_depth=0.3):
    """maps: {name: (array on the grid, kind)} with kind 'depth' or 'risk'. footprint / window / points: file paths. domain: optional bool array
    (e.g. the HEC-RAS model area) that further limits the scored cells. Returns a dict of results."""
    fp = _ring_mask(gio.read_vector(footprint, crs), grid_shape, gt)
    win = _ring_mask(gio.read_vector(window, crs), grid_shape, gt) if window else np.ones(grid_shape, bool)
    ev = win & (domain if domain is not None else True)
    pos, neg = fp & ev, ~fp & ev
    out = dict(cells_scored=int(ev.sum()), footprint_cells=int(pos.sum()), cell_m=abs(gt[1]), maps={})
    pts = []
    if points:
        for f in gio.read_vector(points, crs):
            g = f["geometry"]
            if g["type"] == "Point":
                pts.append((g["coordinates"][0], g["coordinates"][1], int((f["properties"] or {}).get("holdout") or 0)))
    P = np.array(pts) if pts else np.zeros((0, 3))
    if P.size:
        c = np.floor((P[:, 0] - gt[0]) / gt[1]).astype(int)
        r = np.floor((P[:, 1] - gt[3]) / gt[5]).astype(int)
        inside = (r >= 0) & (r < grid_shape[0]) & (c >= 0) & (c < grid_shape[1])
    for name, (arr, kind) in maps.items():
        a = np.asarray(arr, dtype=float)
        covered = np.isfinite(a)
        if kind == "risk" and not (covered & ev).any():          # a complex elsewhere on the DEM: nothing to score
            continue
        a0 = np.where(covered, a, 0.0)
        res = dict(kind=kind, AUC=auc(a0, pos, neg), footprint_covered_pct=100.0 * float((covered & pos).sum()) / max(1, int(pos.sum())))
        yes = a0 >= (wet_depth if kind == "depth" else high)
        res.update(categorical(yes, pos, neg))
        res["threshold"] = wet_depth if kind == "depth" else high
        if P.size:
            v = np.zeros(len(P))
            v[inside] = a0[r[inside], c[inside]]
            hit = v >= res["threshold"]
            res["points_hit_pct"] = 100.0 * float(hit.mean())
            if (P[:, 2] == 1).any():
                res["holdout_points_hit_pct"] = 100.0 * float(hit[P[:, 2] == 1].mean())
            # AUC of the points against the scored background (cells outside the footprint)
            bg = a0[neg]
            if bg.size:
                vv = np.concatenate([v, bg])
                lab_ = np.zeros(vv.size, bool)
                lab_[:v.size] = True
                res["points_AUC"] = auc(vv, lab_, ~lab_)
        out["maps"][name] = res
    return out


def report(res, log=print):
    log("Validation against the observed flood: %d cells of %.0f m scored, %d inside the footprint" % (
        res["cells_scored"], res["cell_m"], res["footprint_cells"]))
    log("  %-28s %6s %6s %6s %6s %9s %10s" % ("map", "AUC", "POD", "FAR", "CSI", "covers", "points hit"))
    for n, m in res["maps"].items():
        log("  %-28s %6.2f %6.2f %6.2f %6.2f %8.0f%% %9s" % (n[:28], m["AUC"], m["POD"], m["FAR"], m["CSI"], m["footprint_covered_pct"],
                                                         ("%.0f%%" % m["points_hit_pct"]) if "points_hit_pct" in m else "-"))


def validate_site(site, maps_dir="", out_dir="", log=print):
    """Score the site's AHP risk maps (every complex) and, if given, a folder of HEC-RAS Max maps against the validation case."""
    from .ras import readback as rb
    fpath, wpath, ppath = site.file("flood_footprint"), site.file("survey_window"), site.file("flood_points")
    if not fpath:
        raise SystemExit("No observed footprint: set flood_footprint (and survey_window) in [validation] of site.toml or in the validation case.")
    risk_dir = os.path.join(site.out_dir, "risk")
    if not os.path.isdir(risk_dir):
        raise SystemExit("No risk maps in %s: run 'fanrisk run' first." % risk_dir)
    # the scoring grid = the DEM grid
    dem, gt, crs, _nd = gio.load_dem(site.file("dem"), site.P["REPROJECT_GEOGRAPHIC"], lambda m: None)
    shape = dem.shape
    X, Y = rb.grid_xy(gt, shape)
    maps = {}
    for fn in sorted(os.listdir(risk_dir)):
        if fn.endswith("_risk_abs.tif"):
            a, agt, _ = gio.read_raster(os.path.join(risk_dir, fn))
            maps["AHP " + fn[4:-13] + " risk"] = (rb.sample_nearest(a, agt, X, Y), "risk")
    if len(maps) > 1:                                     # the complexes together, as one map
        stack = np.stack([m[0] for m in maps.values()])
        maps["AHP all complexes risk"] = (np.where(np.isfinite(stack).any(0), np.nanmax(np.where(np.isfinite(stack), stack, -1), 0), np.nan), "risk")
    domain = None
    if maps_dir:
        (dp, _k), _v, _a = rb.find_tif(maps_dir, "depth", "", False), None, None
        if dp is None:
            raise SystemExit("No Max depth map in %s" % maps_dir)
        d, dgt, _ = gio.read_raster(dp)
        maps["HEC-RAS depth (%s)" % os.path.basename(os.path.dirname(os.path.abspath(maps_dir)))] = (rb.sample_nearest(np.nan_to_num(d), dgt, X, Y), "depth")
        hd = site.handoff_dir()
        if hd:
            domain = rb.points_in_polygon(X, Y, rb.load_ring(hd))
    res = score(shape, gt, crs, fpath, wpath, maps, ppath, domain, high=float(site.P["RISK_CLASS_BREAKS"][2]), wet_depth=float(site.score["wet_depth_m"]))
    res["scored_within"] = "survey window" + (" and the HEC-RAS model area" if domain is not None else "")
    report(res, log)
    out_dir = out_dir or os.path.join(site.out_dir, "validation")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "validation_summary.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    log("Written to %s" % os.path.join(out_dir, "validation_summary.json"))
    return res
