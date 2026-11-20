# -*- coding: utf-8 -*-
"""Score HEC-RAS 2D maximum maps: hazard rating, comparison with observations (footprint, points, widths, depth / speed at a section) and
with the AHP risk map. The numerical part is ras_readback.py unchanged; files are read through fanrisk.gio (no GDAL needed) and every
observation comes from the validation case of the site - nothing about a particular flood is built in.
"""
import json
import math
import os
import re

import numpy as np

from .. import gio


# =============================================================================
# pure functions (numpy only)
# =============================================================================
def grid_xy(gt, shape):
    rows, cols = np.arange(shape[0]), np.arange(shape[1])
    X = np.broadcast_to(gt[0] + (cols + 0.5) * gt[1], shape)
    Y = np.broadcast_to((gt[3] + (rows + 0.5) * gt[5])[:, None], shape)
    return X, Y


def points_in_polygon(X, Y, ring):
    """Even-odd rule, vectorised over the grid."""
    inside = np.zeros(X.shape, dtype=bool)
    ring = np.asarray(ring, dtype=float)
    for (x1, y1), (x2, y2) in zip(ring[:-1], ring[1:]):
        if y1 == y2:
            continue
        cond = (y1 > Y) != (y2 > Y)
        xint = (x2 - x1) * (Y - y1) / (y2 - y1) + x1
        inside ^= cond & (X < xint)
    return inside


def polyline_frame(px, py, X, Y):
    """Chainage s, signed offset d (+ left of travel) and `inside` (nearest point not beyond the ends) of every point."""
    n = len(px)
    seg = np.hypot(np.diff(px), np.diff(py))
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    best = np.full(X.shape, np.inf)
    s_out, d_out, inside = np.zeros(X.shape), np.zeros(X.shape), np.ones(X.shape, dtype=bool)
    for i in range(n - 1):
        ax, ay, dx, dy = px[i], py[i], px[i + 1] - px[i], py[i + 1] - py[i]
        if dx == 0 and dy == 0:
            continue                                    # a repeated vertex has no direction
        t = ((X - ax) * dx + (Y - ay) * dy) / (dx * dx + dy * dy)
        tc = np.clip(t, 0.0, 1.0)
        dist = np.hypot(X - (ax + tc * dx), Y - (ay + tc * dy))
        better = dist < best
        side = np.where(dx * (Y - ay) - dy * (X - ax) >= 0, 1.0, -1.0)
        best = np.where(better, dist, best)
        s_out = np.where(better, cum[i] + tc * seg[i], s_out)
        d_out = np.where(better, side * dist, d_out)
        beyond = ((i == 0) & (t < 0)) | ((i == n - 2) & (t > 1))
        inside = np.where(better, ~beyond, inside)
    return s_out, d_out, inside


def dist_to_polyline(X, Y, line):
    """Distance of every grid point to a polyline [(x, y), ...]."""
    out = np.full(X.shape, np.inf)
    L = np.asarray(line, dtype=float)
    for (x1, y1), (x2, y2) in zip(L[:-1], L[1:]):
        vx, vy = x2 - x1, y2 - y1
        t = np.clip(((X - x1) * vx + (Y - y1) * vy) / max(vx * vx + vy * vy, 1e-12), 0.0, 1.0)
        out = np.minimum(out, np.hypot(X - (x1 + t * vx), Y - (y1 + t * vy)))
    return out


def hazard_rating(depth, vel, df):
    return np.where(depth > 0, depth * (np.nan_to_num(vel) + 0.5) + df, 0.0)


def hazard_class(hr, wet, breaks):
    """0 = dry, 1 = low, 2 = moderate, 3 = significant, 4 = extreme."""
    return np.where(wet, 1 + np.digitize(hr, breaks), 0).astype(np.uint8)


def categorical_scores(mod, obs, domain):
    m, o = mod & domain, obs & domain
    hit, miss, fa = int((m & o).sum()), int((~m & o).sum()), int((m & ~o).sum())
    div = lambda a, b: float(a) / b if b else float("nan")
    return dict(hits=hit, misses=miss, false_alarms=fa, POD=div(hit, hit + miss), FAR=div(fa, hit + fa), CSI=div(hit, hit + miss + fa),
                bias=div(hit + fa, hit + miss), observed_cells=int(o.sum()), modelled_cells=int(m.sum()))


def wet_width_at(s, d, wet, s0, slab, res):
    """Modelled flood width across the centreline at chainage s0: (contiguous width around the channel, total wet width), metres."""
    sel = (np.abs(s - s0) <= slab / 2.0) & np.isfinite(d)
    if not sel.any():
        return float("nan"), float("nan")
    dd, ww = d[sel], wet[sel]
    nslab = max(1.0, float(slab / res))
    total = float(ww.sum()) * res / nslab                              # wet area / slab thickness
    nb = int(math.ceil(2.0 * np.abs(dd).max() / res)) + 3
    edges = np.linspace(-nb / 2.0 * res, nb / 2.0 * res, nb + 1)
    occ = np.histogram(dd[ww], bins=edges)[0] > 0
    if not occ.any():
        return 0.0, total
    centre = (nb // 2)
    idx = np.nonzero(occ)[0]
    start = idx[np.argmin(np.abs(idx - centre))]                        # the wet bin nearest the centreline
    lo = hi = start
    gap = 0
    j = start - 1
    while j >= 0 and gap <= 1:                                           # tolerate a one-cell gap (a dry bar)
        if occ[j]:
            lo, gap = j, 0
        else:
            gap += 1
        j -= 1
    gap = 0
    j = start + 1
    while j < nb and gap <= 1:
        if occ[j]:
            hi, gap = j, 0
        else:
            gap += 1
        j += 1
    return float((hi - lo + 1) * res), total


def sample_nearest(arr, gt, x, y):
    c = np.floor((x - gt[0]) / gt[1]).astype(int)
    r = np.floor((y - gt[3]) / gt[5]).astype(int)
    ok = (r >= 0) & (r < arr.shape[0]) & (c >= 0) & (c < arr.shape[1])
    out = np.full(np.shape(x), np.nan)
    out[ok] = arr[r[ok], c[ok]]
    return out


def rank_tied(a):
    order = np.argsort(a, kind="mergesort")
    sa = a[order]
    _, inv, cnt = np.unique(sa, return_inverse=True, return_counts=True)
    cum = np.cumsum(cnt)
    avg = (cum - cnt + 1 + cum) / 2.0
    r = np.empty(a.size)
    r[order] = avg[inv]
    return r


def spearman(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    ra, rb = rank_tied(a[ok]), rank_tied(b[ok])
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def roc_auc(score, presence, domain):
    ok = domain & np.isfinite(score)
    s, p = score[ok], presence[ok]
    npos, nneg = int(p.sum()), int((~p).sum())
    if npos == 0 or nneg == 0:
        return float("nan")
    r = rank_tied(s)
    return float((r[p].sum() - npos * (npos + 1) / 2.0) / (npos * nneg))


def probe_value(arr, gt, x, y, radius):
    X, Y = grid_xy(gt, arr.shape)
    near = np.hypot(X - x, Y - y) <= max(radius, abs(gt[1]))
    v = arr[near]
    v = v[np.isfinite(v)]
    return float(v.max()) if v.size else float("nan")


def centerline_point(cl, chainage, offset=0.0):
    x = float(np.interp(chainage, cl["s"], cl["x"]))
    y = float(np.interp(chainage, cl["s"], cl["y"]))
    k = int(np.clip(np.searchsorted(cl["s"], chainage), 1, len(cl["s"]) - 1))
    tx, ty = cl["x"][k] - cl["x"][k - 1], cl["y"][k] - cl["y"][k - 1]
    n = math.hypot(tx, ty) or 1.0
    return x - ty / n * offset, y + tx / n * offset


def evaluate(depth, vel, arrival, gt, cl, domain, obs_mask, obs_points, ahp, P):
    """All scores. depth / vel / arrival: arrays (NaN = no water), gt: their geotransform, cl: centreline dict, domain: bool array of the
    model area, obs_mask: observed footprint (bool array or None), obs_points: [(x, y)] or None, ahp: (array, gt) or None. P: settings dict."""
    res = abs(gt[1])
    d0 = np.nan_to_num(depth, nan=0.0)
    wet = (d0 >= P["WET_DEPTH_M"]) & domain
    hr = hazard_rating(d0, vel, P["DEBRIS_FACTOR"])
    hc = hazard_class(hr, wet, P["HAZARD_BREAKS"])
    cell_km2 = res * abs(gt[5]) / 1e6
    out = dict(model=dict(cell_m=float(res), wet_area_km2=float(wet.sum() * cell_km2), max_depth_m=float(np.nanmax(d0[domain])) if domain.any() else float("nan"),
                          max_speed_ms=float(np.nanmax(np.nan_to_num(vel)[domain])) if domain.any() else float("nan"),
                          hazard_class_area_km2={n: float(((hc == k) & domain).sum() * cell_km2) for k, n in enumerate(("dry", "low", "moderate", "significant", "extreme"))}))
    X, Y = grid_xy(gt, depth.shape)
    s, d, inside = polyline_frame(np.asarray(cl["x"]), np.asarray(cl["y"]), X, Y)
    inner = domain.copy()
    for _ in range(2):                                                     # does water reach the edge of the model area?
        p_ = np.pad(inner, 1)
        inner = inner & p_[:-2, 1:-1] & p_[2:, 1:-1] & p_[1:-1, :-2] & p_[1:-1, 2:]
    edge_wet = wet & domain & ~inner
    at_bc = np.zeros(depth.shape, dtype=bool)                              # water standing AT the inflow / outflow lines is expected
    for line in (P.get("BC_LINES") or []):
        at_bc |= dist_to_polyline(X, Y, line) <= 3.0 * abs(gt[1])
    out["model"]["wet_cells_at_boundary_lines"] = int((edge_wet & at_bc).sum())
    edge_wet = edge_wet & ~at_bc
    out["model"]["boundary_lines_known"] = bool(P.get("BC_LINES"))
    out["model"]["wet_cells_at_model_edge"] = int(edge_wet.sum())
    out["model"]["wet_share_at_model_edge_pct"] = 100.0 * float(edge_wet.sum()) / max(1, int(wet.sum()))
    if edge_wet.any():
        out["model"]["edge_chainage_km"] = [float(s[edge_wet].min() / 1000.0), float(s[edge_wet].max() / 1000.0)]
        out["model"]["edge_max_depth_m"] = float(np.nanmax(d0[edge_wet]))
    if obs_mask is not None:
        sc = categorical_scores(wet, obs_mask, domain)
        sc["observed_inside_model_area_pct"] = 100.0 * float((obs_mask & domain).sum()) / max(1, int(obs_mask.sum()))
        out["footprint"] = sc
    if obs_points:
        px, py = np.array([p[0] for p in obs_points]), np.array([p[1] for p in obs_points])
        near = np.array([probe_value(d0, gt, a, b, P["POINT_BUFFER_M"]) for a, b in zip(px, py)])
        inside_pts = sample_nearest(domain.astype(float), gt, px, py) > 0
        used = inside_pts
        out["points"] = dict(n_points=int(len(px)), n_inside_model_area=int(used.sum()),
                             hit_rate=float(np.mean(near[used] >= P["WET_DEPTH_M"])) if used.any() else float("nan"),
                             median_depth_at_points_m=float(np.nanmedian(near[used])) if used.any() else float("nan"))
    widths = []
    for ch, w_obs in P["OBSERVED_WIDTHS"]:
        wc, wt = wet_width_at(s, d, wet & inside, ch, max(2.0 * res, 20.0), res)
        widths.append(dict(chainage_m=float(ch), observed_m=float(w_obs), modelled_contiguous_m=wc, modelled_total_m=wt,
                           error_pct=100.0 * (wc - w_obs) / w_obs if wc == wc else float("nan"),
                           within_tolerance=bool(abs(wc - w_obs) <= P["TOLERANCE"] * w_obs) if wc == wc else False))
    out["widths"] = widths
    probes = []
    for pr in P["PROBES"]:
        x, y = centerline_point(cl, pr["chainage_m"], pr.get("offset_m", 0.0))
        md = probe_value(d0, gt, x, y, P["PROBE_RADIUS_M"])
        mv = probe_value(np.nan_to_num(vel), gt, x, y, P["PROBE_RADIUS_M"])
        probes.append(dict(name=pr["name"], x=x, y=y, depth_obs_m=pr.get("depth_m"), depth_model_m=md,
                           depth_error_pct=100.0 * (md - pr["depth_m"]) / pr["depth_m"] if pr.get("depth_m") else None,
                           depth_ok=bool(pr.get("depth_m") and abs(md - pr["depth_m"]) <= P["TOLERANCE"] * pr["depth_m"]),
                           speed_obs_ms=pr.get("velocity_ms"), speed_model_ms=mv,
                           speed_error_pct=100.0 * (mv - pr["velocity_ms"]) / pr["velocity_ms"] if pr.get("velocity_ms") else None,
                           speed_ok=bool(pr.get("velocity_ms") and abs(mv - pr["velocity_ms"]) <= P["TOLERANCE"] * pr["velocity_ms"])))
    out["probes"] = probes
    if arrival is not None:
        a_ok = np.isfinite(arrival) & domain
        if a_ok.any():
            ends = a_ok & (s > np.nanmax(s[domain]) - 3 * res)
            out["arrival"] = dict(first_water_min=float(np.nanmin(arrival[a_ok])), at_outflow_end_min=float(np.nanmin(arrival[ends])) if ends.any() else float("nan"))
    if ahp is not None:
        ahp_arr, ahp_gt = ahp
        risk = sample_nearest(ahp_arr, ahp_gt, X, Y)
        ev = domain & np.isfinite(risk)
        bin_ = np.digitize(risk[ev], P.get("RISK_BREAKS", (0.2, 0.4, 0.55, 0.7)))
        wv = wet[ev]
        ew = ev & wet
        near = ev & inside & (np.abs(d) <= 300.0)
        out["ahp"] = dict(cells_compared=int(ev.sum()), spearman_risk_vs_hazard_rating=spearman(risk[ev], hr[ev]),
                          spearman_risk_vs_depth=spearman(risk[ev], d0[ev]), roc_auc_risk_vs_wet=roc_auc(risk, wet, ev),
                          spearman_risk_vs_depth_wet_cells_only=spearman(risk[ew], d0[ew]),
                          roc_auc_risk_vs_deep_ge_1m=roc_auc(risk, (d0 >= 1.0) & domain, ev),
                          roc_auc_risk_vs_wet_within_300m_of_centreline=roc_auc(risk, wet, near),
                          wet_cells_share_by_ahp_class_pct=[100.0 * float(wv[bin_ == k].sum()) / max(1, int(wv.sum())) for k in range(5)],
                          wet_share_by_ahp_class_pct=[100.0 * float(wv[bin_ == k].mean()) if (bin_ == k).any() else float("nan") for k in range(5)],
                          area_share_by_ahp_class_pct=[100.0 * float((bin_ == k).mean()) for k in range(5)])
    return out, hc, hr


def file_time(name):
    """(h, m, s) if the file name carries a time step such as 00_17_00 or 00:17:00, else None."""
    import re
    m = re.search(r"(?<!\d)(\d{2})[_ :](\d{2})[_ :](\d{2})(?!\d)", name)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None


def find_tif(folder, key, explicit, allow_snapshot=False):
    """Pick the Max map. Returns (path, kind) with kind 'max', 'snapshot' or 'unlabelled'; stops with a clear message otherwise."""
    if explicit:
        p = explicit if os.path.isabs(explicit) else os.path.join(folder, explicit)
        return p, ("max" if "max" in os.path.basename(p).lower() else ("snapshot" if file_time(os.path.basename(p)) else "unlabelled"))
    c = [f for f in sorted(os.listdir(folder)) if f.lower().endswith((".tif", ".tiff", ".vrt")) and key in f.lower()]
    if not c:
        return None, None
    c.sort(key=lambda f: (not f.lower().endswith(".vrt"), f))          # a .vrt covers all of RAS Mapper's tiles: prefer it
    mx = [f for f in c if "max" in f.lower()]
    if mx:
        return os.path.join(folder, mx[0]), "max"
    if not allow_snapshot:
        raise SystemExit("No MAXIMUM '%s' map in %s. Found only: %s.\nExport 'Depth (Max)' / 'Velocity (Max)' (profile Max), or set [score] allow_snapshot = true "
                         "to look at one time step on purpose." % (key, folder, ", ".join(c)))
    f = c[0]
    return os.path.join(folder, f), ("snapshot" if file_time(f) else "unlabelled")


def load_ring(handoff_dir):
    with open(os.path.join(handoff_dir, "flow_area_polygon.geojson")) as fh:
        return json.load(fh)["features"][0]["geometry"]["coordinates"][0]


def _fmt(v, f="%.2f"):
    return (f % v) if (v is not None and v == v) else "n/a"


def report(o, log=print, wet_depth_m=0.3):
    m = o["model"]
    log("HEC-RAS result: flooded area %.3f km2 (depth >= %.2f m), deepest %.1f m, fastest %.1f m/s" % (m["wet_area_km2"], wet_depth_m, m["max_depth_m"], m["max_speed_ms"]))
    log("  hazard rating classes (km2): " + ", ".join("%s %.3f" % kv for kv in m["hazard_class_area_km2"].items() if kv[0] != "dry"))
    if m.get("wet_cells_at_model_edge"):
        where = "away from the inflow / outflow lines" if m.get("boundary_lines_known") else "(inflow / outflow lines not known, so cells at them are counted too)"
        log("  WARNING: %d flooded cells (%.0f%%), up to %.1f m deep, touch the edge of the model area %s, at chainage %.1f-%.1f km: water is held back "
            "by the model boundary there. Widen the area or extend the outflow line ('fanrisk hecras extend-outflow') and re-run." % (
                m["wet_cells_at_model_edge"], m["wet_share_at_model_edge_pct"], m.get("edge_max_depth_m", float("nan")), where, m["edge_chainage_km"][0], m["edge_chainage_km"][1]))
    elif m.get("boundary_lines_known"):
        log("  model edge: no flooded cell touches the boundary away from the inflow / outflow lines (%d at the lines, as expected)" % m.get("wet_cells_at_boundary_lines", 0))
    if "footprint" in o:
        f = o["footprint"]
        log("vs observed footprint: hit rate (POD) %s, false-alarm ratio %s, critical success index %s, bias %s | %s%% of the observed area lies inside the model area" % (
            _fmt(f["POD"]), _fmt(f["FAR"]), _fmt(f["CSI"]), _fmt(f["bias"]), _fmt(f["observed_inside_model_area_pct"], "%.0f")))
    if "points" in o:
        p = o["points"]
        log("vs observed points: %d of %d inside the model area; %s%% reached by water (median depth there %s m)" % (
            p["n_inside_model_area"], p["n_points"], _fmt(100 * p["hit_rate"], "%.0f"), _fmt(p["median_depth_at_points_m"])))
    for w in o["widths"]:
        log("width at %.0f m: observed %.0f m, modelled %s m (%s%%) %s" % (w["chainage_m"], w["observed_m"], _fmt(w["modelled_contiguous_m"], "%.0f"),
                                                                           _fmt(w["error_pct"], "%+.0f"), "OK" if w["within_tolerance"] else "outside tolerance"))
    for p in o["probes"]:
        log("%s: depth observed %s m, modelled %s m (%s%%) %s | speed observed %s, modelled %s m/s (%s%%) %s" % (
            p["name"], _fmt(p["depth_obs_m"]), _fmt(p["depth_model_m"]), _fmt(p["depth_error_pct"], "%+.0f"), "OK" if p["depth_ok"] else "outside tolerance",
            _fmt(p["speed_obs_ms"]), _fmt(p["speed_model_ms"]), _fmt(p["speed_error_pct"], "%+.0f"), "OK" if p["speed_ok"] else "outside tolerance"))
    if "arrival" in o:
        log("arrival: first water at %.1f min, reaches the outflow end at %s min" % (o["arrival"]["first_water_min"], _fmt(o["arrival"]["at_outflow_end_min"], "%.1f")))
    if "ahp" in o:
        a = o["ahp"]
        log("AHP vs HEC-RAS (%d cells): rank correlation with hazard rating %s, with depth %s; AUC of the AHP risk for 'flooded' %s" % (
            a["cells_compared"], _fmt(a["spearman_risk_vs_hazard_rating"]), _fmt(a["spearman_risk_vs_depth"]), _fmt(a["roc_auc_risk_vs_wet"])))
        log("  share of each AHP class that HEC-RAS floods (very low .. very high, %%): %s" % ", ".join(_fmt(v, "%.0f") for v in a["wet_share_by_ahp_class_pct"]))
        log("  where the flooded cells sit (share of all flooded cells per AHP class, %%): %s | AUC for deep water (>= 1 m) %s, within 300 m of the centreline %s" % (
            ", ".join(_fmt(v, "%.0f") for v in a["wet_cells_share_by_ahp_class_pct"]), _fmt(a["roc_auc_risk_vs_deep_ge_1m"]), _fmt(a["roc_auc_risk_vs_wet_within_300m_of_centreline"])))
        log("  (AHP ranks where a flood CAN go over all feeders; HEC-RAS is ONE flood from ONE feeder on %s-m terrain: expect the AHP to cover the flooded cells, not match them cell by cell.)" % _fmt(m["cell_m"], "%.0f"))


def score_maps(handoff_dir, export_dir, out_dir, settings, ahp_risk_file="", footprint="", points="", allow_snapshot=False, bc_lines=None, log=print):
    """Score one run. bc_lines: [[(x, y), ...], ...] inflow / outflow lines of the HEC-RAS model (None = those of the handoff). settings: dict(WET_DEPTH_M, POINT_BUFFER_M, OBSERVED_WIDTHS [(chainage_m, width_m)], PROBES [dict], PROBE_RADIUS_M,
    TOLERANCE, DEBRIS_FACTOR, HAZARD_BREAKS). Writes hazard_rating.tif, hazard_class.tif, readback_summary.json, readback_metrics.csv and a
    figure into out_dir; returns the result dict."""
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(handoff_dir, "handoff_summary.json")) as fh:
        summ = json.load(fh)
    cl = {k: np.asarray(v, dtype=float) for k, v in summ["centerline"].items()}
    (dpath, dkind), (vpath, vkind), (apath, akind) = (find_tif(export_dir, "depth", "", allow_snapshot), find_tif(export_dir, "vel", "", allow_snapshot),
                                                      find_tif(export_dir, "arrival", "", True))
    td, tv = (file_time(os.path.basename(dpath)) if dpath else None), (file_time(os.path.basename(vpath)) if vpath else None)
    if td and tv and td != tv:
        raise SystemExit("The depth map is from %02d:%02d:%02d but the velocity map from %02d:%02d:%02d - hazard rating needs both from the same moment (use the Max maps)." % (td + tv))
    snapshot = any(k in ("snapshot", "unlabelled") for k in (dkind, vkind) if k)
    if snapshot:
        log("WARNING: these are NOT maximum maps (%s: %s; %s: %s). Everything below describes ONE moment of the flood, not its extent." % (
            os.path.basename(dpath), dkind, os.path.basename(vpath) if vpath else "no velocity", vkind))
    if dpath is None:
        raise SystemExit("No depth GeoTIFF found in %s (the Max maps are exported by 'fanrisk hecras run')." % export_dir)
    depth, gt, crs = gio.read_raster(dpath)
    depth[depth < 0] = np.nan
    log("depth: %s (%d x %d cells of %.1f m)" % (os.path.basename(dpath), depth.shape[1], depth.shape[0], abs(gt[1])))
    X, Y = grid_xy(gt, depth.shape)
    vel = np.zeros_like(depth)
    if vpath:
        v, vgt, _ = gio.read_raster(vpath)
        vel = np.nan_to_num(v if v.shape == depth.shape else sample_nearest(v, vgt, X, Y), nan=0.0)
    else:
        log("WARNING: no velocity raster found: hazard rating uses depth only (speed taken as 0).")
    arrival = None
    if apath:
        arrival, agt, _ = gio.read_raster(apath)
        if arrival.shape != depth.shape:
            arrival = sample_nearest(arrival, agt, X, Y)
        nm = os.path.basename(apath).lower()                   # RAS Mapper writes arrival time in HOURS unless the name says otherwise
        if not re.search(r"\bmin", nm):
            arrival = arrival * 60.0
        log("arrival time: %s (read as %s, reported in minutes)" % (os.path.basename(apath), "minutes" if re.search(r"\bmin", nm) else "hours"))
    domain = points_in_polygon(X, Y, load_ring(handoff_dir))
    if crs.epsg is None and summ.get("epsg"):
        crs = gio.Crs(epsg=summ["epsg"], wkt=gio.utm_wkt(summ["epsg"]))
    obs_mask = obs_pts = None
    if footprint:
        obs_mask = gio.rasterize([f["geometry"] for f in gio.read_vector(footprint, crs)], depth.shape, gt)
    if points:
        obs_pts = []
        for f in gio.read_vector(points, crs):
            g = f["geometry"]
            obs_pts += [tuple(g["coordinates"][:2])] if g["type"] == "Point" else ([tuple(p[:2]) for p in g["coordinates"]] if g["type"] == "MultiPoint" else [])
    ahp = None
    if ahp_risk_file and os.path.exists(ahp_risk_file):
        a_arr, a_gt, _ = gio.read_raster(ahp_risk_file)
        ahp = (a_arr, a_gt)
    else:
        log("(no AHP risk raster given - skipping the AHP comparison)")
    if bc_lines is None:                                       # the lines of the handoff; pass the model's own lines when they were edited
        bc_lines = []
        for fn in ("bc_inflow.geojson", "bc_outflow.geojson"):
            if os.path.exists(os.path.join(handoff_dir, fn)):
                with open(os.path.join(handoff_dir, fn)) as fh:
                    bc_lines.append(json.load(fh)["features"][0]["geometry"]["coordinates"])
    out, hc, hr = evaluate(depth, vel, arrival, gt, cl, domain, obs_mask, obs_pts, ahp, dict(settings, BC_LINES=bc_lines))
    out["files_used"] = dict(depth=os.path.basename(dpath), depth_kind=dkind, velocity=os.path.basename(vpath) if vpath else None, velocity_kind=vkind,
                             arrival=os.path.basename(apath) if apath else None, maximum_maps=not snapshot)
    report(out, log, settings["WET_DEPTH_M"])
    gio.write_raster(os.path.join(out_dir, "hazard_rating.tif"), np.where(domain, hr, np.nan), gt, crs, "float32", gio.NODATA)
    gio.write_raster(os.path.join(out_dir, "hazard_class.tif"), np.where(domain, hc.astype(float), np.nan), gt, crs, "float32", gio.NODATA)
    with open(os.path.join(out_dir, "readback_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    rows = [("flooded_area_km2", out["model"]["wet_area_km2"]), ("max_depth_m", out["model"]["max_depth_m"]), ("max_speed_ms", out["model"]["max_speed_ms"])]
    for sec in ("footprint", "points", "ahp"):
        for k, v in out.get(sec, {}).items():
            if not isinstance(v, list):
                rows.append((sec + "." + k, v))
    for w in out["widths"]:
        rows.append(("width_%dm.modelled_vs_observed" % w["chainage_m"], "%s vs %s" % (_fmt(w["modelled_contiguous_m"], "%.0f"), _fmt(w["observed_m"], "%.0f"))))
    for p in out["probes"]:
        rows.append(("probe_%s.depth_model_vs_obs" % p["name"], "%s vs %s" % (_fmt(p["depth_model_m"]), _fmt(p["depth_obs_m"]))))
    with open(os.path.join(out_dir, "readback_metrics.csv"), "w") as fh:
        fh.write("metric,value\n" + "\n".join("%s,%s" % (k, v) for k, v in rows) + "\n")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = 3 if ahp is not None else 2
        fig, ax = plt.subplots(1, n, figsize=(7 * n, 6))
        ext = (gt[0], gt[0] + gt[1] * depth.shape[1], gt[3] + gt[5] * depth.shape[0], gt[3])
        dd = depth[domain & np.isfinite(depth)]
        ax[0].imshow(np.where(domain & (np.nan_to_num(depth) > 0), depth, np.nan), extent=ext, cmap="Blues", vmin=0, vmax=max(1.0, np.nanpercentile(dd, 98) if dd.size else 1.0))
        if obs_mask is not None:
            ax[0].contour(np.linspace(ext[0], ext[1], depth.shape[1]), np.linspace(ext[3], ext[2], depth.shape[0]), obs_mask.astype(float), [0.5], colors="red")
        ax[0].plot(cl["x"], cl["y"], "k-", lw=0.8)
        ax[0].set_title("maximum depth (m); red = observed footprint")
        ax[1].imshow(np.where(domain, hc.astype(float), np.nan), extent=ext, cmap="YlOrRd", vmin=0, vmax=4)
        ax[1].set_title("hazard rating class (1 low .. 4 extreme)")
        if ahp is not None:
            ax[2].imshow(np.where(domain, sample_nearest(ahp[0], ahp[1], X, Y), np.nan), extent=ext, cmap="RdYlGn_r", vmin=0, vmax=1)
            ax[2].set_title("AHP risk (absolute scale)")
        for a in ax:
            a.set_aspect("equal")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "readback_figure.png"), dpi=80)
        plt.close(fig)
    except Exception as e:                                   # the figure is a bonus
        log("(figure skipped: %s)" % e)
    log("Written to %s" % out_dir)
    return out


def settings_from_site(site):
    """The scoring settings of a Site: thresholds from [score], observations from its validation case (none without one)."""
    sc, case = site.score, (site.case or {})
    return dict(WET_DEPTH_M=float(sc["wet_depth_m"]), POINT_BUFFER_M=float(sc["point_buffer_m"]), PROBE_RADIUS_M=float(sc["probe_radius_m"]),
                TOLERANCE=float(sc["tolerance"]), DEBRIS_FACTOR=float(sc["debris_factor"]), HAZARD_BREAKS=tuple(sc["hazard_breaks"]),
                OBSERVED_WIDTHS=[(float(w["chainage_m"]), float(w["width_m"])) for w in case.get("widths", [])],
                PROBES=[dict(p) for p in case.get("probes", [])])
