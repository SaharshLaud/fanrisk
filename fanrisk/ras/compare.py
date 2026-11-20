# -*- coding: utf-8 -*-
"""Where do the AHP risk map and the HEC-RAS simulation agree, and which places does the flood hit hardest?

    fanrisk compare <site.toml> --run R04

Inputs: the run's Max maps (depth, velocity, arrival time), the AHP risk raster of the same fan complex, the handoff (model area,
centreline). Outputs, in <run>/compare:
    compare_map.png         4 panels for slides: AHP risk classes | HEC-RAS hazard classes | agreement + numbered hotspots | cross-table
    agreement.tif (+ .qml)  1 both high | 2 HEC-RAS only (AHP rates it lower) | 3 AHP only (this flood does not reach it) | 4 both low but wet
    hotspots.geojson / .kml / .csv
                            the most dangerous places of this flood, numbered, with longitude / latitude, chainage, depth, speed,
                            arrival time and the AHP class there. The .kml opens in Google Earth on satellite imagery.
    crosstab.csv            km2 of every AHP class x HEC-RAS hazard class
"""
import csv
import json
import math
import os
import re

import numpy as np

from .. import gio
from . import readback as rb

AHP_NAMES = ("very low", "low", "moderate", "high", "very high")
HR_NAMES = ("dry", "low", "moderate", "significant", "extreme")
AGREE = {1: ("both high", "#b2182b"), 2: ("HEC-RAS only: AHP rates it lower", "#ef8a62"), 3: ("AHP only: this flood does not reach it", "#67a9cf"),
         4: ("wet, both low", "#d9d9d9")}


def _components(mask):
    """4-connected components of a boolean mask -> label array (0 = background), number of labels."""
    lab = np.zeros(mask.shape, dtype=np.int32)
    H, W = mask.shape
    n = 0
    for r0, c0 in zip(*np.nonzero(mask)):
        if lab[r0, c0]:
            continue
        n += 1
        st = [(r0, c0)]
        lab[r0, c0] = n
        while st:
            r, c = st.pop()
            for rr, cc in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
                if 0 <= rr < H and 0 <= cc < W and mask[rr, cc] and not lab[rr, cc]:
                    lab[rr, cc] = n
                    st.append((rr, cc))
    return lab, n


def _hillshade(z, cell):
    gy, gx = np.gradient(np.nan_to_num(z, nan=np.nanmean(z)), cell)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    az, alt = math.radians(315.0), math.radians(45.0)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(hs, 0, 1)


def _write_qml(path_tif):
    items = "\n".join('     <paletteEntry value="%d" label="%s" color="%s" alpha="255"/>' % (k, n, c) for k, (n, c) in AGREE.items())
    txt = ("<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n<qgis version=\"3.34.0\" styleCategories=\"Symbology\">\n <pipe>\n"
           "  <rasterrenderer type=\"paletted\" band=\"1\" opacity=\"0.85\" alphaBand=\"-1\">\n   <rasterTransparency/>\n   <colorPalette>\n"
           + items + "\n   </colorPalette>\n  </rasterrenderer>\n </pipe>\n</qgis>\n")
    with open(os.path.splitext(path_tif)[0] + ".qml", "w") as fh:
        fh.write(txt)


def _kml(path, rows, title):
    pm = []
    for h in rows:
        desc = ("Hazard rating %.1f (%s)&lt;br&gt;max depth %.1f m, max speed %.1f m/s&lt;br&gt;water arrives after %s min&lt;br&gt;"
                "area %.0f m2, %.1f km down the flood path&lt;br&gt;AHP class here: %s" % (
                    h["hazard_rating_max"], h["hazard_class"], h["depth_max_m"], h["speed_max_ms"],
                    "-" if h["arrival_min"] is None else "%.0f" % h["arrival_min"], h["area_m2"], h["chainage_km"], h["ahp_class"]))
        pm.append("<Placemark><name>%d</name><description>%s</description><Point><coordinates>%.6f,%.6f,0</coordinates></Point></Placemark>" % (
            h["rank"], desc, h["lon"], h["lat"]))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>%s</name>\n%s\n</Document></kml>\n'
                 % (title, "\n".join(pm)))


def compare(handoff_dir, maps_dir, ahp_file, out_dir, dem_file="", ahp_breaks=(0.2, 0.4, 0.55, 0.7), hr_breaks=(0.75, 1.25, 2.5), wet_depth=0.3,
            debris_factor=0.5, n_hotspots=8, segment_m=500.0, title="", log=print):
    os.makedirs(out_dir, exist_ok=True)
    summ = json.load(open(os.path.join(handoff_dir, "handoff_summary.json")))
    (dpath, _), (vpath, _), (apath, _) = rb.find_tif(maps_dir, "depth", "", False), rb.find_tif(maps_dir, "vel", "", False), rb.find_tif(maps_dir, "arrival", "", True)
    if dpath is None:
        raise SystemExit("no Max depth map in %s" % maps_dir)
    depth, gt, crs = gio.read_raster(dpath)
    if crs.epsg is None and summ.get("epsg"):
        crs = gio.Crs(epsg=summ["epsg"], wkt=gio.utm_wkt(summ["epsg"]))
    X, Y = rb.grid_xy(gt, depth.shape)
    vel = np.zeros_like(depth)
    if vpath:
        v, vgt, _ = gio.read_raster(vpath)
        vel = np.nan_to_num(v if v.shape == depth.shape else rb.sample_nearest(v, vgt, X, Y))
    arr = None
    if apath:
        a, agt, _ = gio.read_raster(apath)
        arr = (a if a.shape == depth.shape else rb.sample_nearest(a, agt, X, Y))
        if not re.search(r"\bmin", os.path.basename(apath).lower()):
            arr = arr * 60.0
    domain = rb.points_in_polygon(X, Y, rb.load_ring(handoff_dir))
    d0 = np.nan_to_num(depth)
    wet = (d0 >= wet_depth) & domain
    hr = rb.hazard_rating(d0, vel, debris_factor)
    hc = rb.hazard_class(hr, wet, hr_breaks)                       # 0 dry .. 4 extreme
    a_arr, a_gt, _ = gio.read_raster(ahp_file)
    ahp = rb.sample_nearest(a_arr, a_gt, X, Y)
    ahp_ok = domain & np.isfinite(ahp)
    ac = np.where(ahp_ok, np.digitize(np.nan_to_num(ahp), ahp_breaks), -1)   # 0 very low .. 4 very high
    ahp_high = ahp_ok & (ac >= 3)
    ras_high = wet & (hc >= 3)
    agree = np.zeros(depth.shape, dtype=np.uint8)
    agree[wet] = 4
    agree[ahp_high & ~wet] = 3
    agree[ras_high & ~ahp_high] = 2
    agree[ras_high & ahp_high] = 1
    cell_km2 = abs(gt[1] * gt[5]) / 1e6
    area = {AGREE[k][0]: float((agree == k).sum() * cell_km2) for k in AGREE}
    ct = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            ct[i, j] = ((ac == i) & (hc == j) & ahp_ok).sum() * cell_km2
    with open(os.path.join(out_dir, "crosstab.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["AHP class \\ HEC-RAS hazard (km2)"] + list(HR_NAMES))
        for i in range(5):
            w.writerow([AHP_NAMES[i]] + ["%.3f" % v for v in ct[i]])
    gio.write_raster(os.path.join(out_dir, "agreement.tif"), agree, gt, crs, "uint8", 0)
    _write_qml(os.path.join(out_dir, "agreement.tif"))

    # ---- hotspots: connected patches of significant / extreme hazard, ranked by their worst hazard x size ----------
    s_cl, _, _ = rb.polyline_frame(np.asarray(summ["centerline"]["x"]), np.asarray(summ["centerline"]["y"]), X, Y)
    lab, n = _components(ras_high)
    spots = []

    def add(mm):
        if mm.sum() < 2:
            return
        i = np.argmax(np.where(mm, hr, -1))
        r, c = np.unravel_index(i, hr.shape)
        spots.append(dict(cells=mm, r=r, c=c, score=float(hr[r, c]) * math.sqrt(mm.sum())))
    for k in range(1, n + 1):
        m = lab == k
        span = s_cl[m].max() - s_cl[m].min() if m.any() else 0.0
        if span <= 1.2 * segment_m:
            add(m)
            continue
        edges = np.arange(s_cl[m].min(), s_cl[m].max() + segment_m, segment_m)   # a long patch along the channel: one spot per segment
        for a0, a1 in zip(edges[:-1], edges[1:]):
            add(m & (s_cl >= a0) & (s_cl < a1))
    spots = sorted(spots, key=lambda s: -s["score"])[:n_hotspots]
    spots = sorted(spots, key=lambda s: s_cl[s["r"], s["c"]])          # number them down the flood path
    lon, lat = gio.transform_xy([X[s["r"], s["c"]] for s in spots], [Y[s["r"], s["c"]] for s in spots], crs, gio.WGS84) if spots else ([], [])
    rows = []
    for q, s in enumerate(spots):
        r, c, m = s["r"], s["c"], s["cells"]
        rows.append(dict(rank=q + 1, x=float(X[r, c]), y=float(Y[r, c]), lon=float(lon[q]), lat=float(lat[q]), chainage_km=float(s_cl[r, c] / 1000.0),
                         area_m2=float(m.sum() * cell_km2 * 1e6), hazard_rating_max=float(hr[r, c]), hazard_class=HR_NAMES[int(hc[r, c])],
                         depth_max_m=float(d0[m].max()), speed_max_ms=float(vel[m].max()),
                         arrival_min=None if arr is None or not np.isfinite(arr[m]).any() else float(np.nanmin(arr[m])),
                         ahp_risk=None if not np.isfinite(ahp[r, c]) else float(ahp[r, c]),
                         ahp_class=AHP_NAMES[int(ac[r, c])] if ac[r, c] >= 0 else "outside AHP map",
                         agreement=AGREE[int(agree[r, c])][0]))
    if rows:
        with open(os.path.join(out_dir, "hotspots.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        layer = dict(name="hotspots", fields=[(k, "double" if isinstance(v, float) else "integer" if isinstance(v, int) else "string")
                                              for k, v in rows[0].items() if k not in ("x", "y")],
                     rows=[({"type": "Point", "coordinates": [h["x"], h["y"]]}, [h[k] for k in rows[0] if k not in ("x", "y")]) for h in rows])
        gio.write_geojson(os.path.join(out_dir, "hotspots.geojson"), layer, crs)
        _kml(os.path.join(out_dir, "hotspots.kml"), rows, title or "fanrisk hotspots")

    # ---- figure ---------------------------------------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
        ext = (gt[0], gt[0] + gt[1] * depth.shape[1], gt[3] + gt[5] * depth.shape[0], gt[3])
        hs = None
        if dem_file and os.path.exists(dem_file):
            z = rb.sample_nearest(*gio.read_raster(dem_file)[:2], X, Y)
            hs = _hillshade(z, abs(gt[1]))
        ring = np.asarray(rb.load_ring(handoff_dir))
        cl = summ["centerline"]
        fig = plt.figure(figsize=(19, 12))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.28, wspace=0.06)
        ahp_cmap = ListedColormap(["#1a9641", "#a6d96a", "#ffffbf", "#fdae61", "#d7191c"])
        hr_cmap = ListedColormap(["#fee5d9", "#fcae91", "#fb6a4a", "#a50f15"])
        ag_cmap = ListedColormap([AGREE[k][1] for k in (1, 2, 3, 4)])
        panels = [("a  AHP risk (absolute scale)", np.where(ahp_ok, ac, np.nan), ahp_cmap, -0.5, 4.5,
                   [Patch(color=ahp_cmap(i), label=AHP_NAMES[i]) for i in range(5)]),
                  ("b  HEC-RAS hazard rating (depth x speed)", np.where(wet, hc, np.nan), hr_cmap, 0.5, 4.5,
                   [Patch(color=hr_cmap(i), label=HR_NAMES[i + 1]) for i in range(4)]),
                  ("c  Agreement, and the hotspots of this flood", np.where(agree > 0, agree, np.nan), ag_cmap, 0.5, 4.5,
                   [Patch(color=AGREE[k][1], label=AGREE[k][0]) for k in (1, 2, 3, 4)])]
        for j, (ttl, img, cmap, lo, hi, leg) in enumerate(panels):
            ax = fig.add_subplot(gs[0, j])
            if hs is not None:
                ax.imshow(hs, extent=ext, cmap="Greys_r", vmin=0, vmax=1.3, alpha=0.9)
            ax.imshow(img, extent=ext, cmap=cmap, vmin=lo, vmax=hi, alpha=0.85, interpolation="nearest")
            ax.plot(ring[:, 0], ring[:, 1], color="#444", lw=0.7)
            ax.plot(cl["x"], cl["y"], color="#222", lw=0.6, ls="--")
            if j == 2:
                for h in rows:
                    ax.plot(h["x"], h["y"], "o", ms=15, mfc="white", mec="#111", mew=1.3)
                    ax.text(h["x"], h["y"], str(h["rank"]), ha="center", va="center", fontsize=9, weight="bold", color="#111")
            ax.set_xlim(ring[:, 0].min() - 100, ring[:, 0].max() + 100)
            ax.set_ylim(ring[:, 1].min() - 100, ring[:, 1].max() + 100)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(ttl, fontsize=12, loc="left")
            ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2 if len(leg) > 3 else len(leg), fontsize=9, frameon=False)
        ax = fig.add_subplot(gs[1, 0])
        tot = ct.sum(axis=1)
        left = np.zeros(5)
        cols = ["#e0e0e0"] + [hr_cmap(i) for i in range(4)]
        for j in range(5):
            share = np.where(tot > 0, 100 * ct[:, j] / np.where(tot > 0, tot, 1), 0)
            ax.barh(range(5), share, left=left, color=cols[j], edgecolor="white", linewidth=2, label=HR_NAMES[j], height=0.65)
            left += share
        for i in range(5):
            ax.text(102, i, "%.2f km2" % tot[i], va="center", fontsize=9, color="#555")
        ax.set_yticks(range(5))
        ax.set_yticklabels(AHP_NAMES)
        ax.set_ylabel("AHP class")
        ax.set_xlim(0, 120)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("share of the AHP class's area (%)")
        ax.set_title("d  What the simulated flood does in each AHP class", fontsize=12, loc="left")
        ax.legend(title="HEC-RAS hazard", fontsize=8.5, title_fontsize=9, loc="upper center", bbox_to_anchor=(0.45, -0.18), ncol=5, frameon=False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax = fig.add_subplot(gs[1, 1:])
        ax.axis("off")
        ax.set_title("e  Hotspots (numbered down the flood path)", fontsize=12, loc="left")
        hdr = ["#", "longitude, latitude", "km down\nthe path", "max depth\n(m)", "max speed\n(m/s)", "water arrives\n(min)", "HEC-RAS\nhazard", "AHP\nclass", "agreement"]
        cells = [[str(h["rank"]), "%.5f E, %.5f N" % (h["lon"], h["lat"]), "%.1f" % h["chainage_km"], "%.1f" % h["depth_max_m"], "%.1f" % h["speed_max_ms"],
                  "-" if h["arrival_min"] is None else "%.0f" % h["arrival_min"], h["hazard_class"], h["ahp_class"], h["agreement"].split(":")[0]] for h in rows]
        if cells:
            tb = ax.table(cellText=cells, colLabels=hdr, loc="upper left", cellLoc="center", colLoc="center",
                          colWidths=[0.04, 0.2, 0.08, 0.08, 0.08, 0.1, 0.1, 0.11, 0.17])
            tb.auto_set_font_size(False)
            tb.set_fontsize(9.5)
            tb.scale(1, 1.9)
            for (r_, c_), cell in tb.get_celld().items():
                cell.set_edgecolor("#dddddd")
                if r_ == 0:
                    cell.set_facecolor("#f2f2f2")
                    cell.set_text_props(weight="bold")
                elif c_ == 8:
                    k = next((k for k, (nm, _c) in AGREE.items() if nm.startswith(cells[r_ - 1][8])), None)
                    if k:
                        cell.set_facecolor(AGREE[k][1])
                        cell.set_text_props(color="white" if k in (1, 2) else "#222")
        if title:
            fig.suptitle(title, x=0.125, y=0.95, ha="left", fontsize=14, weight="bold")
        plt.savefig(os.path.join(out_dir, "compare_map.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        log("(figure skipped: %s)" % e)

    out = dict(area_km2=area, crosstab_km2=dict(rows=AHP_NAMES, cols=HR_NAMES, values=ct.round(4).tolist()),
               hit_rate_of_ahp_high_pct=100.0 * float((ahp_high & wet).sum()) / max(1, int(ahp_high.sum())),
               share_of_ras_high_rated_high_by_ahp_pct=100.0 * float((ras_high & ahp_high).sum()) / max(1, int(ras_high.sum())), hotspots=rows)
    with open(os.path.join(out_dir, "compare_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    log("Agreement (km2): " + ", ".join("%s %.2f" % kv for kv in area.items()))
    log("  of the places HEC-RAS rates significant / extreme, AHP rates %.0f%% high or very high" % out["share_of_ras_high_rated_high_by_ahp_pct"])
    log("  of the places AHP rates high / very high, this flood wets %.0f%%" % out["hit_rate_of_ahp_high_pct"])
    log("Hotspots (numbered down the flood path):")
    for h in rows:
        log("  %d  %.5f E %.5f N  %4.1f km  depth %.1f m  speed %.1f m/s  arrives %s min  HEC-RAS %s, AHP %s" % (
            h["rank"], h["lon"], h["lat"], h["chainage_km"], h["depth_max_m"], h["speed_max_ms"],
            "-" if h["arrival_min"] is None else "%.0f" % h["arrival_min"], h["hazard_class"], h["ahp_class"]))
    log("Written to %s (compare_map.png, agreement.tif, hotspots.kml / .geojson / .csv, crosstab.csv)" % out_dir)
    return out
