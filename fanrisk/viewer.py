# -*- coding: utf-8 -*-
"""One HTML page that shows a site's results on a web map (opens in any browser; needs internet only for the map library and basemaps).

    fanrisk viewer <site.toml> [--run R04 | --maps <HEC-RAS Max maps folder>] [--out <file.html>]

Layers (each only when its file exists):
    AHP risk classes (all complexes)             runs/<site>/risk/ahp_C*_risk_abs.tif
    HEC-RAS max depth and hazard class           runs/<site>/hecras/runs/<run>/maps, .../readback/hazard_class.tif
    observed flood footprint + damage points     [validation] flood_footprint / flood_points
    buildings and critical facilities hit        runs/<site>/exposure/exposure_buildings.geojson (from 'fanrisk exposure')
    hotspots                                     runs/<site>/hecras/runs/<run>/compare/hotspots.geojson (from 'fanrisk compare')
Rasters are resampled to longitude/latitude (nearest cell) and embedded as PNG images, so the page is one self-contained file.
"""
import base64
import glob
import json
import os
import struct
import zlib

import numpy as np

from . import gio

WGS = gio.Crs(epsg=4326, geographic=True)


def _png(rgba):
    """RGBA uint8 array -> PNG bytes (no imaging library needed)."""
    h, w, _ = rgba.shape
    raw = b"".join(b"\x00" + rgba[r].tobytes() for r in range(h))

    def chunk(t, d):
        c = struct.pack(">I", len(d)) + t + d
        return c + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b"")


def _hex(c):
    c = c.lstrip("#")
    return tuple(int(c[i:i + 2], 16) for i in (0, 2, 4))


def _to_lonlat(arr, gt, crs, res_deg=None):
    """Resample a projected raster to a lon/lat grid. -> (array, [[south, west], [north, east]])."""
    fin = np.isfinite(arr)
    if not fin.any():
        return None, None
    r, c = np.where(fin)
    xs = gt[0] + np.array([c.min(), c.max() + 1, c.min(), c.max() + 1]) * gt[1]
    ys = gt[3] + np.array([r.min(), r.min(), r.max() + 1, r.max() + 1]) * gt[5]
    lon, lat = gio.transform_xy(xs, ys, crs, WGS) if not crs.geographic else (xs, ys)
    w, e, s, n = lon.min(), lon.max(), lat.min(), lat.max()
    res = res_deg or abs(gt[1]) / 111000.0 * 0.5          # half a cell: keeps the cell edges sharp
    nx, ny = int(np.ceil((e - w) / res)), int(np.ceil((n - s) / res))
    while nx * ny > 4e6:
        res *= 1.5
        nx, ny = int(np.ceil((e - w) / res)), int(np.ceil((n - s) / res))
    LON, LAT = np.meshgrid(w + (np.arange(nx) + 0.5) * res, n - (np.arange(ny) + 0.5) * res)
    X, Y = gio.transform_xy(LON.ravel(), LAT.ravel(), WGS, crs) if not crs.geographic else (LON.ravel(), LAT.ravel())
    cc = np.floor((X - gt[0]) / gt[1]).astype(int)
    rr = np.floor((Y - gt[3]) / gt[5]).astype(int)
    ok = (rr >= 0) & (rr < arr.shape[0]) & (cc >= 0) & (cc < arr.shape[1])
    out = np.full(X.shape, np.nan)
    out[ok] = arr[rr[ok], cc[ok]]
    return out.reshape(ny, nx), [[s, w], [n, e]]


def _classes_png(a, breaks, colors, alpha=200):
    """a: values; breaks: upper edges of all but the last class; colors: one per class (None = transparent)."""
    rgba = np.zeros(a.shape + (4,), np.uint8)
    k = np.digitize(np.nan_to_num(a, nan=-1e9), breaks)
    for i, col in enumerate(colors):
        if col is None:
            continue
        m = (k == i) & np.isfinite(a)
        rgba[m, :3] = _hex(col)
        rgba[m, 3] = alpha
    return "data:image/png;base64," + base64.b64encode(_png(rgba)).decode()


def _gj_lonlat(features, crs, props=None):
    def conv(coords):
        a = np.asarray(coords, float)
        if a.ndim == 1:
            lo, la = gio.transform_xy([a[0]], [a[1]], crs, WGS)
            return [round(float(lo[0]), 6), round(float(la[0]), 6)]
        if a.ndim == 2:
            lo, la = gio.transform_xy(a[:, 0], a[:, 1], crs, WGS)
            return [[round(float(x), 6), round(float(y), 6)] for x, y in zip(lo, la)]
        return [conv(c) for c in coords]
    out = []
    for f in features:
        g = f.get("geometry")
        if not g:
            continue
        p = f.get("properties") or {}
        out.append(dict(type="Feature", geometry=dict(type=g["type"], coordinates=conv(g["coordinates"]) if not crs.geographic else g["coordinates"]),
                        properties={k: p.get(k) for k in (props or p.keys())}))
    return dict(type="FeatureCollection", features=out)


def build(site, maps_dir="", run_name="R04", out_path="", log=print):
    from .ras import readback as rb
    _dem, gt_dem, crs, _nd = gio.load_dem(site.file("dem"), site.P["REPROJECT_GEOGRAPHIC"], lambda m: None)
    runs = site.file("runs_dir") or os.path.join(site.out_dir, "hecras", "runs")
    run_dir = os.path.dirname(os.path.abspath(maps_dir)) if maps_dir else os.path.join(runs, run_name)
    md = maps_dir or os.path.join(run_dir, "maps")
    images, vectors, legends, bounds = [], [], [], []

    # AHP risk, all complexes on the DEM grid
    fs = sorted(glob.glob(os.path.join(site.out_dir, "risk", "ahp_C*_risk_abs.tif")))
    if fs:
        best = np.full(_dem.shape, np.nan)
        X, Y = rb.grid_xy(gt_dem, _dem.shape)
        for f in fs:
            a, g, _ = gio.read_raster(f)
            s = rb.sample_nearest(a, g, X, Y)
            best = np.where(np.isfinite(s) & ~(s <= best), s, best)
        ll, bb = _to_lonlat(best, gt_dem, crs)
        if ll is not None:
            br = list(site.P["RISK_CLASS_BREAKS"])
            cols = ["#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"]
            names = ["very low", "low", "moderate", "high", "very high"]
            # very low and low are left transparent: they cover most of the valley and hide the imagery
            images.append(dict(name="AHP risk, moderate and above (all fan complexes)", url=_classes_png(ll, br, [None, None] + cols[2:], 160), bounds=bb, on=True))
            legends.append(dict(title="AHP risk", items=list(zip(names[2:], cols[2:]))))
            bounds.append(bb)
    # HEC-RAS
    if os.path.isdir(md):
        p, _k = rb.find_tif(md, "depth", "", False)
        if p:
            a, g, _ = gio.read_raster(p)
            ll, bb = _to_lonlat(np.where(np.nan_to_num(a) > 0.01, a, np.nan), g, crs)
            if ll is not None:
                cols = ["#c6dbef", "#6baed6", "#2171b5", "#08306b"]
                images.append(dict(name="HEC-RAS max depth (%s)" % os.path.basename(run_dir), url=_classes_png(ll, [0.3, 1.0, 2.0], cols, 200), bounds=bb, on=False))
                legends.append(dict(title="Max depth", items=list(zip(["< 0.3 m", "0.3-1 m", "1-2 m", "> 2 m"], cols))))
                bounds.append(bb)
    hp = os.path.join(run_dir, "readback", "hazard_class.tif")
    if os.path.exists(hp):
        a, g, _ = gio.read_raster(hp)
        ll, bb = _to_lonlat(np.where(np.nan_to_num(a) >= 1, a, np.nan), g, crs)
        if ll is not None:
            cols = ["#d1e5f0", "#92c5de", "#f4a582", "#b2182b"]
            images.append(dict(name="HEC-RAS hazard class (%s)" % os.path.basename(run_dir), url=_classes_png(ll, [1.5, 2.5, 3.5], cols, 210), bounds=bb, on=True))
            legends.append(dict(title="Hazard (depth x velocity)", items=list(zip(["low", "moderate", "significant", "extreme"], cols))))
    # observed flood
    if site.file("flood_footprint"):
        vectors.append(dict(name="2010 observed flood footprint", kind="outline", on=True,
                            data=_gj_lonlat(gio.read_vector(site.file("flood_footprint"), crs), crs, [])))
    if site.file("flood_points"):
        vectors.append(dict(name="2010 damage points (Bhatt et al. 2011)", kind="dots", on=False,
                            data=_gj_lonlat(gio.read_vector(site.file("flood_points"), crs), crs, [])))
    # exposure
    ep = os.path.join(site.out_dir, "exposure", "exposure_buildings.geojson")
    if os.path.exists(ep):
        doc = json.load(open(ep))
        feats = doc["features"]
        b = [f for f in feats if f["properties"].get("kind") == "building"]
        c = [f for f in feats if f["properties"].get("kind") != "building"]
        keep = ["kind", "name", "ahp_risk", "hec_class", "hec_depth_m", "in_2010_footprint", "area_m2"]
        vectors.append(dict(name="Buildings hit (%d)" % len(b), kind="buildings", on=True, data=_gj_lonlat(b, crs, keep)))
        if c:
            vectors.append(dict(name="Critical facilities hit (%d)" % len(c), kind="facilities", on=True, data=_gj_lonlat(c, crs, keep)))
    hs = os.path.join(run_dir, "compare", "hotspots.geojson")
    if os.path.exists(hs):
        vectors.append(dict(name="Hotspots (AHP vs HEC-RAS)", kind="hotspots", on=False, data=_gj_lonlat(gio.read_vector(hs, crs), crs, None)))
    if not images and not vectors:
        raise SystemExit("Nothing to show yet: run 'fanrisk run' (and HEC-RAS / exposure) first.")
    if not bounds:
        bounds = [[[0, 0], [0, 0]]]
    s = min(b[0][0] for b in bounds)
    w = min(b[0][1] for b in bounds)
    n = max(b[1][0] for b in bounds)
    e = max(b[1][1] for b in bounds)
    fp = [v for v in vectors if v["kind"] == "outline"]
    if fp and fp[0]["data"]["features"]:                  # open on the observed flood when there is one
        pts = np.array([p for f in fp[0]["data"]["features"] for ring in f["geometry"]["coordinates"] for p in (ring if f["geometry"]["type"] == "Polygon" else ring[0])])
        s, w, n, e = pts[:, 1].min(), pts[:, 0].min(), pts[:, 1].max(), pts[:, 0].max()
        d = max(n - s, e - w)
        s, w, n, e = s - d, w - d, n + d * 0.6, e + d * 0.6
    out_path = out_path or os.path.join(site.out_dir, "viewer.html")
    html = TEMPLATE.replace("__TITLE__", site.name).replace("__DATA__", json.dumps(dict(images=images, vectors=vectors, legends=legends, view=[[s, w], [n, e]])))
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    log("viewer: %s (%.1f MB) - %d raster layers, %d vector layers. Open it in a browser." % (out_path, os.path.getsize(out_path) / 1e6, len(images), len(vectors)))
    return out_path


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Flood risk viewer</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
 html,body{margin:0;height:100%;font:13px/1.4 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:#fff;color:#1a1a1a}
 #map{position:absolute;top:44px;bottom:0;left:0;right:0}
 header{height:44px;display:flex;align-items:center;gap:12px;padding:0 16px;border-bottom:1px solid #ddd;box-sizing:border-box}
 header b{font-size:14px} header span{color:#666;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
 .legend{background:rgba(255,255,255,.94);padding:8px 10px;border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.25);max-width:190px}
 .legend h4{margin:6px 0 3px;font-size:12px} .legend h4:first-child{margin-top:0}
 .legend i{display:inline-block;width:14px;height:10px;margin-right:6px;border-radius:2px;vertical-align:middle}
 .leaflet-popup-content{font-size:12px} .leaflet-popup-content td{padding:1px 6px 1px 0}
 #op{position:absolute;bottom:18px;left:10px;z-index:999;background:rgba(255,255,255,.94);padding:6px 8px;border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.25)}
</style></head><body>
<header><b>__TITLE__</b><span>AHP screening, HEC-RAS 2D simulation, the 2010 flood and who is exposed</span></header>
<div id="map"></div>
<div id="op">Raster opacity <input id="opacity" type="range" min="0" max="100" value="85"></div>
<script>
const D = __DATA__;
const HZ = ["dry","low","moderate","significant","extreme"];
const map = L.map('map', {preferCanvas:true});
const sat = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {maxZoom:19, attribution:'Imagery &copy; Esri'});
// OpenStreetMap's own tile servers refuse pages opened from a file (no referrer), so the street map comes from Esri too
const osm = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}', {maxZoom:19, attribution:'Map &copy; Esri, HERE, Garmin, &copy; OpenStreetMap contributors'});
sat.addTo(map);
const overlays = {}, rasters = [];
D.images.forEach(im => { const l = L.imageOverlay(im.url, im.bounds, {opacity:.85}); rasters.push(l); overlays[im.name] = l; if (im.on) l.addTo(map); });
function fmt(v, d){ return (v===null||v===undefined) ? '-' : (typeof v==='number' ? v.toFixed(d) : v); }
function popup(p){
  return '<table>' + (p.name ? '<tr><td colspan=2><b>'+p.name+'</b></td></tr>' : '') +
    '<tr><td>type</td><td>'+(p.kind||'')+'</td></tr>' +
    '<tr><td>HEC-RAS hazard</td><td>'+(p.hec_class===null||p.hec_class===undefined ? '-' : HZ[Math.round(p.hec_class)])+'</td></tr>' +
    '<tr><td>max depth</td><td>'+fmt(p.hec_depth_m,2)+' m</td></tr>' +
    '<tr><td>AHP risk</td><td>'+fmt(p.ahp_risk,2)+'</td></tr>' +
    '<tr><td>in 2010 flood</td><td>'+(p.in_2010_footprint ? 'yes' : 'no')+'</td></tr></table>';
}
const hzCol = ["#999","#4393c3","#92c5de","#f4a582","#b2182b"];
D.vectors.forEach(v => {
  let l;
  if (v.kind === 'outline') l = L.geoJSON(v.data, {style:{color:'#000', weight:2, fill:false, dashArray:'5 3'}});
  else if (v.kind === 'dots') l = L.geoJSON(v.data, {pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:2.5,color:'#6a3d9a',weight:0,fillOpacity:.9})});
  else if (v.kind === 'buildings') l = L.geoJSON(v.data, {pointToLayer:(f,ll)=>{const p=f.properties,h=Math.round(p.hec_class||0);
      return L.circleMarker(ll,{radius:3.5,color:'#fff',weight:.6,fillColor:h>=1?hzCol[h]:(p.ahp_risk>=0.55?'#fd8d3c':'#555'),fillOpacity:.95});},
      onEachFeature:(f,x)=>x.bindPopup(popup(f.properties))});
  else if (v.kind === 'facilities') l = L.geoJSON(v.data, {pointToLayer:(f,ll)=>L.marker(ll,{icon:L.divIcon({className:'',html:'<div style="width:14px;height:14px;background:#2166ac;border:2px solid #fff;transform:rotate(45deg);box-shadow:0 0 3px #000"></div>',iconSize:[14,14]})}),
      onEachFeature:(f,x)=>x.bindPopup(popup(f.properties))});
  else l = L.geoJSON(v.data, {pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:7,color:'#000',weight:1.5,fillColor:'#ffff33',fillOpacity:.9}),
      onEachFeature:(f,x)=>x.bindPopup('<pre style="margin:0">'+JSON.stringify(f.properties,null,1)+'</pre>')});
  overlays[v.name] = l; if (v.on) l.addTo(map);
});
L.control.layers({"Satellite (Esri)":sat, "Street map (Esri)":osm}, overlays, {collapsed:false}).addTo(map);
const lg = L.control({position:'bottomright'});
lg.onAdd = () => { const d = L.DomUtil.create('div','legend');
  d.innerHTML = D.legends.map(g => '<h4>'+g.title+'</h4>' + g.items.map(i => '<div><i style="background:'+i[1]+'"></i>'+i[0]+'</div>').join('')).join('') +
    '<h4>Buildings hit</h4><div><i style="background:#b2182b"></i>by HEC-RAS hazard</div><div><i style="background:#fd8d3c"></i>AHP high only</div>' +
    '<h4>Other</h4><div><i style="background:#2166ac"></i>critical facility</div><div><i style="border:2px dashed #000;height:6px;width:10px"></i>2010 flood</div>';
  return d; };
lg.addTo(map);
map.fitBounds(D.view);
L.control.scale({imperial:false}).addTo(map);
document.getElementById('opacity').oninput = e => rasters.forEach(r => r.setOpacity(e.target.value/100));
</script></body></html>
"""
