# -*- coding: utf-8 -*-
"""Geo file input / output for fanrisk.

Two backends, chosen automatically (set the environment variable FANRISK_IO=basic to force the second one):

    "rasterio"  rasterio (+ geopandas when installed): any raster GDAL reads, reprojection, GeoPackage and shapefile output.
    "basic"     tifffile + pure Python: GeoTIFFs that are already in a projected (UTM) CRS, GeoJSON vectors in / out. No reprojection of
                rasters. Enough to run the whole pipeline on a UTM DEM; used by the tests.

Geometries travel as GeoJSON-like dicts ({"type": "Polygon", "coordinates": ...}); rasters as numpy arrays plus a GDAL-order geotransform
(x0, dx, 0, y0, 0, -dy) and a Crs.
"""
import json
import math
import os
import re

import numpy as np

NODATA = -9999.0


class Crs(object):
    """A coordinate reference system: WKT text and / or EPSG code."""

    def __init__(self, wkt="", epsg=None, geographic=False, geokeys=None):
        self.wkt, self.epsg, self.geographic, self.geokeys = wkt or "", (int(epsg) if epsg else None), bool(geographic), geokeys

    def authid(self):
        return "EPSG:%d" % self.epsg if self.epsg else ""

    def same(self, other):
        if self.epsg and other.epsg:
            return self.epsg == other.epsg
        return bool(self.wkt) and self.wkt == other.wkt

    def __repr__(self):
        return "Crs(%s)" % (self.authid() or (self.wkt[:40] + "..."))


WGS84 = Crs(epsg=4326, geographic=True)
_BACKEND = [None]


def backend():
    if _BACKEND[0] is None:
        want = os.environ.get("FANRISK_IO", "").lower()
        name = "basic"
        if want != "basic":
            try:
                import rasterio  # noqa: F401
                name = "rasterio"
            except Exception:
                name = "basic"
        if name == "basic":
            try:
                import tifffile  # noqa: F401
            except Exception:
                raise RuntimeError("fanrisk needs a raster library: install rasterio (pip install rasterio geopandas) or, at least, tifffile.")
        _BACKEND[0] = name
    return _BACKEND[0]


def have(module):
    try:
        __import__(module)
        return True
    except Exception:
        return False


def utm_epsg_for(lon, lat):
    return (32600 if lat >= 0 else 32700) + int((lon + 180.0) // 6.0) + 1


# =============================================================================
# pure-Python UTM <-> WGS84 (Krueger series; used by the basic backend)
# =============================================================================
_A, _F, _K0 = 6378137.0, 1.0 / 298.257223563, 0.9996
_N = _F / (2.0 - _F)
_AA = _A / (1.0 + _N) * (1.0 + _N ** 2 / 4.0 + _N ** 4 / 64.0)
_AL = (_N / 2.0 - 2.0 * _N ** 2 / 3.0 + 5.0 * _N ** 3 / 16.0, 13.0 * _N ** 2 / 48.0 - 3.0 * _N ** 3 / 5.0, 61.0 * _N ** 3 / 240.0)
_BE = (_N / 2.0 - 2.0 * _N ** 2 / 3.0 + 37.0 * _N ** 3 / 96.0, _N ** 2 / 48.0 + _N ** 3 / 15.0, 17.0 * _N ** 3 / 480.0)
_DE = (2.0 * _N - 2.0 * _N ** 2 / 3.0 - 2.0 * _N ** 3, 7.0 * _N ** 2 / 3.0 - 8.0 * _N ** 3 / 5.0, 56.0 * _N ** 3 / 15.0)


def _utm_zone(epsg):
    if epsg and 32601 <= epsg <= 32660:
        return epsg - 32600, False
    if epsg and 32701 <= epsg <= 32760:
        return epsg - 32700, True
    return None, None


def _utm_fwd(lon, lat, zone, south):
    lon, lat = np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    phi, lam = np.radians(lat), np.radians(lon) - math.radians(zone * 6.0 - 183.0)
    e = 2.0 * math.sqrt(_N) / (1.0 + _N)
    t = np.sinh(np.arctanh(np.sin(phi)) - e * np.arctanh(e * np.sin(phi)))
    xi, eta = np.arctan2(t, np.cos(lam)), np.arctanh(np.sin(lam) / np.sqrt(1.0 + t * t))
    E, Nn = eta.copy(), xi.copy()
    for j, a in enumerate(_AL, 1):
        E = E + a * np.cos(2 * j * xi) * np.sinh(2 * j * eta)
        Nn = Nn + a * np.sin(2 * j * xi) * np.cosh(2 * j * eta)
    return 500000.0 + _K0 * _AA * E, (1.0e7 if south else 0.0) + _K0 * _AA * Nn


def _utm_inv(x, y, zone, south):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    xi, eta = (y - (1.0e7 if south else 0.0)) / (_K0 * _AA), (x - 500000.0) / (_K0 * _AA)
    xi1, eta1 = xi.copy(), eta.copy()
    for j, b in enumerate(_BE, 1):
        xi1 = xi1 - b * np.sin(2 * j * xi) * np.cosh(2 * j * eta)
        eta1 = eta1 - b * np.cos(2 * j * xi) * np.sinh(2 * j * eta)
    chi = np.arcsin(np.sin(xi1) / np.cosh(eta1))
    phi = chi.copy()
    for j, d in enumerate(_DE, 1):
        phi = phi + d * np.sin(2 * j * chi)
    lam = math.radians(zone * 6.0 - 183.0) + np.arctan2(np.sinh(eta1), np.cos(xi1))
    return np.degrees(lam), np.degrees(phi)


def transform_xy(xs, ys, src, dst):
    """Arrays of x, y from Crs `src` to Crs `dst` (x = easting / longitude)."""
    xs, ys = np.atleast_1d(np.asarray(xs, dtype=float)), np.atleast_1d(np.asarray(ys, dtype=float))
    if src.same(dst) or xs.size == 0:
        return xs, ys
    zs, ss = _utm_zone(src.epsg)
    zd, sd = _utm_zone(dst.epsg)
    if src.epsg == 4326 and zd:                         # the common cases need no library
        return _utm_fwd(xs, ys, zd, sd)
    if dst.epsg == 4326 and zs:
        return _utm_inv(xs, ys, zs, ss)
    if zs and zd:
        lon, lat = _utm_inv(xs, ys, zs, ss)
        return _utm_fwd(lon, lat, zd, sd)
    if backend() == "rasterio":
        from rasterio.warp import transform
        x2, y2 = transform(_rio_crs(src), _rio_crs(dst), list(xs), list(ys))
        return np.asarray(x2, dtype=float), np.asarray(y2, dtype=float)
    raise RuntimeError("cannot transform coordinates from %r to %r without rasterio (install it, or give the data in WGS84 / UTM)" % (src, dst))


def utm_wkt(epsg):
    z, south = _utm_zone(epsg)
    if not z:
        return ""
    return ('PROJCS["WGS 84 / UTM zone %d%s",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],'
            'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",%d],'
            'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",%d],UNIT["metre",1],'
            'AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","%d"]]' % (z, "S" if south else "N", z * 6 - 183, 10000000 if south else 0, epsg))


# =============================================================================
# rasters
# =============================================================================
def _rio_crs(crs):
    from rasterio.crs import CRS
    if crs.epsg:
        return CRS.from_epsg(crs.epsg)
    return CRS.from_wkt(crs.wkt)


def _crs_from_rio(c):
    if c is None:
        return Crs()
    try:
        epsg = c.to_epsg()
    except Exception:
        epsg = None
    return Crs(wkt=c.to_wkt(), epsg=epsg, geographic=bool(c.is_geographic))


def _basic_read(path):
    import tifffile
    if path.lower().endswith(".vrt"):                   # RAS Mapper writes a .vrt next to one .tif: follow it
        txt = open(path, encoding="utf-8", errors="ignore").read()
        srcs = re.findall(r"<SourceFilename[^>]*>([^<]+)</SourceFilename>", txt)
        if len(set(srcs)) != 1:
            raise RuntimeError("%s mosaics %d files: install rasterio to read it" % (os.path.basename(path), len(set(srcs))))
        path = os.path.join(os.path.dirname(path), srcs[0])
    with tifffile.TiffFile(path) as t:
        p = t.pages[0]
        try:
            arr = p.asarray()
        except ValueError as e:                          # e.g. LZW without the imagecodecs package: Pillow can usually decode it
            try:
                from PIL import Image
                Image.MAX_IMAGE_PIXELS = None
                with Image.open(path) as im:
                    arr = np.array(im)
            except Exception:
                raise RuntimeError("%s: %s - install rasterio (or imagecodecs) to read this compression" % (os.path.basename(path), e))
        if arr.ndim == 3:
            arr = arr[..., 0] if arr.shape[-1] <= 4 else arr[0]
        tg = p.tags
        if "ModelPixelScaleTag" not in tg or "ModelTiepointTag" not in tg:
            raise RuntimeError("%s has no simple georeferencing (pixel scale + tie point): install rasterio to read it" % os.path.basename(path))
        sc, tp = tg["ModelPixelScaleTag"].value, tg["ModelTiepointTag"].value
        nd = None
        if "GDAL_NODATA" in tg:
            try:
                nd = float(str(tg["GDAL_NODATA"].value).strip().strip("\x00"))
            except ValueError:
                nd = None
        keys = tuple(tg["GeoKeyDirectoryTag"].value) if "GeoKeyDirectoryTag" in tg else None
    epsg, geographic, point = None, False, False
    if keys:
        for k in range(4, len(keys) - 3, 4):
            kid, loc, _cnt, val = keys[k:k + 4]
            if loc != 0:
                continue
            if kid == 1024:
                geographic = (val == 2)
            elif kid == 1025:
                point = (val == 2)
            elif kid == 3072 and 0 < val < 32767:
                epsg = val
            elif kid == 2048 and epsg is None and 0 < val < 32767:
                epsg = val
    x0, y0 = tp[3] - tp[0] * sc[0], tp[4] + tp[1] * sc[1]
    if point:
        x0, y0 = x0 - 0.5 * sc[0], y0 + 0.5 * sc[1]
    crs = Crs(wkt=utm_wkt(epsg) if epsg else "", epsg=epsg, geographic=geographic, geokeys=keys)
    return arr, (x0, sc[0], 0.0, y0, 0.0, -sc[1]), crs, nd


def read_raster_raw(path):
    """-> (array as stored, geotransform, Crs, declared nodata or None). Band 1 only."""
    if not os.path.exists(path):
        raise FileNotFoundError("raster not found: %s" % path)
    if backend() == "rasterio":
        import rasterio
        with rasterio.open(path) as ds:
            t = ds.transform
            return ds.read(1), (t.c, t.a, t.b, t.f, t.d, t.e), _crs_from_rio(ds.crs), ds.nodata
    return _basic_read(path)


def read_raster(path):
    """-> (float64 array with NaN for nodata, geotransform, Crs)."""
    arr, gt, crs, nd = read_raster_raw(path)
    a = arr.astype(np.float64)
    if nd is not None and nd == nd:
        a[np.isclose(a, nd)] = np.nan
    a[~np.isfinite(a)] = np.nan
    return a, gt, crs


def write_raster(path, arr, gt, crs, dtype="float32", nodata=None):
    a = np.asarray(arr)
    if nodata is not None and a.dtype.kind == "f":
        a = np.where(np.isfinite(a), a, nodata)
    a = a.astype(dtype)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if backend() == "rasterio":
        import rasterio
        from rasterio.transform import Affine
        kw = dict(driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, dtype=dtype, transform=Affine.from_gdal(*gt), compress="deflate")
        if crs is not None and (crs.epsg or crs.wkt):
            kw["crs"] = _rio_crs(crs)
        if nodata is not None:
            kw["nodata"] = nodata
        with rasterio.open(path, "w", **kw) as dst:
            dst.write(a, 1)
        return path
    import tifffile
    ex = [(33550, 12, 3, (abs(gt[1]), abs(gt[5]), 0.0), True), (33922, 12, 6, (0.0, 0.0, 0.0, gt[0], gt[3], 0.0), True)]
    keys = None
    if crs is not None:
        if crs.epsg:
            keys = ((1, 1, 0, 3, 1024, 0, 1, 2, 1025, 0, 1, 1, 2048, 0, 1, crs.epsg) if crs.geographic else
                    (1, 1, 0, 3, 1024, 0, 1, 1, 1025, 0, 1, 1, 3072, 0, 1, crs.epsg))
        elif crs.geokeys:
            keys = tuple(crs.geokeys)
    if keys:
        ex.append((34735, "H", len(keys), keys, True))
    if nodata is not None:
        ex.append((42113, "s", 0, repr(float(nodata)) if dtype.startswith("float") else str(int(nodata)), True))
    tifffile.imwrite(path, a, extratags=ex, compression="deflate")
    return path


def _sample(a, gt_src, X, Y, alg):
    cf, rf = (X - gt_src[0]) / gt_src[1] - 0.5, (Y - gt_src[3]) / gt_src[5] - 0.5
    H, W = a.shape
    if alg == "nearest":
        r, c = np.rint(rf).astype(np.int64), np.rint(cf).astype(np.int64)
        ok = (r >= 0) & (r < H) & (c >= 0) & (c < W)
        out = np.full(X.shape, np.nan)
        out[ok] = a[r[ok], c[ok]]
        return out
    r0, c0 = np.floor(rf).astype(np.int64), np.floor(cf).astype(np.int64)
    fr, fc = rf - r0, cf - c0
    ok = (rf >= -0.5) & (rf <= H - 0.5) & (cf >= -0.5) & (cf <= W - 0.5)
    r0c, r1c, c0c, c1c = np.clip(r0, 0, H - 1), np.clip(r0 + 1, 0, H - 1), np.clip(c0, 0, W - 1), np.clip(c0 + 1, 0, W - 1)
    v = np.stack([a[r0c, c0c], a[r0c, c1c], a[r1c, c0c], a[r1c, c1c]])
    w = np.stack([(1 - fr) * (1 - fc), (1 - fr) * fc, fr * (1 - fc), fr * fc])
    w = np.where(np.isfinite(v), w, 0.0)
    ws = w.sum(axis=0)
    out = np.where(ws > 0, (np.nan_to_num(v) * w).sum(axis=0) / np.where(ws > 0, ws, 1.0), np.nan)
    out[~ok] = np.nan
    return out


def raster_on_grid(path, shape, gt, crs, alg="bilinear"):
    """A raster file resampled onto the grid (shape, gt, crs): float64 with NaN outside / nodata. alg: bilinear | nearest."""
    H, W = shape
    if backend() == "rasterio":
        import rasterio
        from rasterio.transform import Affine
        from rasterio.warp import reproject, Resampling
        with rasterio.open(path) as ds:
            if ds.crs is None:
                raise RuntimeError("%s has no CRS" % os.path.basename(path))
            dst = np.full((H, W), np.nan, dtype=np.float64)
            reproject(rasterio.band(ds, 1), dst, src_transform=ds.transform, src_crs=ds.crs, src_nodata=ds.nodata, dst_transform=Affine.from_gdal(*gt),
                      dst_crs=_rio_crs(crs), dst_nodata=np.nan, resampling=Resampling.bilinear if alg == "bilinear" else Resampling.nearest)
        return dst
    a, g, c = read_raster(path)
    if c.epsg and crs.epsg and c.epsg != crs.epsg:
        raise RuntimeError("%s is in %s but the DEM is in %s: reproject it, or install rasterio so fanrisk can" % (os.path.basename(path), c.authid(), crs.authid()))
    cols, rows = np.arange(W), np.arange(H)
    X = np.broadcast_to(gt[0] + (cols + 0.5) * gt[1], (H, W))
    Y = np.broadcast_to((gt[3] + (rows + 0.5) * gt[5])[:, None], (H, W))
    return _sample(a, g, X, Y, alg)


def load_dem(path, reproject_geographic=True, log=print):
    """-> (array as stored, geotransform, Crs, declared nodata). A DEM in degrees is warped to UTM (needs rasterio)."""
    arr, gt, crs, nd = read_raster_raw(path)
    if abs(gt[2]) > 1e-9 or abs(gt[4]) > 1e-9 or gt[5] >= 0:
        raise RuntimeError("Only north-up, non-rotated rasters are supported (warp the DEM first).")
    looks_geo = crs.geographic or (not crs.epsg and not crs.wkt and abs(gt[1]) < 0.01 and abs(gt[0]) <= 180 and abs(gt[3]) <= 90)
    if not looks_geo:
        if not crs.epsg and not crs.wkt:
            log("WARNING: the DEM has no CRS; it is taken as projected, in metres.")
        return arr, gt, crs, nd
    if not reproject_geographic:
        raise RuntimeError("The DEM is in a geographic CRS (degrees); reproject it to UTM or set REPROJECT_GEOGRAPHIC = true.")
    if backend() != "rasterio":
        raise RuntimeError("The DEM is in degrees. Reproject it to UTM first, or install rasterio so fanrisk can do it (pip install rasterio).")
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling
    with rasterio.open(path) as ds:
        b = ds.bounds
        lon, lat = 0.5 * (b.left + b.right), 0.5 * (b.bottom + b.top)
        epsg = utm_epsg_for(lon, lat)
        res_m = 0.5 * (abs(gt[1]) * 111320.0 * math.cos(math.radians(lat)) + abs(gt[5]) * 110540.0)
        log("DEM is in degrees - warping to EPSG:%d at %.1f m ..." % (epsg, res_m))
        dst_crs = rasterio.crs.CRS.from_epsg(epsg)
        tr, w, h = calculate_default_transform(ds.crs or "EPSG:4326", dst_crs, ds.width, ds.height, *b, resolution=res_m)
        dst = np.full((h, w), NODATA, dtype=np.float32)
        reproject(rasterio.band(ds, 1), dst, src_transform=ds.transform, src_crs=ds.crs or "EPSG:4326", src_nodata=ds.nodata, dst_transform=tr,
                  dst_crs=dst_crs, dst_nodata=NODATA, resampling=Resampling.bilinear)
    return dst, (tr.c, tr.a, tr.b, tr.f, tr.d, tr.e), Crs(wkt=dst_crs.to_wkt(), epsg=epsg), NODATA


# =============================================================================
# mask <-> geometry
# =============================================================================
def _trace_loops(mask):
    """Boundary loops of a boolean mask as lists of (col, row) grid corners; the mask is on the right of the walking direction
    (so in row-down coordinates outer loops have positive shoelace area, holes negative). Diagonal neighbours are separate parts."""
    m = np.pad(np.asarray(mask, dtype=bool), 1)
    r, c = np.nonzero(m[1:-1, 1:-1])
    nxt = {}

    def add(sel, a, b):
        for rr, cc in zip(r[sel].tolist(), c[sel].tolist()):
            nxt.setdefault((a[0] + cc, a[1] + rr), []).append((b[0] + cc, b[1] + rr))
    add(~m[r, c + 1], (0, 0), (1, 0))                  # top edge, walked left -> right
    add(~m[r + 1, c + 2], (1, 0), (1, 1))              # right edge, downwards
    add(~m[r + 2, c + 1], (1, 1), (0, 1))              # bottom edge, right -> left
    add(~m[r + 1, c], (0, 1), (0, 0))                  # left edge, upwards
    loops = []
    for start in list(nxt.keys()):
        while nxt.get(start):
            a = start
            b = nxt[a].pop()
            loop = [a]
            while b != start:
                loop.append(b)
                outs = nxt[b]
                if len(outs) == 1:
                    a, b = b, outs.pop()
                else:                                   # pinch point (two cells touching at a corner): turn right = stay with this part
                    dx, dy = b[0] - a[0], b[1] - a[1]
                    want = (b[0] - dy, b[1] + dx)
                    k = outs.index(want) if want in outs else 0
                    a, b = b, outs.pop(k)
            keep = [loop[i] for i in range(len(loop))
                    if (loop[i][0] - loop[i - 1][0]) * (loop[(i + 1) % len(loop)][1] - loop[i][1]) != (loop[i][1] - loop[i - 1][1]) * (loop[(i + 1) % len(loop)][0] - loop[i][0])]
            loops.append(keep)
    return loops


def _area2(loop):
    a = np.asarray(loop, dtype=float)
    x, y = a[:, 0], a[:, 1]
    return float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def _in_ring(px, py, ring):
    inside = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        if (y1 > py) != (y2 > py) and px < x1 + (py - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def polygonize(mask, gt):
    """Boolean mask -> GeoJSON geometry (Polygon / MultiPolygon) in world coordinates, or None when the mask is empty."""
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return None
    if backend() == "rasterio":
        from rasterio.features import shapes
        from rasterio.transform import Affine
        polys = [g["coordinates"] for g, v in shapes(mask.astype(np.uint8), mask=mask, connectivity=8, transform=Affine.from_gdal(*gt)) if v == 1]
        polys = [[[[float(x), float(y)] for x, y in ring] for ring in p] for p in polys]
    else:
        loops = _trace_loops(mask)
        outers = sorted([l for l in loops if _area2(l) > 0], key=lambda l: abs(_area2(l)))
        holes = [l for l in loops if _area2(l) < 0]
        parts = [[o] for o in outers]
        for h in holes:
            (ax, ay), (bx, by) = h[0], h[1]             # the mask lies on the right of every edge: step half a cell in from the first one
            L = float(math.hypot(bx - ax, by - ay))
            ux, uy = (bx - ax) / L, (by - ay) / L
            px, py = ax + 0.5 * ux - 0.5 * uy, ay + 0.5 * uy + 0.5 * ux
            for part in parts:                          # smallest outer first
                if _in_ring(px, py, part[0]):
                    part.append(h)
                    break
        w = lambda loop: [[gt[0] + x * gt[1], gt[3] + y * gt[5]] for x, y in loop + [loop[0]]]
        polys = [[w(l) for l in part] for part in parts]
    if not polys:
        return None
    return {"type": "Polygon", "coordinates": polys[0]} if len(polys) == 1 else {"type": "MultiPolygon", "coordinates": polys}


def smooth_geometry(geom, iterations=2, offset=0.25):
    """Chaikin corner cutting of every polygon ring (cosmetic: takes the staircase out of raster outlines)."""
    def ring_(r):
        p = np.asarray(r[:-1] if (len(r) > 1 and r[0] == r[-1]) else r, dtype=float)
        if len(p) < 4:
            return r
        for _ in range(iterations):
            q = np.roll(p, -1, axis=0)
            p = np.column_stack([(1 - offset) * p + offset * q, offset * p + (1 - offset) * q]).reshape(-1, 2)
        p = np.round(p, 3)
        return p.tolist() + [p[0].tolist()]
    if geom is None:
        return None
    if geom["type"] == "Polygon":
        return {"type": "Polygon", "coordinates": [ring_(r) for r in geom["coordinates"]]}
    if geom["type"] == "MultiPolygon":
        return {"type": "MultiPolygon", "coordinates": [[ring_(r) for r in poly] for poly in geom["coordinates"]]}
    return geom


def _poly_mask(rings, shape, gt):
    H, W = shape
    out = np.zeros(shape, dtype=bool)
    for ring in rings:                                   # even-odd: holes flip back
        a = np.asarray(ring, dtype=float)[:, :2]
        if len(a) < 4:
            continue
        c0 = max(0, int(math.floor((a[:, 0].min() - gt[0]) / gt[1])))
        c1 = min(W, int(math.ceil((a[:, 0].max() - gt[0]) / gt[1])) + 1)
        r0 = max(0, int(math.floor((a[:, 1].max() - gt[3]) / gt[5])))
        r1 = min(H, int(math.ceil((a[:, 1].min() - gt[3]) / gt[5])) + 1)
        if c1 <= c0 or r1 <= r0:
            continue
        X = (gt[0] + (np.arange(c0, c1) + 0.5) * gt[1])[None, :]
        Y = (gt[3] + (np.arange(r0, r1) + 0.5) * gt[5])[:, None]
        ins = np.zeros((r1 - r0, c1 - c0), dtype=bool)
        for (x1, y1), (x2, y2) in zip(a[:-1], a[1:]):
            if y1 == y2:
                continue
            cond = (y1 > Y) != (y2 > Y)
            xs = x1 + (Y - y1) * (x2 - x1) / (y2 - y1)
            ins ^= cond & (X < xs)
        out[r0:r1, c0:c1] ^= ins
    return out


def rasterize(geoms, shape, gt):
    """GeoJSON geometries (grid CRS) -> boolean array: polygons by cell centre, lines and points by the cells they touch."""
    geoms = [g for g in geoms if g]
    out = np.zeros(shape, dtype=bool)
    if not geoms:
        return out
    if backend() == "rasterio":
        from rasterio.features import rasterize as rz
        from rasterio.transform import Affine
        return rz([(g, 1) for g in geoms], out_shape=shape, transform=Affine.from_gdal(*gt), fill=0, dtype="uint8").astype(bool)
    H, W = shape

    def pts_(coords):
        a = np.asarray(coords, dtype=float).reshape(-1, np.asarray(coords).shape[-1])[:, :2]
        c, r = np.floor((a[:, 0] - gt[0]) / gt[1]).astype(int), np.floor((a[:, 1] - gt[3]) / gt[5]).astype(int)
        ok = (r >= 0) & (r < H) & (c >= 0) & (c < W)
        out[r[ok], c[ok]] = True

    def line_(coords):
        a = np.asarray(coords, dtype=float)[:, :2]
        for p, q in zip(a[:-1], a[1:]):
            n = max(2, int(math.hypot(*(q - p)) / (0.5 * abs(gt[1]))) + 1)
            pts_(np.column_stack([np.linspace(p[0], q[0], n), np.linspace(p[1], q[1], n)]))
    for g in geoms:
        t, co = g["type"], g["coordinates"]
        if t == "Polygon":
            out |= _poly_mask(co, shape, gt)
        elif t == "MultiPolygon":
            for poly in co:
                out |= _poly_mask(poly, shape, gt)
        elif t == "LineString":
            line_(co)
        elif t == "MultiLineString":
            for l in co:
                line_(l)
        elif t == "Point":
            pts_([co])
        elif t == "MultiPoint":
            pts_(co)
    return out


# =============================================================================
# vectors
# =============================================================================
def _map_coords(obj, fn):
    if isinstance(obj, (list, tuple)) and obj and isinstance(obj[0], (int, float)):
        return fn([obj])[0]
    if isinstance(obj, (list, tuple)) and obj and isinstance(obj[0], (list, tuple)) and obj[0] and isinstance(obj[0][0], (int, float)):
        return fn(obj)
    return [_map_coords(o, fn) for o in obj]


def _reproject_geom(g, src, dst):
    def fn(pts):
        a = np.asarray(pts, dtype=float)
        x, y = transform_xy(a[:, 0], a[:, 1], src, dst)
        return [[float(u), float(v)] for u, v in zip(x, y)]
    return {"type": g["type"], "coordinates": _map_coords(g["coordinates"], fn)}


def _geojson_crs(doc):
    name = (((doc.get("crs") or {}).get("properties") or {}).get("name") or "")
    m = re.search(r"EPSG:{1,2}(\d+)", name)
    if m:
        e = int(m.group(1))
        return Crs(epsg=e, geographic=(e == 4326), wkt=utm_wkt(e))
    return WGS84


def read_vector(path, dst_crs=None):
    """A vector file -> [dict(geometry=GeoJSON geometry in dst_crs, properties=dict)]. Without geopandas only GeoJSON is read."""
    if not os.path.exists(path):
        raise FileNotFoundError("vector file not found: %s" % path)
    if backend() == "rasterio" and have("geopandas"):
        import geopandas as gpd
        from shapely.geometry import mapping
        gdf = gpd.read_file(path)
        if dst_crs is not None and gdf.crs is not None:
            gdf = gdf.to_crs(dst_crs.authid() or dst_crs.wkt)
        cols = [c for c in gdf.columns if c != gdf.geometry.name]
        out = []
        for _, row in gdf.iterrows():
            g = row[gdf.geometry.name]
            if g is None or g.is_empty:
                continue
            out.append(dict(geometry=json.loads(json.dumps(mapping(g))), properties={c: row[c] for c in cols}))
        return out
    if not path.lower().endswith((".geojson", ".json")):
        raise RuntimeError("%s: only GeoJSON can be read without geopandas (pip install geopandas), or save the layer as GeoJSON" % os.path.basename(path))
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    src = _geojson_crs(doc)
    feats = doc["features"] if doc.get("type") == "FeatureCollection" else [doc]
    out = []
    for f in feats:
        g = f.get("geometry")
        if not g or not g.get("coordinates"):
            continue
        if g["type"] == "GeometryCollection":
            continue
        if dst_crs is not None and not src.same(dst_crs):
            g = _reproject_geom(g, src, dst_crs)
        out.append(dict(geometry=g, properties=f.get("properties") or {}))
    return out


def _clean(v):
    if v is None:
        return None
    if isinstance(v, (np.floating, float)):
        return None if v != v else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def write_geojson(path, layer, crs):
    feats = [dict(type="Feature", properties={n: _clean(v) for (n, _t), v in zip(layer["fields"], attrs)}, geometry=g) for g, attrs in layer["rows"] if g]
    doc = dict(type="FeatureCollection", name=layer["name"], features=feats)
    if crs is not None and crs.epsg:
        doc["crs"] = dict(type="name", properties=dict(name="urn:ogc:def:crs:EPSG::%d" % crs.epsg))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh)
    return path


def write_layers(out_dir, gpkg_name, layers, crs, log=print):
    """layers: [dict(name, fields=[(name, 'integer'|'double'|'string')], rows=[(geometry, [values])])].
    Written into one GeoPackage when geopandas is installed, otherwise one GeoJSON per layer in <out_dir>/vectors. Returns the paths written."""
    layers = [l for l in layers if l["rows"]]
    if backend() == "rasterio" and have("geopandas"):
        try:
            import geopandas as gpd
            from shapely.geometry import shape
            gpkg = os.path.join(out_dir, gpkg_name)
            if os.path.exists(gpkg):
                os.remove(gpkg)
            cast = dict(integer=lambda v: None if v is None else int(v), double=lambda v: None if v is None else float(v), string=lambda v: None if v is None else str(v))
            for l in layers:
                rows = [(g, a) for g, a in l["rows"] if g]
                data = {n: [cast[t](_clean(a[i])) for _g, a in rows] for i, (n, t) in enumerate(l["fields"])}
                gdf = gpd.GeoDataFrame(data, geometry=[shape(g) for g, _a in rows], crs=(crs.authid() or crs.wkt or None) if crs else None)
                gdf.to_file(gpkg, layer=l["name"], driver="GPKG")
            return [gpkg]
        except Exception as e:
            log("WARNING: the GeoPackage could not be written (%s: %s) - writing GeoJSON files instead." % (type(e).__name__, e))
    vdir = os.path.join(out_dir, "vectors")
    return [write_geojson(os.path.join(vdir, l["name"] + ".geojson"), l, crs) for l in layers]


def write_shapefile(path, geom, crs, label):
    """One geometry with a 'name' field as an ESRI shapefile (RAS Mapper imports these). Needs geopandas; returns False without it."""
    if not (backend() == "rasterio" and have("geopandas")):
        return False
    import geopandas as gpd
    from shapely.geometry import shape
    gdf = gpd.GeoDataFrame(dict(name=[label]), geometry=[shape(geom)], crs=(crs.authid() or crs.wkt or None) if crs else None)
    gdf.to_file(path)
    return True


def geom_area(g):
    """Planar area of a polygon geometry (outer rings minus holes)."""
    def ring_area(r):
        a = np.asarray(r, dtype=float)
        return 0.5 * abs(float(np.dot(a[:-1, 0], a[1:, 1]) - np.dot(a[1:, 0], a[:-1, 1])))
    polys = [g["coordinates"]] if g["type"] == "Polygon" else (g["coordinates"] if g["type"] == "MultiPolygon" else [])
    return sum(ring_area(p[0]) - sum(ring_area(h) for h in p[1:]) for p in polys if p)


def first_line(features):
    """First LineString of read_vector() output as [(x, y), ...]."""
    for f in features:
        g = f["geometry"]
        co = g["coordinates"] if g["type"] == "LineString" else (g["coordinates"][0] if g["type"] == "MultiLineString" and g["coordinates"] else None)
        if co and len(co) >= 2:
            return [(float(p[0]), float(p[1])) for p in co]
    return None
