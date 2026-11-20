# -*- coding: utf-8 -*-
"""Download a DEM for any area: Copernicus GLO-30 (30 m, free, no account) from the AWS open-data bucket, mosaicked and clipped.

    python tools/get_dem.py <south> <west> <north> <east> <out.tif>
    e.g. python tools/get_dem.py 30.93 78.70 31.12 78.92 data/dharali/dem.tif

The result is in degrees (EPSG:4326); 'fanrisk run' warps it to UTM by itself. Needs rasterio.
Licence: Copernicus DEM (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the
European Union and ESA; free for any use with this attribution. https://registry.opendata.aws/copernicus-dem/
"""
import math
import os
import sys
import tempfile
import urllib.error
import urllib.request

BASE = "https://copernicus-dem-30m.s3.amazonaws.com"


def tile_name(lat0, lon0):
    ns, ew = ("N" if lat0 >= 0 else "S"), ("E" if lon0 >= 0 else "W")
    return "Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM" % (ns, abs(lat0), ew, abs(lon0))


def main(a):
    if len(a) < 6:
        raise SystemExit(__doc__)
    s, w, n, e = map(float, a[1:5])
    out = a[5]
    if not (s < n and w < e):
        raise SystemExit("give the box as south west north east (south < north, west < east)")
    import rasterio
    from rasterio.merge import merge
    tmp = tempfile.mkdtemp(prefix="fanrisk_dem_")
    paths = []
    for lat0 in range(math.floor(s), math.floor(n - 1e-9) + 1):
        for lon0 in range(math.floor(w), math.floor(e - 1e-9) + 1):
            t = tile_name(lat0, lon0)
            url = "%s/%s/%s.tif" % (BASE, t, t)
            dst = os.path.join(tmp, t + ".tif")
            print("downloading %s ..." % t)
            try:
                urllib.request.urlretrieve(url, dst)
            except urllib.error.HTTPError as ex:
                print("  not available (%s) - sea or a tile not yet public; skipped" % ex.code)
                continue
            print("  %.1f MB" % (os.path.getsize(dst) / 1e6))
            paths.append(dst)
    if not paths:
        raise SystemExit("No tiles downloaded for this box.")
    srcs = [rasterio.open(p) for p in paths]
    arr, tr = merge(srcs, bounds=(w, s, e, n), nodata=-32767.0)
    prof = srcs[0].profile.copy()
    for x in srcs:
        x.close()
    prof.update(driver="GTiff", height=arr.shape[1], width=arr.shape[2], transform=tr, count=1, nodata=-32767.0, compress="deflate",
                tiled=False, blockxsize=None, blockysize=None)
    prof.pop("blockxsize", None)
    prof.pop("blockysize", None)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with rasterio.open(out, "w", **prof) as d:
        d.write(arr[0], 1)
    print("written %s: %d x %d cells, box %.3f-%.3f N, %.3f-%.3f E" % (out, arr.shape[2], arr.shape[1], s, n, w, e))


if __name__ == "__main__":
    main(sys.argv)
