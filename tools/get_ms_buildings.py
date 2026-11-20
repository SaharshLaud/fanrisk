# -*- coding: utf-8 -*-
"""Download Microsoft Global ML Building Footprints for one area (OSM misses most houses in Ladakh).

    python tools/get_ms_buildings.py <south> <west> <north> <east> [out folder]
    e.g. python tools/get_ms_buildings.py 34.086 77.539 34.134 77.630 data

Finds the zoom-9 quadkey tile(s) covering the box in Microsoft's dataset-links.csv and saves each as
<out>/ms_buildings_<quadkey>.csv.gz (gzipped GeoJSON lines, one building polygon per line, lon/lat). Licence: ODbL.
Source: https://github.com/microsoft/GlobalMLBuildingFootprints
"""
import math
import os
import sys
import urllib.request

LINKS = ["https://bfppub.blob.core.windows.net/%24web/2026-08-13/dataset-links.csv",
         "https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv"]


def quadkey(lat, lon, z=9):
    n = 2 ** z
    x = int((lon + 180.0) / 360.0 * n)
    s = math.sin(math.radians(lat))
    y = int((0.5 - math.log((1 + s) / (1 - s)) / (4 * math.pi)) * n)
    return "".join(str(((x >> i) & 1) + 2 * ((y >> i) & 1)) for i in range(z - 1, -1, -1))


def main(a):
    if len(a) < 5:
        raise SystemExit(__doc__)
    s, w, n, e = map(float, a[1:5])
    out = a[5] if len(a) > 5 else "data"
    os.makedirs(out, exist_ok=True)
    qks = sorted({quadkey(la, lo) for la in (s, n) for lo in (w, e)})
    print("tiles needed (zoom 9 quadkeys): %s" % ", ".join(qks))
    txt = None
    for u in LINKS:
        try:
            print("reading the index %s ..." % u)
            txt = urllib.request.urlopen(u, timeout=120).read().decode("utf-8", "ignore")
            break
        except Exception as ex:
            print("  failed: %s" % ex)
    if not txt:
        raise SystemExit("Could not read Microsoft's dataset-links.csv - check the internet connection.")
    rows = [r.split(",") for r in txt.splitlines()]
    head = [h.strip().lower() for h in rows[0]]
    iq, iu, il = head.index("quadkey"), head.index("url"), head.index("location")
    got = 0
    for r in rows[1:]:
        if len(r) > iu and r[iq].strip() in qks:
            dst = os.path.join(out, "ms_buildings_%s_%s.csv.gz" % (r[iq].strip(), r[il].strip().replace(" ", "")))
            if os.path.exists(dst) and os.path.getsize(dst) > 0:
                print("  have %s" % dst)
                got += 1
                continue
            print("  downloading %s (%s) ..." % (r[iq], r[il]))
            urllib.request.urlretrieve(r[iu].strip(), dst)
            print("  -> %s, %.1f MB" % (dst, os.path.getsize(dst) / 1e6))
            got += 1
    if not got:
        raise SystemExit("No tile for quadkeys %s in the index." % ", ".join(qks))
    print("done: %d file(s) in %s" % (got, out))


if __name__ == "__main__":
    main(sys.argv)
