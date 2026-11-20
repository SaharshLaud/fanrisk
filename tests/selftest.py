# -*- coding: utf-8 -*-
"""fanrisk self-test:  python tests/selftest.py [--dem some_dem.tif] [--work folder]
Checks the file layer with whichever backend is installed and, with --dem, runs the whole pipeline on that DEM in a scratch folder."""
import argparse
import json
import os
import sys
import tempfile
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fanrisk import gio, config, engine, hydrology, __version__  # noqa: E402

FAILS = []


def check(name, fn):
    try:
        r = fn()
        print("  [PASS] %s%s" % (name, "" if r in (None, True) else " - %s" % (r,)))
    except Exception as e:
        FAILS.append(name)
        print("  [FAIL] %s: %s: %s" % (name, type(e).__name__, e))
        traceback.print_exc()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem")
    ap.add_argument("--work")
    a = ap.parse_args()
    work = a.work or tempfile.mkdtemp(prefix="fanrisk_test_")
    os.makedirs(work, exist_ok=True)
    print("fanrisk %s self-test | python %s | backend %s | geopandas %s | work folder %s" % (__version__, sys.version.split()[0], gio.backend(), gio.have("geopandas"), work))
    gt, crs = (500000.0, 30.0, 0.0, 3800000.0, 0.0, -30.0), gio.Crs(epsg=32643)
    rng = np.random.default_rng(1)

    def t_utm():
        x, y = gio.transform_xy([-79.387139], [43.642566], gio.WGS84, gio.Crs(epsg=32617))      # a published reference point (630084 E, 4833438 N)
        assert abs(x[0] - 630084) < 2 and abs(y[0] - 4833438) < 2, (x, y)
        lo, la = gio.transform_xy(x, y, gio.Crs(epsg=32617), gio.WGS84)
        assert abs(lo[0] + 79.387139) < 1e-7 and abs(la[0] - 43.642566) < 1e-7
    check("UTM <-> longitude / latitude", t_utm)

    def t_raster():
        z = rng.random((40, 50)).astype("float32")
        z[0, 0] = np.nan
        p = gio.write_raster(os.path.join(work, "t.tif"), z, gt, crs, "float32", -9999.0)
        b, g2, c2 = gio.read_raster(p)
        assert np.allclose(z, b, equal_nan=True) and np.allclose(g2, gt) and c2.epsg == 32643, (g2, c2)
        r = gio.raster_on_grid(p, (40, 50), gt, crs)
        assert np.allclose(r[2:-2, 2:-2], z[2:-2, 2:-2], atol=1e-5)
        h = gio.raster_on_grid(p, (20, 25), (gt[0], 60.0, 0, gt[3], 0, -60.0), crs, "nearest")
        assert h.shape == (20, 25) and np.isfinite(h).sum() > 400
    check("GeoTIFF write / read / resample", t_raster)

    def t_poly():
        for k in range(20):
            m = rng.random((30, 30)) > 0.45
            g = gio.polygonize(m, gt)
            assert np.array_equal(gio.rasterize([g], m.shape, gt), m), k
            assert abs(gio.geom_area(g) - m.sum() * 900.0) < 1e-3
        assert gio.smooth_geometry(g)["type"] == g["type"]
    check("mask -> polygons -> mask (20 random masks)", t_poly)

    def t_vec():
        m = np.zeros((30, 30), bool)
        m[5:20, 5:25] = True
        g = gio.polygonize(m, gt)
        layers = [dict(name="polys", fields=[("id", "integer"), ("name", "string"), ("v", "double")], rows=[(g, [1, "a", 2.5])]),
                  dict(name="pts", fields=[("id", "integer")], rows=[({"type": "Point", "coordinates": [500100.0, 3799900.0]}, [7])])]
        paths = gio.write_layers(work, "t.gpkg", layers, crs)
        gj = gio.write_geojson(os.path.join(work, "polys.geojson"), layers[0], crs)
        f = gio.read_vector(gj, crs)
        assert len(f) == 1 and np.array_equal(gio.rasterize([f[0]["geometry"]], m.shape, gt), m)
        f4 = gio.read_vector(gj, gio.Crs(epsg=32644))                  # reprojection of a vector
        assert abs(f4[0]["geometry"]["coordinates"][0][0][0] - 500000) > 1000
        ok = gio.write_shapefile(os.path.join(work, "t.shp"), g, crs, "x")
        return "%s; shapefile %s" % (os.path.basename(paths[0]), "written" if ok else "skipped (no geopandas)")
    check("vector layers write / read", t_vec)

    def t_geo():
        if gio.backend() != "rasterio":
            return "skipped (needs rasterio)"
        z = (1000 + 50 * rng.random((60, 60))).astype("float32")
        p = gio.write_raster(os.path.join(work, "geo.tif"), z, (77.5, 0.0003, 0, 34.2, 0, -0.0003), gio.WGS84, "float32", -9999.0)
        arr, g2, c2, nd = gio.load_dem(p, True, lambda m: None)
        assert c2.epsg == 32643 and 20 < g2[1] < 40 and np.sum(arr > 900) > 2000, (c2, g2)
        return "DEM in degrees warped to %s at %.1f m" % (c2.authid(), g2[1])
    check("DEM in degrees -> UTM", t_geo)

    def t_cn():
        q, cn = engine.runoff_depth_mm(60, dict(RUNOFF_METHOD="scs_cn", CURVE_NUMBER=80, CN_AMC=3, CN_IA_RATIO=0.2))
        assert abs(q - 36.2) < 0.3 and abs(cn - 90.2) < 0.2, (q, cn)
        t = hydrology.cn_from_landcover(np.array([[10.0, 60.0, 50.0, np.nan]]), np.array([[1.0, 4.0, 12.0, 2.0]]))
        assert t[0, 0] == 36 and t[0, 1] == 94 and t[0, 2] == 92 and np.isnan(t[0, 3])
    check("SCS curve-number runoff and land-cover table", t_cn)

    def t_cfg():
        for n in config.available("profiles") + config.available("cases"):
            pass
        d = os.path.join(work, "cfg")
        os.makedirs(d, exist_ok=True)
        gio.write_raster(os.path.join(d, "dem.tif"), rng.random((10, 10)).astype("float32"), gt, crs, "float32", -9999.0)
        open(os.path.join(d, "site.toml"), "w").write('[site]\ndem = "dem.tif"\nvalidation_case = "leh_2010"\n[storm]\nrainfall_mm = 50\nduration_min = 20\n[hecras]\ninflow = "observed"\n')
        s = config.load_site(d)
        assert s.P["RAS_OBSERVED_HYDROGRAPH"]["peak_m3s"] == 1070 and s.case["probes"][0]["depth_m"] == 5.2
        open(os.path.join(d, "site.toml"), "w").write('[site]\ndem = "dem.tif"\n[storm]\nrainfal_mm = 50\n')
        try:
            config.load_site(d)
            raise AssertionError("a misspelt key was accepted")
        except config.ConfigError as e:
            assert "rainfall_mm" in str(e)
    check("site file, profiles, validation case, typo detection", t_cfg)

    if a.dem:
        def t_run():
            from fanrisk import pipeline
            d = os.path.join(work, "site")
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, "site.toml"), "w").write('[site]\nname = "self-test"\ndem = %s\nregion_profile = "arid_mountain"\n[stages]\nrisk = true\nhecras_handoff = true\n'
                                                          '[storm]\nrainfall_mm = 50\nduration_min = 20\n[targets]\nmode = "top"\n' % json.dumps(os.path.abspath(a.dem)))
            r = pipeline.run(config.load_site(d), echo=False)
            assert not r["errors"], r["errors"]
            return "%d fans in %d complexes; handoff %s; outputs in %s" % (r["n_fans"], r["n_complexes"], "written" if r["handoff"] else "not written (no fan found)", r["output_dir"])
        check("whole pipeline on %s" % os.path.basename(a.dem), t_run)
    print("SELF-TEST %s" % ("FAILED: " + ", ".join(FAILS) if FAILS else "PASSED"))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
