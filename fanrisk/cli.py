# -*- coding: utf-8 -*-
"""Command line of fanrisk:  python -m fanrisk <command> [site.toml]"""
import argparse
import os
import shutil
import sys

from . import __version__, config


def _site(a):
    return config.load_site(a.site)


def cmd_init(a):
    os.makedirs(a.folder, exist_ok=True)
    dst = os.path.join(a.folder, "site.toml")
    if os.path.exists(dst):
        raise SystemExit("%s already exists - not overwritten" % dst)
    shutil.copy(os.path.join(config.HERE, "templates", "site.toml"), dst)
    print("Wrote %s\nPut your DEM next to it, edit the file, then:  python -m fanrisk check %s" % (dst, dst))


def cmd_profiles(a):
    for kind, label in (("profiles", "Region profiles"), ("cases", "Validation cases")):
        print(label + ":")
        for n in config.available(kind):
            d = config.load_toml(os.path.join(config.HERE, kind, n + ".toml"))
            h = d.get("profile") or d.get("case") or {}
            print("  %-16s %s" % (n, h.get("description") or h.get("name") or ""))
            if h.get("use_when"):
                print("  %-16s use when: %s" % ("", h["use_when"]))


def cmd_check(a):
    from . import gio
    s = _site(a)
    print("fanrisk %s | site '%s' | file backend: %s" % (__version__, s.name, gio.backend()))
    print("  [ OK ] site file reads; region profile '%s'%s" % (s.profile.get("name"), "; validation case '%s'" % s.case.get("name") if s.case else ""))
    arr, gt, crs, nd = gio.read_raster_raw(s.file("dem"))
    geo = crs.geographic or (not crs.epsg and abs(gt[1]) < 0.01)
    print("  [ OK ] DEM %s: %d x %d cells, %s%s" % (os.path.basename(s.file("dem")), arr.shape[1], arr.shape[0], crs.authid() or "no EPSG code",
                                                 " (degrees: will be reprojected to UTM)" if geo else ", cell %.1f m" % abs(gt[1])))
    if geo and gio.backend() != "rasterio":
        print("  [FAIL] a DEM in degrees needs rasterio (pip install rasterio geopandas)")
    for k, p in sorted(s.files.items()):
        if k != "dem":
            print("  [ OK ] %-20s %s" % (k, p))
    print("  stages: fans, %s%s" % ("risk" if s.stages["risk"] else "(no risk)", ", HEC-RAS handoff" if s.stages["hecras_handoff"] else ""))
    rep = config.settings_report(s)
    if rep:
        print("  settings that differ from the method defaults:")
        print("\n".join("  " + r for r in rep))
    for n in s.notes:
        print("  [NOTE] " + n)
    for mod, why in (("geopandas", "GeoPackage / shapefile output, non-GeoJSON vector input"), ("ras_commander", "running HEC-RAS"), ("matplotlib", "score figures")):
        print("  [%s] %s (%s)" % (" OK " if gio.have(mod) else "miss", mod, why))
    print("CHECK OK")


def cmd_run(a):
    from . import pipeline
    r = pipeline.run(_site(a))
    if r["errors"]:
        raise SystemExit(1)


def cmd_hecras(a):
    s = _site(a)
    from .ras import geometry
    if a.action == "fit":
        hd = s.handoff_dir()
        if not hd:
            raise SystemExit("No handoff yet: run 'fanrisk run' with [stages] hecras_handoff = true first.")
        if not a.template:                                # the bundled template: one 2D area with 'inflow' / 'outflow' lines (HEC-RAS 7.0.2)
            a.template = os.path.join(config.HERE, "templates", "hecras", "leh2d.prj")
            print("Using the bundled HEC-RAS template (%s). Give --template <your .prj> to use your own." % a.template)
        dst = a.out or os.path.join(s.out_dir, "hecras", "model")
        f = geometry.fit_project(a.template, hd, dst, base_plan=str(s.hecras["base_plan"]).zfill(2), spacing=a.spacing, version=str(s.hecras["version"]),
                                  mannings_n=s.hecras.get("mannings_n"))
        print('Now set in site.toml:  [hecras] project = "%s"   then:  python -m fanrisk hecras smoke %s' % (f["prj"].replace("\\", "/"), a.site))
        return
    if a.action == "extend-outflow":
        if not s.file("project"):
            raise SystemExit("[hecras] project is missing in site.toml")
        g = geometry.geometry_file(os.path.dirname(s.file("project")))
        shutil.copy2(g, g + ".before_extend")
        geometry.extend_outflow(g, float(a.metres))
        print("Geometry updated (previous version kept as %s). Re-run:  python -m fanrisk hecras run %s" % (os.path.basename(g) + ".before_extend", a.site))
        return
    from .ras import batch
    batch.configure(s, a.runs.split(",") if a.runs else None)
    batch.main(["fanrisk", a.action])


def cmd_score(a):
    s = _site(a)
    from .ras import readback, geometry
    hd = s.handoff_dir()
    if not hd:
        raise SystemExit("No handoff found in %s" % os.path.join(s.out_dir, "hecras"))
    import json
    cid = json.load(open(os.path.join(hd, "handoff_summary.json")))["complex"]
    bc = None
    if a.geometry:
        bc = list(geometry.read_bc_lines(a.geometry).values())
    readback.score_maps(hd, a.maps, a.out or os.path.join(a.maps, "..", "readback"), readback.settings_from_site(s),
                        ahp_risk_file=os.path.join(s.out_dir, "risk", "ahp_C%d_risk_abs.tif" % cid), footprint=s.file("flood_footprint"),
                        points=s.file("flood_points"), allow_snapshot=bool(s.score["allow_snapshot"]), bc_lines=bc)


def cmd_exposure(a):
    from . import exposure
    s = _site(a)
    osm = a.osm or os.path.join(os.path.dirname(os.path.abspath(s.file("dem"))), "map.osm")
    exposure.run(s, osm, maps_dir=a.maps or "", run_name=a.run, out_dir=a.out or "", buildings=a.buildings or None)


def cmd_viewer(a):
    from . import viewer
    viewer.build(_site(a), maps_dir=a.maps or "", run_name=a.run, out_path=a.out or "")


def cmd_validate(a):
    from . import validate
    validate.validate_site(_site(a), maps_dir=a.maps or "", out_dir=a.out or "")


def cmd_compare(a):
    s = _site(a)
    from .ras import compare
    import json
    hd = s.handoff_dir()
    if not hd:
        raise SystemExit("No handoff found in %s" % os.path.join(s.out_dir, "hecras"))
    cid = json.load(open(os.path.join(hd, "handoff_summary.json")))["complex"]
    runs = s.file("runs_dir") or os.path.join(s.out_dir, "hecras", "runs")
    maps = a.maps or os.path.join(runs, a.run, "maps")
    if not os.path.isdir(maps):
        raise SystemExit("No maps folder %s: give --run R04 (a finished run) or --maps <folder>" % maps)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(maps)), "compare")
    compare.compare(hd, maps, os.path.join(s.out_dir, "risk", "ahp_C%d_risk_abs.tif" % cid), out,
                    dem_file=os.path.join(s.out_dir, "risk", "ahp_C%d_dem_clip.tif" % cid), ahp_breaks=tuple(s.P["RISK_CLASS_BREAKS"]),
                    hr_breaks=tuple(s.score["hazard_breaks"]), wet_depth=float(s.score["wet_depth_m"]), debris_factor=float(s.score["debris_factor"]),
                    title="%s - %s: AHP risk map vs HEC-RAS simulation" % (s.name, a.run if not a.maps else os.path.basename(maps)))


def main(argv=None):
    p = argparse.ArgumentParser(prog="fanrisk", description="Alluvial-fan flash-flood risk mapping from a DEM (version %s)." % __version__)
    sub = p.add_subparsers(dest="cmd")
    q = sub.add_parser("init", help="make a folder with a site.toml template")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_init)
    sub.add_parser("profiles", help="list region profiles and validation cases").set_defaults(fn=cmd_profiles)
    for name, fn, hlp in (("check", cmd_check, "check the site file, inputs and installed libraries"), ("run", cmd_run, "fans -> risk maps -> HEC-RAS handoff")):
        q = sub.add_parser(name, help=hlp)
        q.add_argument("site")
        q.set_defaults(fn=fn)
    q = sub.add_parser("hecras", help="HEC-RAS: fit | check | smoke | run | maps | score | summary | extend-outflow")
    q.add_argument("action", choices=["fit", "check", "smoke", "run", "maps", "score", "summary", "extend-outflow"])
    q.add_argument("site")
    q.add_argument("--template", help="fit: an existing HEC-RAS .prj to copy and re-target (default: the bundled template)")
    q.add_argument("--out", help="fit: folder for the new model")
    q.add_argument("--spacing", type=float, help="fit: computation-point spacing in metres")
    q.add_argument("--metres", default=500, help="extend-outflow: length to add")
    q.add_argument("--runs", help="comma-separated run names, e.g. R01,R04")
    q.set_defaults(fn=cmd_hecras)
    q = sub.add_parser("score", help="score a folder of exported Max maps")
    q.add_argument("site")
    q.add_argument("maps")
    q.add_argument("--out")
    q.add_argument("--geometry", help="the model's .g01 (its boundary lines are then used for the edge check)")
    q.set_defaults(fn=cmd_score)
    q = sub.add_parser("validate", help="score the risk maps (and HEC-RAS maps) against the observed flood of the validation case")
    q.add_argument("site")
    q.add_argument("--maps", help="a folder of HEC-RAS Max maps to score as well")
    q.add_argument("--out")
    q.set_defaults(fn=cmd_validate)
    q = sub.add_parser("exposure", help="buildings, roads and critical facilities (OpenStreetMap) in the risk / HEC-RAS / observed flood maps")
    q.add_argument("site")
    q.add_argument("--osm", help="an .osm XML file (default: map.osm next to the DEM)")
    q.add_argument("--buildings", nargs="*", help="Microsoft building files (default: ms_buildings_*.csv.gz next to the DEM)")
    q.add_argument("--run", default="R04", help="HEC-RAS run whose hazard map to use (default R04)")
    q.add_argument("--maps", help="or: any folder of HEC-RAS Max maps")
    q.add_argument("--out")
    q.set_defaults(fn=cmd_exposure)
    q = sub.add_parser("viewer", help="one HTML web map of the results (risk, HEC-RAS depth / hazard, observed flood, exposure)")
    q.add_argument("site")
    q.add_argument("--run", default="R04", help="HEC-RAS run to show (default R04)")
    q.add_argument("--maps", help="or: any folder of HEC-RAS Max maps")
    q.add_argument("--out", help="the .html file to write (default runs/<site>/viewer.html)")
    q.set_defaults(fn=cmd_viewer)
    q = sub.add_parser("compare", help="where the AHP map and a HEC-RAS run agree, with numbered hotspots (PNG, GeoTIFF, KML)")
    q.add_argument("site")
    q.add_argument("--run", default="R04", help="run name in the runs folder (default R04)")
    q.add_argument("--maps", help="or: any folder of HEC-RAS Max maps")
    q.add_argument("--out")
    q.set_defaults(fn=cmd_compare)
    a = p.parse_args(argv)
    if not a.cmd:
        p.print_help()
        return 0
    try:
        a.fn(a)
    except config.ConfigError as e:
        raise SystemExit("SITE FILE PROBLEM: %s" % e)
    return 0
