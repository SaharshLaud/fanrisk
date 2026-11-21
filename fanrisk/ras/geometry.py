# -*- coding: utf-8 -*-
"""HEC-RAS 2D model files as plain text: fit an existing project to a new handoff without opening RAS Mapper.

    fit_project()      copy a template project and re-target it: perimeter, computation points and boundary lines from the handoff,
                       outflow friction slope, shallow-water equation set, projection, and (with ras-commander) the terrain.
    extend_outflow()   lengthen the outflow boundary line along the perimeter (when water ponds against the model edge).
    read_bc_lines()    the boundary lines of a geometry file (used by scoring).

HEC-RAS's own geometry preprocessor builds the mesh from the computation points on the first run. Only the text files are edited; the
template project is never changed. Works on a project with ONE 2D flow area and an 'inflow' and an 'outflow' boundary line.
"""
import datetime
import glob
import json
import math
import os
import re
import shutil

import numpy as np


def _read(p):
    with open(p, encoding="latin-1", newline="") as fh:
        raw = fh.read()
    eol = "\r\n" if "\r\n" in raw else "\n"
    L = raw.split(eol)
    if L and L[-1] == "":
        L = L[:-1]
    return L, eol


def _write(p, L, eol):
    with open(p, "w", encoding="latin-1", newline="") as fh:
        fh.write(eol.join(L) + eol)


def f16(v):
    return ("%.9f" % v)[:16].rjust(16)


def geometry_file(project_dir):
    g = sorted(p for p in glob.glob(os.path.join(project_dir, "*.g[0-9][0-9]")))
    if not g:
        raise RuntimeError("no geometry file (*.g01) in %s" % project_dir)
    return g[0]


def _perimeter(L):
    i0 = next(i for i, l in enumerate(L) if l.startswith("Storage Area Surface Line="))
    n0 = int(L[i0].split("=")[1])
    return i0, n0, np.array([[float(l[0:16]), float(l[16:32])] for l in L[i0 + 1:i0 + 1 + n0]])


def _bc_blocks(L):
    """name -> (index of 'BC Line Name', index of 'BC Line Arc', number of points, points)."""
    out = {}
    for ib, l in enumerate(L):
        if not l.startswith("BC Line Name="):
            continue
        name = l.split("=", 1)[1].strip()
        ia = next(i for i in range(ib, len(L)) if L[i].startswith("BC Line Arc="))
        na = int(L[ia].split("=")[1])
        P = []
        for k in range(int(math.ceil(na / 2.0))):
            s = L[ia + 1 + k]
            P += [(float(s[c:c + 16]), float(s[c + 16:c + 32])) for c in range(0, len(s.rstrip()), 32) if s[c:c + 16].strip()]
        out[name] = (ib, ia, na, P[:na])
    return out


def read_bc_lines(g_path):
    """{name: [(x, y), ...]} of every boundary-condition line in a geometry file."""
    L, _ = _read(g_path)
    return {k: v[3] for k, v in _bc_blocks(L).items()}


def _dp(P, tol):
    keep = np.zeros(len(P), bool)
    keep[0] = keep[-1] = True
    st = [(0, len(P) - 1)]
    while st:
        i, j = st.pop()
        if j <= i + 1:
            continue
        A, B = P[i], P[j]
        v = B - A
        L2 = v @ v
        Q = P[i + 1:j] - A
        t = np.clip(Q @ v / L2, 0, 1) if L2 > 0 else np.zeros(len(Q))
        d = np.hypot(*(Q - np.outer(t, v)).T)
        k = int(d.argmax())
        if d[k] > tol:
            keep[i + 1 + k] = True
            st += [(i, i + 1 + k), (i + 1 + k, j)]
    return P[keep]


def _simple(R):
    A, B = R[:-1], R[1:]
    n = len(A)
    ccw = lambda a, b, c: (c[..., 1] - a[..., 1]) * (b[..., 0] - a[..., 0]) > (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
    for i in range(n):
        j = np.arange(i + 2, n)
        j = j[~((i == 0) & (j == n - 1))]
        if j.size and np.any((ccw(A[i][None], A[j], B[j]) != ccw(B[i][None], A[j], B[j])) & (ccw(A[i][None], B[i][None], A[j]) != ccw(A[i][None], B[i][None], B[j]))):
            return False
    return True


def _area(R):
    return 0.5 * float(np.sum(R[:-1, 0] * R[1:, 1] - R[1:, 0] * R[:-1, 1]))


def fit_geometry(g_path, handoff_dir, spacing, log=print):
    """Rewrite the geometry file in place for the handoff's flow area: perimeter (simplified), a regular grid of computation points at
    `spacing` metres, and the inflow / outflow lines laid along the perimeter. Returns facts about the result."""
    L, EOL = _read(g_path)
    ring0 = np.array(json.load(open(os.path.join(handoff_dir, "flow_area_polygon.geojson")))["features"][0]["geometry"]["coordinates"][0], float)
    open_ring = ring0[:-1]
    k = int(np.argmax(np.hypot(*(open_ring - open_ring[0]).T)))
    for tol in (25.0, 18.0, 12.0, 7.0, 3.0, 0.0):
        a = _dp(open_ring[:k + 1], tol) if tol > 0 else open_ring[:k + 1]
        b = _dp(np.vstack([open_ring[k:], open_ring[:1]]), tol) if tol > 0 else np.vstack([open_ring[k:], open_ring[:1]])
        R = np.vstack([a[:-1], b])
        if _simple(R):
            break
    else:
        raise RuntimeError("the flow-area polygon of the handoff crosses itself: fix flow_area_polygon.geojson first")
    i0, n0, old = _perimeter(L)
    if np.sign(_area(R)) != np.sign(_area(np.vstack([old, old[:1]]))):
        R = R[::-1]
    log("perimeter: %d vertices (from %d), tolerance %.0f m, area %.2f km2" % (len(R), len(ring0), tol, abs(_area(R)) / 1e6))
    x0, x1, y0, y1 = R[:, 0].min(), R[:, 0].max(), R[:, 1].min(), R[:, 1].max()
    xs = np.arange(math.floor(x0 / spacing) * spacing + spacing / 2, x1, spacing)
    ys = np.arange(math.ceil(y1 / spacing) * spacing - spacing / 2, y0, -spacing)
    X, Y = np.meshgrid(xs, ys)
    ins = np.zeros(X.shape, bool)
    for (a1, b1), (a2, b2) in zip(R[:-1], R[1:]):
        if b1 != b2:
            ins ^= ((b1 > Y) != (b2 > Y)) & (X < (a2 - a1) * (Y - b1) / (b2 - b1) + a1)
    dmin = np.full(X.shape, np.inf)
    for (a1, b1), (a2, b2) in zip(R[:-1], R[1:]):
        vx, vy = a2 - a1, b2 - b1
        t = np.clip(((X - a1) * vx + (Y - b1) * vy) / (vx * vx + vy * vy), 0, 1)
        dmin = np.minimum(dmin, np.hypot(X - (a1 + t * vx), Y - (b1 + t * vy)))
    sel = ins & (dmin >= 0.5 * spacing)
    pts = np.column_stack([X[sel], Y[sel]])
    log("computation points: %d at %.0f m spacing" % (len(pts), spacing))

    def bc_path(line):
        A, B = np.array(line[0]), np.array(line[-1])
        Rr = R[:-1]
        n = len(Rr)
        ia, ib = int(np.argmin(np.hypot(*(Rr - A).T))), int(np.argmin(np.hypot(*(Rr - B).T)))
        walk = lambda i, j: [Rr[(i + s) % n] for s in range((j - i) % n + 1)]
        p1, p2 = walk(ia, ib), walk(ib, ia)
        ln = lambda p: sum(math.hypot(p[q + 1][0] - p[q][0], p[q + 1][1] - p[q][1]) for q in range(len(p) - 1))
        return np.array(p1 if ln(p1) <= ln(p2) else p2)

    bcs, facts = {}, dict(perimeter_vertices=len(R), area_km2=abs(_area(R)) / 1e6, points=len(pts), spacing_m=spacing)
    for nm, f in (("inflow", "bc_inflow"), ("outflow", "bc_outflow")):
        p = bc_path(json.load(open(os.path.join(handoff_dir, f + ".geojson")))["features"][0]["geometry"]["coordinates"])
        bcs[nm] = p
        facts[nm + "_m"] = float(sum(np.hypot(*np.diff(p, axis=0).T)))
        log("BC line %-7s %d vertices, %.0f m long" % (nm, len(p), facts[nm + "_m"]))
    have = set(_bc_blocks(L))
    if not {"inflow", "outflow"} <= {h.lower() for h in have}:
        raise RuntimeError("the template geometry needs boundary lines named 'inflow' and 'outflow' (found: %s)" % (", ".join(sorted(have)) or "none"))

    out, i = [], 0
    now = datetime.datetime.now().strftime("%d%b%Y %H:%M:%S")
    while i < len(L):
        l = L[i]
        if l.startswith("Viewing Rectangle="):
            out.append("Viewing Rectangle= %s , %s , %s , %s " % (x0 - 100, x1 + 100, y1 + 100, y0 - 100))
            i += 1
            continue
        if l.startswith("Storage Area Surface Line="):
            out.append("Storage Area Surface Line= %d " % len(R))
            out += [f16(a) + f16(b) for a, b in R]
            i += 1 + n0
            continue
        if l.startswith("Storage Area Point Generation Data="):
            out.append("Storage Area Point Generation Data=,,%g,%g" % (spacing, spacing))
            i += 1
            continue
        if l.startswith("Storage Area 2D Points="):
            n = int(l.split("=")[1])
            out.append("Storage Area 2D Points= %d " % len(pts))
            out += ["".join(f16(v) for xy in pts[q:q + 2] for v in xy) for q in range(0, len(pts), 2)]
            i += 1 + int(math.ceil(n / 2.0))
            continue
        if l.startswith("Storage Area 2D PointsPerimeterTime="):
            out.append("Storage Area 2D PointsPerimeterTime=" + now)
            i += 1
            continue
        if l.startswith("BreakLine Name="):          # breaklines only steer RAS Mapper's point generator: dropped
            i += 1
            while i < len(L) and not L[i].startswith(("BreakLine Name=", "BC Line Name=", "LCMann")):
                i += 1
            continue
        if l.startswith("BC Line Name=") and l.split("=", 1)[1].strip().lower() in bcs:
            p = bcs[l.split("=", 1)[1].strip().lower()]
            cum = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
            mid = np.array([np.interp(cum[-1] / 2, cum, p[:, 0]), np.interp(cum[-1] / 2, cum, p[:, 1])])
            out += [l, L[i + 1], "BC Line Start Position= %s , %s " % (p[0][0], p[0][1]), "BC Line Middle Position= %s , %s " % (mid[0], mid[1]),
                    "BC Line End Position= %s , %s " % (p[-1][0], p[-1][1]), "BC Line Arc= %d " % len(p)]
            out += ["".join(f16(v) for xy in p[q:q + 2] for v in xy) for q in range(0, len(p), 2)]
            j = i + 1
            while not L[j].startswith("BC Line Text Position="):
                j += 1
            i = j
            continue
        out.append(l)
        i += 1
    _write(g_path, out, EOL)
    return facts


def extend_outflow(g_path, extra_m, name="outflow", log=print):
    """Lengthen a boundary line by `extra_m` along the perimeter, beyond its end, following the perimeter round its corners (at most one
    vertex is inserted) and stopping short of any other boundary line. The perimeter shape and the computation points do not change. Returns (old length, new length)."""
    L, EOL = _read(g_path)
    i0, n0, Rarr = _perimeter(L)
    R = [tuple(r) for r in Rarr]
    blocks = _bc_blocks(L)
    key = next((k for k in blocks if k.lower() == name.lower()), None)
    if key is None:
        raise RuntimeError("no boundary line named '%s' (found: %s)" % (name, ", ".join(blocks)))
    ib, ia, na, P = blocks[key]
    near = lambda p: min(range(n0), key=lambda j: math.hypot(R[j][0] - p[0], R[j][1] - p[1]))
    idx = [near(p) for p in P]
    step = 1 if (idx[1] - idx[0]) % n0 == 1 else -1
    taken = set()                                          # perimeter vertices used by the OTHER boundary lines: never run into them
    for k2, (_, _, _, P2) in blocks.items():
        if k2 != key:
            taken.update(near(p) for p in P2)
    e, left, extra_idx = idx[-1], float(extra_m), []
    while True:                                            # walk along the perimeter, corner by corner, until the length is used up
        j = (e + step) % n0
        (xa, ya), (xb, yb) = R[e], R[j]
        seg = math.hypot(xb - xa, yb - ya)
        if j in taken or j in idx:
            seg_cap = max(0.0, seg - 1.0)                  # stop just short of the next line
            if left > seg_cap:
                log("%s line: the next boundary line is %.0f m further on - extension cut from %.0f m to %.0f m" % (
                    key, seg, extra_m, extra_m - left + seg_cap))
                left = seg_cap
            break
        if left < seg - 1:
            break
        extra_idx.append(j)
        left -= seg
        e = j
    if left < 1.0 and not extra_idx:
        raise RuntimeError("no room to extend the %s line along the perimeter" % key)
    if left >= 1.0:
        t = left / seg
        nv = (xa + t * (xb - xa), ya + t * (yb - ya))
        ins = max(e, j) if abs(e - j) == 1 else 0             # the wrap-around edge (last vertex <-> first): append at the end
        R2 = R[:ins] + [nv] + R[ins:] if ins else R + [nv]
        newP = [R[q] for q in idx] + [R[q] for q in extra_idx] + [nv]
    else:                                                  # the length ends on a corner: no new vertex
        R2 = list(R)
        newP = [R[q] for q in idx] + [R[q] for q in extra_idx]
    cum = [0.0]
    for q in range(1, len(newP)):
        cum.append(cum[-1] + math.hypot(newP[q][0] - newP[q - 1][0], newP[q][1] - newP[q - 1][1]))
    h = cum[-1] / 2
    q = next(q for q in range(1, len(cum)) if cum[q] >= h)
    u = (h - cum[q - 1]) / (cum[q] - cum[q - 1])
    mid = (newP[q - 1][0] + u * (newP[q][0] - newP[q - 1][0]), newP[q - 1][1] + u * (newP[q][1] - newP[q - 1][1]))
    out = L[:i0] + ["Storage Area Surface Line= %d " % len(R2)] + [f16(x) + f16(y) for x, y in R2] + L[i0 + 1 + n0:]
    d = len(R2) - n0
    ib, ia = ib + d, ia + d
    blk = ["BC Line Start Position= %s , %s " % newP[0], "BC Line Middle Position= %s , %s " % mid, "BC Line End Position= %s , %s " % newP[-1],
           "BC Line Arc= %d " % len(newP)] + ["".join(f16(v) for xy in newP[m:m + 2] for v in xy) for m in range(0, len(newP), 2)]
    js = next(i for i in range(ib, len(out)) if out[i].startswith("BC Line Start Position="))
    je = ia + 1 + int(math.ceil(na / 2.0))
    out = out[:js] + blk + out[je:]
    _write(g_path, out, EOL)
    old_len = cum[len(idx) - 1]
    log("%s line: %.0f m -> %.0f m (%d -> %d vertices, %d perimeter corner(s) passed); perimeter %d -> %d vertices" % (
        key, old_len, cum[-1], na, len(newP), len(extra_idx), n0, len(R2)))
    return old_len, cum[-1]


def _set_line(path, prefix, new_line):
    L, eol = _read(path)
    hit = [i for i, l in enumerate(L) if l.startswith(prefix)]
    for i in hit:
        L[i] = new_line
    if hit:
        _write(path, L, eol)
    return len(hit)


def fit_project(template_prj, handoff_dir, dst_dir, base_plan="01", spacing=None, swe=True, version="7.0.2", terrain=True, mannings_n=None, log=print):
    """Copy the template HEC-RAS project (text files and land-cover folder, no results) to dst_dir and fit it to the handoff.
    Returns the new .prj path. The mesh itself is built by HEC-RAS on the first run ('fanrisk hecras smoke')."""
    src = os.path.dirname(os.path.abspath(template_prj))
    stem = os.path.splitext(os.path.basename(template_prj))[0]
    if os.path.exists(dst_dir) and os.listdir(dst_dir):
        raise RuntimeError("%s already exists and is not empty: choose another folder (nothing is overwritten)" % dst_dir)
    os.makedirs(dst_dir, exist_ok=True)
    with open(os.path.join(handoff_dir, "handoff_summary.json")) as fh:
        summ = json.load(fh)
    plan_file = "%s.p%s" % (stem, base_plan)
    if not os.path.exists(os.path.join(src, plan_file)):
        raise RuntimeError("plan %s not found in the template project" % plan_file)
    pl, _ = _read(os.path.join(src, plan_file))
    geom = next(l.split("=")[1].strip() for l in pl if l.startswith("Geom File="))
    flow = next(l.split("=")[1].strip() for l in pl if l.startswith("Flow File="))
    for ext in ("prj", "p" + base_plan, geom, flow, "rasmap", "dss", "x" + geom[1:]):
        p = os.path.join(src, "%s.%s" % (stem, ext))
        if os.path.exists(p):
            shutil.copy2(p, os.path.join(dst_dir, os.path.basename(p)))
    for folder in os.listdir(src):
        if os.path.isdir(os.path.join(src, folder)) and folder.lower().startswith("land"):
            shutil.copytree(os.path.join(src, folder), os.path.join(dst_dir, folder))
    prj = os.path.join(dst_dir, stem + ".prj")
    L, eol = _read(prj)                                   # keep only the base plan, its geometry and its flow file
    L = [l for l in L if not ((l.startswith("Plan File=") and l.split("=")[1].strip() != "p" + base_plan) or
                              (l.startswith("Geom File=") and l.split("=")[1].strip() != geom) or
                              (l.startswith("Unsteady File=") and l.split("=")[1].strip() != flow))]
    _write(prj, L, eol)
    g_path = os.path.join(dst_dir, "%s.%s" % (stem, geom))
    if spacing is None:
        spacing = max(float(summ.get("settings", {}).get("RAS_CELL_SIZE_M", 10.0)), 0.5 * float(summ.get("terrain_cell_m", 20.0)))
    facts = fit_geometry(g_path, handoff_dir, float(spacing), log)
    slope = float(summ["outflow_friction_slope"])
    n = _set_line(os.path.join(dst_dir, "%s.%s" % (stem, flow)), "Friction Slope=", "Friction Slope=%.4f" % slope)
    log("outflow friction slope %.4f written to %s.%s (%d line)" % (slope, stem, flow, n))
    if swe:
        _set_line(os.path.join(dst_dir, plan_file), "UNET D2 Equation=", "UNET D2 Equation= 1 ")
        log("equation set: full shallow-water equations (SWE-ELM)")

    # roughness: the template's land-cover table carries no Manning's n values, so HEC-RAS uses the 2D area's base n everywhere.
    Lg, eg = _read(g_path)
    old_n = next((l.split("=", 1)[1].strip() for l in Lg if l.startswith("Storage Area Mannings=")), "")
    if mannings_n is not None:
        _set_line(g_path, "Storage Area Mannings=", "Storage Area Mannings=%g" % float(mannings_n))
    n_used = float(mannings_n) if mannings_n is not None else (float(old_n) if old_n else float("nan"))
    log("roughness: one Manning's n = %g for the whole 2D area%s (per-class land-cover roughness is not applied)" % (
        n_used, "" if mannings_n is not None else " (from the template)"))
    epsg = summ.get("epsg")
    proj = os.path.join(dst_dir, "projection.prj")
    src_prj = os.path.join(handoff_dir, "flow_area_polygon.prj")
    if os.path.exists(src_prj):
        shutil.copy2(src_prj, proj)
    elif epsg:
        from .. import gio
        with open(proj, "w") as fh:
            fh.write(gio.utm_wkt(int(epsg)))
    rm = os.path.join(dst_dir, stem + ".rasmap")
    if os.path.exists(rm) and os.path.exists(proj):
        with open(rm, encoding="utf-8", errors="ignore", newline="") as fh:
            txt = fh.read()
        txt = re.sub(r'<RASProjectionFilename Filename="[^"]*"', '<RASProjectionFilename Filename=".\\\\projection.prj"', txt)
        with open(rm, "w", encoding="utf-8", newline="") as fh:
            fh.write(txt)
    facts.update(prj=prj, terrain="not built", mannings_n=n_used)
    if terrain:
        tdir = os.path.join(dst_dir, "Terrain")
        os.makedirs(tdir, exist_ok=True)
        tif = os.path.join(tdir, "terrain_handoff.tif")
        shutil.copy2(os.path.join(handoff_dir, "terrain_burned.tif" if os.path.exists(os.path.join(handoff_dir, "terrain_burned.tif")) else "terrain.tif"), tif)
        try:
            from pathlib import Path
            from ras_commander.terrain import RasTerrain
            RasTerrain.create_terrain_hdf(input_rasters=[Path(tif)], output_hdf=Path(tdir) / "Terrain.hdf", projection_prj=Path(proj), units="Meters",
                                          hecras_version=version)
            facts["terrain"] = "built"
            log("terrain: Terrain\\Terrain.hdf built from the handoff terrain")
        except Exception as e:
            log("WARNING: the terrain could not be built here (%s: %s). Build it in RAS Mapper from %s (name it 'Terrain')." % (type(e).__name__, e, tif))
    log("model written to %s - nothing in the template was changed" % dst_dir)
    return facts
