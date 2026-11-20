# -*- coding: utf-8 -*-
"""Max depth / velocity maps straight from a computed plan's results HDF, without RAS Mapper.

Used when RAS Mapper's export cannot run (e.g. Windows security blocks ras-commander's RasStoreMapHelper.exe). Reads, from <project>.p<NN>.hdf:
    Geometry/2D Flow Areas/<area>/Cells Center Coordinate, Cells Minimum Elevation, Faces Cell Indexes
    Results/.../Summary Output/2D Flow Areas/<area>/Maximum Water Surface, Maximum Face Velocity
and writes, on the model terrain's grid:
    Depth (Max).tif     max water surface (inverse-distance from the 3 nearest cell centres) minus the terrain; NaN where dry
    Velocity (Max).tif  per cell: the largest of the max |velocity| of its faces. HDF face velocities are normal to the face, so a face
                        parallel to the flow shows ~0; the cell's largest face value is the closest simple match to RAS Mapper's
                        point speed (the mean of the faces is also written, as speed_cellmean_alt.tif, and checked)
Both are maxima over the whole run (the profile 'Max'), like RAS Mapper's Depth (Max) / Velocity (Max).

    python -m fanrisk.ras.hdfmaps <plan .hdf> <terrain .tif> <out folder> [<RAS Mapper maps folder to compare with>]
"""
import os
import sys

import numpy as np

from .. import gio

WET_M = 0.001


def _datasets(h5):
    out = {}

    def visit(name, obj):
        if hasattr(obj, "shape"):
            out[name] = obj
    h5.visititems(visit)
    return out


def _row0(a):
    a = np.asarray(a)
    if a.ndim == 2 and a.shape[0] == 2 and a.shape[1] != 2:
        return a[0]
    if a.ndim == 2 and a.shape[1] == 2:
        return a[:, 0]
    return a.ravel()


def read_mesh_max(hdf_path, log=print):
    """-> list of dicts per 2D area: xy (n,2), zmin (n,), wse (n,), vcell (n,) - only real (non-ghost) cells."""
    import h5py
    areas = []
    with h5py.File(hdf_path, "r") as h:
        ds = _datasets(h)
        names = sorted({k.split("/")[2] for k in ds if k.startswith("Geometry/2D Flow Areas/") and k.endswith("/Cells Center Coordinate")})
        if not names:
            raise RuntimeError("no 2D flow area geometry in %s" % hdf_path)
        for a in names:
            g = "Geometry/2D Flow Areas/%s/" % a
            xy = np.asarray(ds[g + "Cells Center Coordinate"][()], float)
            zmin = np.asarray(ds[g + "Cells Minimum Elevation"][()], float) if g + "Cells Minimum Elevation" in ds else np.full(len(xy), np.nan)
            wk = [k for k in ds if k.startswith("Results/") and "Summary Output" in k and k.endswith("/2D Flow Areas/%s/Maximum Water Surface" % a)]
            vk = [k for k in ds if k.startswith("Results/") and "Summary Output" in k and k.endswith("/2D Flow Areas/%s/Maximum Face Velocity" % a)]
            if not wk:
                raise RuntimeError("no 'Maximum Water Surface' results for area '%s' in %s - has the plan been computed?" % (a, hdf_path))
            wse = _row0(ds[wk[0]][()]).astype(float)
            n = min(len(xy), len(wse), len(zmin))
            xy, zmin, wse = xy[:n], zmin[:n], wse[:n]
            vcell = np.full(n, np.nan)
            vmean = np.full(n, np.nan)
            if vk and g + "Faces Cell Indexes" in ds:
                fv = np.abs(_row0(ds[vk[0]][()]).astype(float))
                fc = np.asarray(ds[g + "Faces Cell Indexes"][()], int)[:len(fv)]
                s, c, m = np.zeros(n), np.zeros(n), np.full(n, -np.inf)
                for j in (0, 1):
                    i = fc[:, j]
                    ok = (i >= 0) & (i < n) & np.isfinite(fv)
                    np.add.at(s, i[ok], fv[ok])
                    np.add.at(c, i[ok], 1)
                    np.maximum.at(m, i[ok], fv[ok])
                vmean = np.where(c > 0, s / np.maximum(c, 1), np.nan)
                vcell = np.where(c > 0, m, np.nan)
            real = np.isfinite(zmin) & np.isfinite(wse)
            log("   HDF: area '%s': %d cells (%d real), max WSE read%s" % (a, n, int(real.sum()), ", face velocities read" if vk else ", NO face velocities"))
            areas.append(dict(name=a, xy=xy[real], zmin=zmin[real], wse=wse[real], vcell=vcell[real], vmean=vmean[real]))
    return areas


def maps_from_hdf(hdf_path, terrain_tif, out_dir, log=print):
    from scipy.spatial import cKDTree
    z, gt, crs = gio.read_raster(terrain_tif)
    H, W = z.shape
    X = gt[0] + (np.arange(W) + 0.5) * gt[1]
    Y = gt[3] + (np.arange(H) + 0.5) * gt[5]
    XX, YY = np.meshgrid(X, Y)
    P = np.c_[XX.ravel(), YY.ravel()]
    depth = np.full(z.size, np.nan)
    vel = np.full(z.size, np.nan)
    vel2 = np.full(z.size, np.nan)
    best = np.full(z.size, np.inf)
    for A in read_mesh_max(hdf_path, log):
        if len(A["xy"]) < 3:
            continue
        tree = cKDTree(A["xy"])
        d, i = tree.query(P, k=3)
        spacing = float(np.median(tree.query(A["xy"], k=2)[0][:, 1]))
        near = d[:, 0] <= 0.75 * spacing + 0.5 * abs(gt[1])               # inside the mesh (the pixel's own cell is the nearest centre)
        wet_cell = (A["wse"] - A["zmin"]) > WET_M
        w = 1.0 / np.maximum(d, 1e-3) ** 2
        wse = (A["wse"][i] * w).sum(1) / w.sum(1)
        dep = wse - z.ravel()
        use = near & (d[:, 0] < best) & wet_cell[i[:, 0]] & np.isfinite(dep) & (dep > WET_M)
        depth[use] = dep[use]
        vel[use] = A["vcell"][i[use, 0]]
        vel2[use] = A.get("vmean", A["vcell"])[i[use, 0]]
        best[near] = np.minimum(best[near], d[near, 0])
    depth, vel = depth.reshape(z.shape), vel.reshape(z.shape)
    os.makedirs(out_dir, exist_ok=True)
    fd = os.path.join(out_dir, "Depth (Max).tif")
    fv = os.path.join(out_dir, "Velocity (Max).tif")
    gio.write_raster(fd, depth, gt, crs, "float32", gio.NODATA)
    gio.write_raster(fv, vel, gt, crs, "float32", gio.NODATA)
    gio.write_raster(os.path.join(out_dir, "speed_cellmean_alt.tif"), vel2.reshape(z.shape), gt, crs, "float32", gio.NODATA)
    with open(os.path.join(out_dir, "maps_from_hdf.txt"), "w") as fh:
        fh.write("Depth (Max) and Velocity (Max) made by fanrisk.ras.hdfmaps from %s on the grid of %s (RAS Mapper was not used).\n" % (hdf_path, terrain_tif))
    wet = np.isfinite(depth)
    log("   maps from the results HDF: wet area %.3f km2, max depth %.2f m, max speed %.2f m/s -> %s" % (
        wet.sum() * abs(gt[1] * gt[5]) / 1e6, np.nanmax(depth) if wet.any() else 0, np.nanmax(vel) if np.isfinite(vel).any() else 0, out_dir))
    return [fd, fv]


def compare(out_dir, ref_dir, log=print):
    """Agreement of these maps with RAS Mapper's Max maps of the same run (a check of this exporter)."""
    from . import readback as rb
    res = {}
    for key, fn in (("depth", "Depth (Max).tif"), ("velocity", "Velocity (Max).tif"), ("velocity (cell mean)", "speed_cellmean_alt.tif")):
        if not os.path.exists(os.path.join(out_dir, fn)):
            continue
        a, agt, _ = gio.read_raster(os.path.join(out_dir, fn))
        rp, _k = rb.find_tif(ref_dir, key.split()[0], "", True)
        r, rgt, _ = gio.read_raster(rp)
        Xg, Yg = rb.grid_xy(agt, a.shape)
        rs = rb.sample_nearest(r, rgt, Xg, Yg)
        wa, wr = np.isfinite(a) & (np.nan_to_num(a) >= (0.3 if key == "depth" else 0)), np.isfinite(rs) & (np.nan_to_num(rs) >= (0.3 if key == "depth" else 0))
        both = np.isfinite(a) & np.isfinite(rs)
        d = (a - rs)[both]
        res[key] = dict(bias=float(np.mean(d)) if d.size else np.nan, mae=float(np.mean(np.abs(d))) if d.size else np.nan,
                        wet_iou=float((wa & wr).sum() / max(1, (wa | wr).sum())))
        log("   check vs RAS Mapper %-20s: mean diff %+.2f, mean |diff| %.2f%s" % (key, res[key]["bias"], res[key]["mae"],
                                                                          (", wet (>=0.3 m) overlap IoU %.2f" % res[key]["wet_iou"]) if key == "depth" else ""))
    return res


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    maps_from_hdf(sys.argv[1], sys.argv[2], sys.argv[3])
    if len(sys.argv) > 4:
        compare(sys.argv[3], sys.argv[4])
