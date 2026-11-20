# -*- coding: utf-8 -*-
"""fanrisk numerical engine (numpy only).

Alluvial-fan delineation after Norini et al. (2016, Geomorphology 273, 134-149), AHP flash-flood hazard / risk on the fans and their
run-out corridors, feeder hydrology, and the HEC-RAS 2D handoff geometry. This is the tested core of norini_alluvial_fans_v5.py with the
QGIS glue removed; the only change in behaviour is that runoff can come from an SCS curve number (RUNOFF_METHOD = "scs_cn") as well as
from a runoff coefficient. Functions take the settings dict `P` (see defaults.py): nothing here reads files or knows about a place.
"""
import re
import math
import time
import heapq
from types import SimpleNamespace

import numpy as np

_T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - _T0, msg))
    try:                                  # keep the QGIS UI alive during long loops
        QCoreApplication.processEvents()  # noqa: F821  (imported in the QGIS section)
    except Exception:
        pass

# =============================================================================
# 2. NUMERICAL CORE  (pure numpy - no QGIS needed, so it can be unit-tested)
# =============================================================================


class Grid(object):
    """North-up raster geometry. Pixel (r, c) has its centre at x0+(c+.5)dx, y0-(r+.5)dy."""

    def __init__(self, gt, shape):
        self.x0, self.dx, _, self.y0, _, ndy = gt
        self.dy = -ndy
        self.H, self.W = shape

    def xy_to_rc(self, x, y):             # fractional row/col where integer = pixel centre
        return (self.y0 - y) / self.dy - 0.5, (x - self.x0) / self.dx - 0.5

    def rc_to_xy(self, r, c):
        return self.x0 + (c + 0.5) * self.dx, self.y0 - (r + 0.5) * self.dy

    def cell_of(self, x, y):
        return (int(math.floor((self.y0 - y) / self.dy)),
                int(math.floor((x - self.x0) / self.dx)))


def odd(n):
    n = max(1, int(n))
    return n if n % 2 == 1 else n + 1


def border_cells(valid):
    """Valid cells that touch an invalid cell or the raster edge (8-neighbourhood)."""
    H, W = valid.shape
    vp = np.pad(valid, 1, constant_values=False)
    inner = np.ones_like(valid)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr or dc:
                inner &= vp[1 + dr:1 + dr + H, 1 + dc:1 + dc + W]
    return valid & ~inner


def erode(valid, n):
    v = valid.copy()
    for _ in range(int(n)):
        v &= ~border_cells(v)
    return v


def fill_nan(z):
    f = np.isfinite(z)
    return np.where(f, z, np.nanmean(z[f]) if f.any() else 0.0)


def disk_mean(a, R):
    """NaN-aware mean over a circular window of radius R cells (the paper's 'convolution filter')."""
    a = np.asarray(a, dtype=np.float64)
    R = int(R)
    if R < 1:
        return a.copy()
    fin = np.isfinite(a)
    H, W = a.shape
    s = np.pad(np.where(fin, a, 0.0), R)
    n = np.pad(fin.astype(np.float64), R)
    cs = np.pad(np.cumsum(s, axis=1), ((0, 0), (1, 0)))
    cn = np.pad(np.cumsum(n, axis=1), ((0, 0), (1, 0)))
    tot = np.zeros((H, W))
    cnt = np.zeros((H, W))
    for dr in range(-R, R + 1):
        w = int(math.floor(math.sqrt(R * R - dr * dr) + 1e-9))
        rows = slice(R + dr, R + dr + H)
        hi = slice(R + w + 1, R + w + 1 + W)
        lo = slice(R - w, R - w + W)
        tot += cs[rows, hi] - cs[rows, lo]
        cnt += cn[rows, hi] - cn[rows, lo]
    with np.errstate(invalid="ignore", divide="ignore"):
        out = tot / cnt
    out[(cnt == 0) | ~fin] = np.nan
    return out


def horn_slope(z, dx, dy):
    """Slope in degrees (Horn 3x3). z must not contain NaN."""
    zp = np.pad(z, 1, mode="edge")
    a, b, c = zp[:-2, :-2], zp[:-2, 1:-1], zp[:-2, 2:]
    d, f = zp[1:-1, :-2], zp[1:-1, 2:]
    g, h, i = zp[2:, :-2], zp[2:, 1:-1], zp[2:, 2:]
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8.0 * dx)
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8.0 * dy)
    return np.degrees(np.arctan(np.hypot(dzdx, dzdy)))


def surface_curvature(z, dx, dy, k=1):
    """Surface (total) curvature, ArcGIS sign convention: + convex, - concave, x100.
    k = baseline in cells; larger k measures valley-scale rather than pixel-scale curvature."""
    k = max(1, int(k))
    zp = np.pad(z, k, mode="edge")
    z5 = zp[k:-k, k:-k]
    z4, z6 = zp[k:-k, :-2 * k], zp[k:-k, 2 * k:]
    z2, z8 = zp[:-2 * k, k:-k], zp[2 * k:, k:-k]
    D = ((z4 + z6) / 2.0 - z5) / (k * dx) ** 2
    E = ((z2 + z8) / 2.0 - z5) / (k * dy) ** 2
    return -2.0 * (D + E) * 100.0


def bilinear(z, rf, cf):
    """Bilinear sample at fractional (row, col) (integer = pixel centre). NaN outside / on NaN."""
    H, W = z.shape
    r0 = np.floor(rf).astype(np.int64)
    c0 = np.floor(cf).astype(np.int64)
    fr, fc = rf - r0, cf - c0
    inside = (r0 >= 0) & (c0 >= 0) & (r0 < H - 1) & (c0 < W - 1)
    r0c = np.clip(r0, 0, H - 2)
    c0c = np.clip(c0, 0, W - 2)
    out = (z[r0c, c0c] * (1 - fr) * (1 - fc) + z[r0c, c0c + 1] * (1 - fr) * fc +
           z[r0c + 1, c0c] * fr * (1 - fc) + z[r0c + 1, c0c + 1] * fr * fc)
    return np.where(inside, out, np.nan)


def smooth1d(x, w):
    w = odd(w)
    if w <= 1:
        return x.copy()
    xp = np.pad(x, w // 2, mode="edge")
    return np.convolve(xp, np.ones(w) / w, mode="valid")


def flood(allowed, seeds, conn8=True):
    """Cells of `allowed` connected to any of `seeds` [(r, c), ...] (iterative DFS)."""
    H, W = allowed.shape
    Wp = W + 2
    a = np.zeros((H + 2, Wp), dtype=bool)
    a[1:-1, 1:-1] = allowed
    al = bytearray(a.ravel().tobytes())
    out = bytearray(len(al))
    offs = (-Wp - 1, -Wp, -Wp + 1, -1, 1, Wp - 1, Wp, Wp + 1) if conn8 else (-Wp, -1, 1, Wp)
    stack = []
    for r, c in seeds:
        i = (r + 1) * Wp + c + 1
        if al[i]:
            al[i] = 0
            out[i] = 1
            stack.append(i)
    while stack:
        i = stack.pop()
        for o in offs:
            n = i + o
            if al[n]:
                al[n] = 0
                out[n] = 1
                stack.append(n)
    return np.frombuffer(bytes(out), dtype=np.uint8).reshape(H + 2, Wp)[1:-1, 1:-1].astype(bool)


def fill_holes(mask):
    border = np.zeros_like(mask)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    bg = ~mask
    seeds = [tuple(s) for s in np.argwhere(bg & border)]
    return ~flood(bg, seeds, conn8=False)


def priority_flood(dem, valid):
    """Depression filling + flow routing in one pass (Priority-Flood, Barnes et al. 2014).
    Returns filled DEM (2-D), receiver (flat index of the downstream cell, -1 = outlet) and `order`
    (flat indices in pop order: every receiver appears before its donors)."""
    H, W = dem.shape
    Hp, Wp = H + 2, W + 2
    vp = np.zeros((Hp, Wp), dtype=bool)
    vp[1:-1, 1:-1] = valid
    zp = np.full((Hp, Wp), np.inf)
    zp[1:-1, 1:-1] = np.where(valid, dem, np.inf)
    seed_rc = np.argwhere(border_cells(valid))
    z = zp.ravel().tolist()
    closed = bytearray((~vp).ravel().tobytes())
    recv = [-1] * (Hp * Wp)
    heap = []
    cnt = 0
    for r, c in seed_rc.tolist():
        i = (r + 1) * Wp + c + 1
        closed[i] = 1
        heap.append((z[i], cnt, i))
        cnt += 1
    heapq.heapify(heap)
    order = []
    offs = (-Wp - 1, -Wp, -Wp + 1, -1, 1, Wp - 1, Wp, Wp + 1)
    push, pop = heapq.heappush, heapq.heappop
    while heap:
        zc, _, c = pop(heap)
        order.append(c)
        for o in offs:
            n = c + o
            if closed[n]:
                continue
            closed[n] = 1
            zn = z[n]
            if zn < zc:                   # depression / flat: raise to the spill level
                zn = zc
                z[n] = zc
            recv[n] = c
            cnt += 1
            push(heap, (zn, cnt, n))
    filled = np.array(z, dtype=np.float64).reshape(Hp, Wp)[1:-1, 1:-1].copy()
    filled[~valid] = np.nan
    rp = np.array(recv, dtype=np.int64).reshape(Hp, Wp)[1:-1, 1:-1]
    ok = rp >= 0
    recv_u = np.where(ok, (rp // Wp - 1) * W + (rp % Wp - 1), -1).ravel()
    op = np.array(order, dtype=np.int64)
    order_u = (op // Wp - 1) * W + (op % Wp - 1)
    return filled, recv_u, order_u


def flow_accumulation(recv, order):
    """Number of cells draining through each cell (itself included)."""
    n = recv.size
    acc = [0] * n
    ol = order.tolist()
    for c in ol:
        acc[c] = 1
    rl = recv.tolist()
    for c in reversed(ol):                # donors are always processed before their receiver
        r = rl[c]
        if r >= 0:
            acc[r] += acc[c]
    return np.array(acc, dtype=np.float64)


def build_donors(recv):
    """CSR lists of upstream neighbours: donors[starts[i]:starts[i+1]] drain into cell i."""
    n = recv.size
    src = np.flatnonzero(recv >= 0)
    dst = recv[src]
    o = np.argsort(dst, kind="stable")
    starts = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(np.bincount(dst, minlength=n), out=starts[1:])
    return starts, src[o]


def upstream_cells(cell, starts, donors):
    """Watershed analysis: all cells that drain through `cell` (flat index)."""
    parts = [np.array([cell], dtype=np.int64)]
    front = parts[0]
    while front.size:
        s, e = starts[front], starts[front + 1]
        cnt = e - s
        tot = int(cnt.sum())
        if tot == 0:
            break
        rep = np.repeat(np.arange(front.size), cnt)
        off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        front = donors[s[rep] + off]
        parts.append(front)
    return np.concatenate(parts)


def build_main_donor(recv, acc, mask):
    """For every cell: its upstream neighbour inside `mask` with the largest accumulation (-1 = none)."""
    out = np.full(recv.size, -1, dtype=np.int64)
    idx = np.flatnonzero(mask)
    r = recv[idx]
    k = r >= 0
    idx, r = idx[k], r[k]
    k = mask[r]
    idx, r = idx[k], r[k]
    if idx.size == 0:
        return out
    o = np.lexsort((acc[idx], r))
    idx, r = idx[o], r[o]
    last = np.ones(idx.size, dtype=bool)
    last[:-1] = r[1:] != r[:-1]
    out[r[last]] = idx[last]
    return out


def walk(start, nxt, values, nsteps):
    """Follow `nxt` (main_donor = upstream, recv = downstream) nsteps cells from every start cell.
    Returns (sum of finite values, number of finite values, number of cells visited)."""
    cur = start.copy()
    alive = np.ones(start.size, dtype=bool)
    tot = np.zeros(start.size)
    nval = np.zeros(start.size)
    ncell = np.zeros(start.size)
    for _ in range(int(nsteps)):
        v = values[cur]
        g = alive & np.isfinite(v)
        tot[g] += v[g]
        nval[g] += 1
        ncell[alive] += 1
        nx = nxt[cur]
        alive &= nx >= 0
        cur = np.where(alive, nx, cur)
    return tot, nval, ncell



# ---- DEM preparation and fan surface ----------------------------------------
def clean_dem(arr, declared=None, extra=(), zero_mode="auto"):
    """-> (dem with NaN for every nodata cell, note). Handles NaN, the declared nodata, common fill values and
    zero-filled margins (typical after reprojection) that border much higher terrain."""
    a = np.array(arr, dtype=np.float64)
    bad = ~np.isfinite(a) | (np.abs(a) > 1e30)
    for v in ([] if declared is None else [float(declared)]) + [float(v) for v in extra] + [-32768.0, -9999.0]:
        if abs(v) < 1e30:
            bad |= np.isclose(a, v, rtol=0.0, atol=1e-6)
    notes = []
    zero = (a == 0.0) & ~bad
    if zero_mode is not False and zero.any():
        if zero_mode is True:
            bad |= zero
            notes.append("%d zero cells treated as nodata" % int(zero.sum()))
        else:
            seeds = np.argwhere(zero & border_cells(~bad))
            if seeds.size:
                region = flood(zero, [tuple(s) for s in seeds], conn8=True)
                ring = (disk_mean(region.astype(float), 1) > 0) & ~region & ~bad
                if region.sum() >= 20 and ring.any() and np.median(a[ring]) > 50.0:
                    bad |= region
                    notes.append("%d zero-filled margin cells treated as nodata" % int(region.sum()))
    a[bad] = np.nan
    return a, "; ".join(notes)


def closing(z, r):
    """Grayscale closing (dilate then erode) with a disk of radius r cells: fills channels and gullies narrower
    than about 2r cells while leaving broad valleys alone. NaN stays NaN."""
    r = int(r)
    if r < 1:
        return z.copy()
    H, W = z.shape
    fin = np.isfinite(z)
    offs = [(dr, dc) for dr in range(-r, r + 1) for dc in range(-r, r + 1) if dr * dr + dc * dc <= r * r]

    def reduce(a, fn, pad_val):
        p = np.pad(a, r, constant_values=pad_val)
        out = a.copy()
        for dr, dc in offs:
            out = fn(out, p[r + dr:r + dr + H, r + dc:r + dc + W])
        return out

    dil = reduce(np.where(fin, z, -np.inf), np.maximum, -np.inf)
    dil = np.where(np.isfinite(dil), dil, np.inf)
    ero = reduce(dil, np.minimum, np.inf)
    return np.where(fin, ero, np.nan)


def otsu_threshold(v, lo=0.0, hi=60.0, nbins=120):
    """Threshold separating a bimodal distribution (here: smoothed slope of valley floors vs mountains)."""
    v = v[np.isfinite(v)]
    h, e = np.histogram(np.clip(v, lo, hi), bins=nbins, range=(lo, hi))
    p = h / max(h.sum(), 1)
    w0 = np.cumsum(p)
    mu = np.cumsum(p * (e[:-1] + e[1:]) / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        sb = (mu[-1] * w0 - mu) ** 2 / (w0 * (1.0 - w0))
    return float(e[1:][np.argmax(np.nan_to_num(sb))])


def mscale(m, cells, step):
    """metres, but never less than `cells` pixels (keeps the rules meaningful on coarse DEMs)."""
    return max(float(m), cells * step)


def walk_stat(start, nxt, values, nsteps, fn):
    """Follow `nxt` nsteps cells from every start cell and reduce the values with fn (np.maximum / np.minimum)."""
    cur = start.copy()
    alive = np.ones(start.size, dtype=bool)
    out = values[cur].astype(np.float64)
    for _ in range(int(nsteps)):
        nx = nxt[cur]
        alive &= nx >= 0
        cur = np.where(alive, nx, cur)
        out = np.where(alive, fn(out, values[cur]), out)
    return out


# ---- apex detection: valley mouth ("openness") --------------------------------
def openness_cells(zs, grid, cells, L, ndir=16, blocked_deg=8.0, chunk=3000):
    """Number of compass directions (of ndir) in which the terrain does NOT rise above blocked_deg within L.
    Confined valley: few open directions; fan / plain: most. Terrain outside the DEM counts as blocked."""
    W = grid.W
    step = min(grid.dx, grid.dy)
    dist = np.arange(1, max(2, int(L / step)) + 1) * step
    th = np.arange(ndir) * 2 * np.pi / ndir
    out = np.zeros(cells.size)
    zflat = zs.ravel()
    for s in range(0, cells.size, chunk):
        cs = cells[s:s + chunk]
        x, y = grid.rc_to_xy(cs // W, cs % W)
        X = x[:, None, None] + np.cos(th)[None, :, None] * dist[None, None, :]
        Y = y[:, None, None] + np.sin(th)[None, :, None] * dist[None, None, :]
        rf, cf = grid.xy_to_rc(X, Y)
        Z = bilinear(zs, rf.reshape(-1), cf.reshape(-1)).reshape(X.shape)
        ang = np.degrees(np.arctan((Z - zflat[cs][:, None, None]) / dist[None, None, :]))
        ang = np.where(np.isfinite(ang), ang, 90.0)
        out[s:s + chunk] = (ang.max(axis=2) < blocked_deg).sum(axis=1)
    return out


def detect_apexes_openness(ctx, P):
    """Fan apex = first point on a channel where the surrounding terrain has opened up to OPEN_FRACTION of what it
    will be 2 km downstream, after being walled in upstream (the channel leaves its confined valley)."""
    step = ctx.step
    sc = np.flatnonzero(ctx.streams_flat)
    if sc.size == 0:
        return []
    md = ctx.main_donor_all
    L = mscale(P["OPEN_RADIUS_M"], 20, step)
    O = np.zeros(ctx.recv.size)
    O[sc] = openness_cells(ctx.zsurf, ctx.grid, sc, L, 16, P["OPEN_BLOCK_DEG"])
    k = max(1, int(round(mscale(P["OPEN_SMOOTH_M"], 3, step) / step)))
    up, dn = walk(sc, md, O, k), walk(sc, ctx.recv, O, k)
    Os = np.zeros_like(O)
    Os[sc] = (up[0] + dn[0] - O[sc]) / np.maximum(up[1] + dn[1] - 1, 1)
    npl = max(2, int(round(mscale(P["OPEN_PLATEAU_M"], 40, step) / step)))
    plat = np.zeros_like(O)
    plat[sc] = walk_stat(sc, ctx.recv, Os, npl, np.maximum)
    omin = np.zeros_like(O)
    omin[sc] = walk_stat(sc, md, Os, npl, np.minimum)
    u = md[sc]
    thr = P["OPEN_FRACTION"] * plat[sc]
    u_os = np.where(u >= 0, Os[np.maximum(u, 0)], np.inf)         # headwater cells cannot be a crossing
    hit = (Os[sc] >= thr) & (u_os < thr) & ((plat[sc] - omin[sc]) >= P["OPEN_MIN_RISE"]) & (plat[sc] >= 8.0)
    W = ctx.grid.W
    out = []
    for c in sc[hit].tolist():
        why = "" if ctx.core_flat[c] else "near DEM edge"
        out.append(dict(row=c // W, col=c % W, cell=c, source="auto", method="openness", name="", system="",
                        mean_curv=float("nan"), feeder_len_m=float("nan"),
                        status=("candidate" if not why else "rejected: " + why)))
    return out



# ---- apex detection: the paper's steep-zone exits (section 3.1) -------------------
def detect_apexes_paper(ctx, P):
    """Fan apexes = points where a selected feeder channel leaves the steep sub-area."""
    step = ctx.step
    net = ctx.network_flat
    idx = np.flatnonzero(net)
    r = ctx.recv[idx]
    has_r = r >= 0
    nxt_in = np.zeros(idx.size, dtype=bool)
    nxt_in[has_r] = net[r[has_r]]
    exits = idx[has_r & ~nxt_in]                       # last channel cell inside the steep zone
    if exits.size == 0:
        return []
    md = build_main_donor(ctx.recv, ctx.acc, net)
    ns = max(1, int(round(P["TRUNK_LENGTH_M"] / step)))
    tot, nval, ncell = walk(exits, md, ctx.curv_flat, ns)
    mean_curv = tot / np.maximum(nval, 1)
    feeder_len = ncell * step
    out = []
    W = ctx.grid.W
    for j, c in enumerate(exits.tolist()):
        why = ""
        if not ctx.core_flat[c]:
            why = "near DEM edge"
        elif feeder_len[j] < P["MIN_FEEDER_LENGTH_M"]:
            why = "feeder channel too short"
        elif mean_curv[j] > P["CURV_THRESHOLD"]:
            why = "feeder channel not concave (curvature)"
        out.append(dict(row=c // W, col=c % W, cell=c, source="auto", name="", system="",
                        mean_curv=float(mean_curv[j]), feeder_len_m=float(feeder_len[j]),
                        status=("candidate" if not why else "rejected: " + why)))
    return out



# ---- known apexes ----------------------------------------------------------
def snap_manual(manual, grid, streams, acc, P, log):
    """Snap each manual apex to the NEAREST SIZEABLE stream cell within SNAP_RADIUS_M (>= 30 % of the largest
    accumulation in that radius, so a point beside a big river does not land on a tiny creek)."""
    H, W = grid.H, grid.W
    rad = int(math.ceil(P["SNAP_RADIUS_M"] / min(grid.dx, grid.dy)))
    st2, ac2 = streams.reshape(H, W), acc.reshape(H, W)
    out = []
    for m in manual:
        r, c = grid.cell_of(m["x"], m["y"])
        if not (0 <= r < H and 0 <= c < W):
            log("  manual apex %s lies outside the DEM - skipped" % m["name"])
            continue
        cand = dict(row=r, col=c, cell=r * W + c, source="manual", method="manual", name=m["name"],
                    system=m.get("system", ""), mean_curv=float("nan"), feeder_len_m=float("nan"), status="candidate")
        r0, r1, c0, c1 = max(0, r - rad), min(H, r + rad + 1), max(0, c - rad), min(W, c + rad + 1)
        rr, cc = np.mgrid[r0:r1, c0:c1]
        d2 = (rr - r) ** 2 + (cc - c) ** 2
        ok = st2[r0:r1, c0:c1] & (d2 <= rad * rad)
        if ok.any():
            aw = ac2[r0:r1, c0:c1]
            ok &= aw >= 0.3 * aw[ok].max()
        if not ok.any():
            cand["status"] = "rejected: no stream cell within %.0f m" % P["SNAP_RADIUS_M"]
            log("  manual apex %s: no stream cell within %.0f m (lower FLOW_ACC_THRESHOLD_KM2 or raise SNAP_RADIUS_M)"
                % (m["name"], P["SNAP_RADIUS_M"]))
        else:
            j = np.unravel_index(np.argmin(np.where(ok, d2, 1 << 30)), d2.shape)
            cand["row"], cand["col"] = int(r0 + j[0]), int(c0 + j[1])
            cand["cell"] = cand["row"] * W + cand["col"]
        out.append(cand)
    return out


# ---- radial profiles and fan footprint (paper section 3.2) --------------------
def analyse_profile(z, step, P, thr_toe):
    """One radial profile (index 0 = apex) on the fan surface. Returns (end index, reason, smoothed profile).
    The toe is a SUSTAINED drop of the along-profile slope (long window) below thr_toe; a profile also ends when it
    climbs again (opposite wall, neighbouring ridge)."""
    n = z.size
    bad = ~np.isfinite(z)
    nv = n if not bad.any() else int(np.argmax(bad))
    zs_full = np.full(n, np.nan)
    if nv < 6:
        return 0, "nodata", zs_full
    zs = smooth1d(z[:nv], odd(round(mscale(P["PROFILE_SMOOTH_M"], 3, step) / step)))
    zs_full[:nv] = zs
    g = -np.gradient(zs, step)                                 # + where the profile descends
    sl = np.degrees(np.arctan(smooth1d(g, odd(round(mscale(P["SLOPE_WINDOW_M"], 8, step) / step)))))
    cands = []
    hit = np.flatnonzero(zs - np.minimum.accumulate(zs) > P["RISE_TOLERANCE_M"])
    if hit.size:
        cands.append((int(np.argmin(zs[:hit[0] + 1])), "rise"))
    ig = int(math.ceil(mscale(P["GRACE_M"], 6, step) / step))
    run = max(1, int(round(mscale(P["MIN_RUN_M"], 6, step) / step)))
    if nv - run + 1 > ig:
        cs = np.convolve((sl < thr_toe).astype(int), np.ones(run, dtype=int), mode="valid")
        ii = np.arange(cs.size)
        ok = np.flatnonzero((cs == run) & (ii >= ig) & (zs[0] - zs[:cs.size] >= P["MIN_RADIAL_RELIEF_M"]))
        if ok.size:
            cands.append((int(ok[0]), "toe"))
    if cands:
        i_end, why = min(cands, key=lambda t: t[0])
    else:
        i_end, why = nv - 1, ("edge" if nv < n else "maxlen")
    relief = zs[0] - zs[i_end]
    if relief < P["MIN_RADIAL_RELIEF_M"] or i_end * step < mscale(P["MIN_RADIAL_LENGTH_M"], 3, step):
        return 0, why + "-short", zs_full
    if math.degrees(math.atan(relief / (i_end * step))) < 0.75 * thr_toe:
        return 0, why + "-flat", zs_full
    return i_end, why, zs_full


def radial_analysis(zsurf, grid, ax, ay, P):
    step = min(grid.dx, grid.dy)
    K = max(8, int(round(360.0 / P["RADIAL_STEP_DEG"])))
    nmax = max(6, int(P["MAX_FAN_LENGTH_M"] // step))
    ang = np.arange(K) * (2 * np.pi / K)
    dist = np.arange(nmax + 1) * step
    rf, cf = grid.xy_to_rc(ax + np.cos(ang)[:, None] * dist[None, :], ay + np.sin(ang)[:, None] * dist[None, :])
    Z = bilinear(zsurf, rf, cf)
    if not np.isfinite(Z[0, 0]):
        return None
    z_apex = float(Z[0, 0])
    # the fan's own reference slope: steep end of the mean slope over the first REF_LEN_M, over all directions
    iref = max(2, min(nmax, int(round(mscale(P["REF_LEN_M"], 20, step) / step))))
    drop = (Z[:, 0] - Z[:, iref]) / (iref * step)
    drop = np.sort(drop[np.isfinite(drop)])
    # steepest-descent sector (top eighth of the directions): at a confined valley mouth only a narrow sector descends
    s_ref = math.degrees(math.atan(max(float(drop[-max(3, K // 8):].mean()), 0.0))) if drop.size else 0.0
    thr_toe = float(np.clip(P["TOE_FRACTION"] * s_ref, P["TOE_MIN_DEG"], P["TOE_MAX_DEG"]))
    iend = np.zeros(K, dtype=int)
    why = [""] * K
    ZS = np.full((K, nmax + 1), np.nan)
    for k in range(K):
        iend[k], why[k], ZS[k] = analyse_profile(Z[k], step, P, thr_toe)
    mw = int(P["RADIAL_MEDIAN_WINDOW"])
    if mw >= 3:                                                # kill single stray radials
        h = mw // 2
        iend_s = np.floor(np.median(np.stack([np.roll(iend, s) for s in range(-h, h + 1)]), axis=0)).astype(int)
    else:
        iend_s = iend.copy()
    iend_s = np.minimum(iend_s, np.maximum(0, np.isfinite(ZS).sum(axis=1) - 1))
    reason = [why[k] if iend_s[k] == iend[k] else "median" for k in range(K)]
    live = [k for k in range(K) if iend_s[k] > 0]
    return dict(ang=ang, D=iend_s * step, Zend=ZS[np.arange(K), iend_s], reason=reason, z_apex=z_apex, step=step,
                n_fan=len(live), s_ref=s_ref, thr_toe=thr_toe,
                n_trunc=sum(1 for k in live if reason[k] == "maxlen"), n_edge=sum(1 for k in live if reason[k] == "edge"))


def build_fan(zsurf, slope, grid, ax, ay, ra, P):
    """Semi-conical surface from apex + end points, intersected with the surface (rules 2 and 3 of the paper)."""
    D, Zend, z_apex = ra["D"], ra["Zend"], ra["z_apex"]
    K = D.size
    dth = 2 * np.pi / K
    rmax = float(D.max())
    if rmax <= 0:
        return None
    ra_f, ca_f = grid.xy_to_rc(ax, ay)
    r_a, c_a = int(round(ra_f)), int(round(ca_f))
    hd = float(P["HEAD_RECOVERY_M"])
    rad = int(math.ceil(max(rmax, hd) / min(grid.dx, grid.dy))) + 3
    r0, r1 = max(0, r_a - rad), min(grid.H, r_a + rad + 1)
    c0, c1 = max(0, c_a - rad), min(grid.W, c_a + rad + 1)
    rr, cc = np.mgrid[r0:r1, c0:c1]
    ddx = grid.x0 + (cc + 0.5) * grid.dx - ax
    ddy = grid.y0 - (rr + 0.5) * grid.dy - ay
    dist = np.hypot(ddx, ddy)
    t = np.mod(np.arctan2(ddy, ddx), 2 * np.pi) / dth
    k0 = np.floor(t).astype(int) % K
    k1 = (k0 + 1) % K
    f = t - np.floor(t)
    Dth = D[k0] * (1 - f) + D[k1] * f                          # fan radius in this direction
    Zth = Zend[k0] * (1 - f) + Zend[k1] * f                    # elevation at the fan toe in this direction
    foot = (Dth > 0) & (dist <= Dth)
    relief = z_apex - Zth
    cone = z_apex - relief * np.clip(dist / np.where(Dth > 0, Dth, 1.0), 0.0, 1.0)
    tol = P["CONE_TOLERANCE_M"] + P["CONE_TOLERANCE_FRAC"] * relief
    dw = zsurf[r0:r1, c0:c1]
    m = foot & np.isfinite(dw) & (dw <= cone + tol)            # fan = where the surface lies on/below the cone
    if P["MASK_SMOOTH_CELLS"] >= 1:
        m = (disk_mean(m.astype(float), P["MASK_SMOOTH_CELLS"]) >= 0.5) & foot
    if hd > 0:                                                 # head recovery: the fan continues up to its true apex
        s_tan = math.tan(math.radians(max(ra["s_ref"], 0.5)))
        s_surf = horn_slope(np.where(np.isfinite(dw), dw, np.nanmean(dw)), grid.dx, grid.dy)
        head = ((dist <= hd) & np.isfinite(dw) & (dw > z_apex - P["CONE_TOLERANCE_M"]) &
                (dw <= z_apex + 2.0 * s_tan * dist + P["CONE_TOLERANCE_M"]) & (s_surf <= max(3.0 * ra["s_ref"], 12.0)))
        m |= disk_mean(head.astype(float), 1) >= 0.5
    sr, sc = r_a - r0, c_a - c0
    if not (0 <= sr < m.shape[0] and 0 <= sc < m.shape[1]):
        return None
    if not m[sr, sc]:                                          # snap seed to the nearest fan cell
        pts = np.argwhere(m)
        if pts.size == 0:
            return None
        j = np.argmin(((pts - [sr, sc]) ** 2).sum(axis=1))
        if ((pts[j] - [sr, sc]) ** 2).sum() > 16:
            return None
        sr, sc = int(pts[j][0]), int(pts[j][1])
    m = fill_holes(flood(m, [(sr, sc)]))
    m &= np.isfinite(dw)
    if not m.any():
        return None
    return dict(mask=m, cone=np.where(m, cone, np.nan), win=(r0, r1, c0, c1),
                area_km2=float(m.sum() * grid.dx * grid.dy / 1e6),
                mean_slope=float(np.nanmean(slope[r0:r1, c0:c1][m])))


# ---- apex refinement along the stem -------------------------------------------
def eval_apex(ctx, cell, P):
    """(radial analysis, fan) for an apex at a flat cell index; fan is None if no acceptable fan starts there."""
    W = ctx.grid.W
    ax, ay = ctx.grid.rc_to_xy(cell // W, cell % W)
    ra = radial_analysis(ctx.zsurf, ctx.grid, ax, ay, P)
    if ra is None or ra["s_ref"] < P["MIN_FAN_SLOPE_DEG"] or ra["n_fan"] < P["MIN_FAN_RAYS"]:
        return ra, None
    return ra, build_fan(ctx.zsurf, ctx.slope, ctx.grid, ax, ay, ra, P)


def refine_apex(ctx, cell, P):
    """The detectors only know where the terrain opens up; the fan itself knows where it starts. Walk UP the main
    channel in strides and keep moving the apex while every step still gains at least APEX_REFINE_GAIN of fan area
    (the signature of a lost fan head and flanks); stop after two steps that no longer pay off. At a sharp mountain
    front the area climbs steeply to the true apex; where the valley opens gradually it is flat and the apex stays."""
    step = ctx.step
    stride = max(1, int(round(mscale(P["APEX_REFINE_STEP_M"], 3, step) / step)))
    up = [cell]
    for _ in range(int(P["APEX_REFINE_UP_M"] / step)):
        n = ctx.main_donor_all[up[-1]]
        if n < 0:
            break
        up.append(n)
    ra, fan = eval_apex(ctx, cell, P)
    cur = (cell, ra, fan)
    prev = fan["area_km2"] if fan is not None else 0.0
    fails = 0
    for i in range(stride, len(up), stride):
        c = up[i]
        if not ctx.core_flat[c]:
            break
        ra2, fan2 = eval_apex(ctx, c, P)
        a = fan2["area_km2"] if fan2 is not None else 0.0
        if a > prev * (1.0 + P["APEX_REFINE_GAIN"]) and a > 0:
            cur, fails = (c, ra2, fan2), 0
        else:
            fails += 1
            if fails >= 2:
                break
        prev = a
    return cur


# ---- coalescence: fans -> fan complexes ----------------------------------------
def merge_fans(fans, ctx, P):
    """Group fans that overlap or touch into complexes (one polygon per bajada / flood plain) and attach the union of
    all the catchments feeding them. Every fan belongs to exactly one complex."""
    grid = ctx.grid
    n = len(fans)
    g = max(1, int(round(P["MERGE_GAP_M"] / min(grid.dx, grid.dy))))
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    if P["MERGE_FANS"]:
        for i in range(n):
            a = fans[i]["win"]
            for j in range(i + 1, n):
                b = fans[j]["win"]
                if find(i) == find(j) or a[1] + g < b[0] or b[1] + g < a[0] or a[3] + g < b[2] or b[3] + g < a[2]:
                    continue
                r0, r1, c0, c1 = min(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), max(a[3], b[3])
                mi = np.zeros((r1 - r0, c1 - c0), dtype=bool)
                mj = np.zeros_like(mi)
                mi[a[0] - r0:a[1] - r0, a[2] - c0:a[3] - c0] = fans[i]["mask"]
                mj[b[0] - r0:b[1] - r0, b[2] - c0:b[3] - c0] = fans[j]["mask"]
                ov = int((mi & mj).sum())
                linked = ov > 0 and ov >= P["MERGE_MIN_OVERLAP"] * min(int(mi.sum()), int(mj.sum()))
                if not linked:
                    linked = bool(((disk_mean(mi.astype(float), g) > 0) & mj).any())
                if linked:
                    parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    dem_flat = ctx.dem.ravel()
    out = []
    for members in groups.values():
        wins = [fans[i]["win"] for i in members]
        r0, r1 = min(w[0] for w in wins), max(w[1] for w in wins)
        c0, c1 = min(w[2] for w in wins), max(w[3] for w in wins)
        u = np.zeros((r1 - r0, c1 - c0), dtype=bool)
        for i in members:
            w = fans[i]["win"]
            u[w[0] - r0:w[1] - r0, w[2] - c0:w[3] - c0] |= fans[i]["mask"]
        if len(members) > 1:                                   # close the thin gaps between neighbouring fans
            up = np.pad(u, g)
            dil = disk_mean(up.astype(float), g) > 0
            closed = (disk_mean(dil.astype(float), g) >= 0.999)[g:-g, g:-g] | u
            u = fill_holes(closed & ctx.valid[r0:r1, c0:c1])
        cells = np.unique(np.concatenate([fans[i]["catch_cells"] for i in members]))
        zmax = float(np.nanmax(dem_flat[cells]))
        zapex = min(fans[i]["cand"]["z"] for i in members)
        area_c = cells.size * ctx.cell_km2
        flags = sorted(set(f for i in members for f in fans[i]["flags"].split(",") if f))
        out.append(dict(members=sorted(members), mask=u, win=(r0, r1, c0, c1), area_km2=float(u.sum() * ctx.cell_km2),
                        catch_cells=cells, catch_area_km2=area_c, catch_relief_m=zmax - zapex,
                        melton=(zmax - zapex) / 1000.0 / math.sqrt(max(area_c, 1e-9)),
                        catch_complete=not bool(ctx.border_flat[cells].any()), flags=",".join(flags), iou_ref=None))
    out.sort(key=lambda c: -c["area_km2"])
    for k, cx in enumerate(out):
        cx["cid"] = k + 1
        for i in cx["members"]:
            fans[i]["complex_id"] = k + 1
        cx["fan_ids"] = ",".join(str(fans[i]["fid"]) for i in cx["members"])
    return out


def score_reference(fans, complexes, reference, log):
    """IoU / recall of every reference polygon against the best-matching fan and the best-matching complex."""
    for ref in reference:
        rm, rarea = ref["mask"], int(ref["mask"].sum())
        if rarea == 0:
            log("  reference %s is empty on this DEM grid" % ref["name"])
            continue
        best = {}
        for kind, items in (("fan", fans), ("complex", complexes)):
            b = (0.0, None, 0.0)
            for it in items:
                r0, r1, c0, c1 = it["win"]
                inter = int((it["mask"] & rm[r0:r1, c0:c1]).sum())
                iou = inter / float(int(it["mask"].sum()) + rarea - inter)
                if it.get("iou_ref") is None or iou > it["iou_ref"]:
                    it["iou_ref"] = iou
                if iou > b[0]:
                    b = (iou, it, inter / float(rarea))
            best[kind] = b
        for kind in ("fan", "complex"):
            iou, it, rec = best[kind]
            if it is None:
                log("  reference %-8s (%.2f km2): no overlapping %s" % (ref["name"], rarea * ref["cell_km2"], kind))
            else:
                log("  reference %-8s (%.2f km2): best %-7s IoU %.2f | recall %3.0f%% | %s %.2f km2" % (
                    ref["name"], rarea * ref["cell_km2"], kind, iou, 100 * rec, kind, it["area_km2"]))


# ---- the whole algorithm ---------------------------------------------------------
def run_pipeline(dem, gt, P, manual=None, reference=None, log=log):
    """dem: 2-D float array (NaN = nodata) in a projected CRS (metres); gt: GDAL geotransform;
    manual: [dict(name, system, x, y)] apex coordinates in the DEM CRS; reference: [dict(name, mask)] hand-traced fans
    rasterised on the DEM grid. Returns a dict of results."""
    P = dict(P)
    H, W = dem.shape
    grid = Grid(gt, (H, W))
    step = min(grid.dx, grid.dy)
    cell_km2 = grid.dx * grid.dy / 1e6
    valid = np.isfinite(dem)
    method, mode = P["APEX_METHOD"].lower(), P["APEX_MODE"].lower()

    log("Hydrological analysis: sink filling + flow directions (%d x %d cells of %.1f m) ..." % (H, W, step))
    filled, recv, order = priority_flood(dem, valid)
    log("Flow accumulation ...")
    acc = flow_accumulation(recv, order)
    acc_km2 = P["FLOW_ACC_THRESHOLD_KM2"] or max(0.5, 400 * cell_km2)
    streams = (acc >= acc_km2 / cell_km2) & valid.ravel()
    log("  drainage network: %d channel cells (threshold %.0f cells = %.2f km2)" %
        (streams.sum(), acc_km2 / cell_km2, acc_km2))

    log("Fan surface (channels filled) and slope ...")
    zsurf = disk_mean(closing(dem, int(np.clip(round(P["CHANNEL_FILL_M"] / step), 1, 6))), 1)
    zf = fill_nan(dem)
    slope = horn_slope(zf, grid.dx, grid.dy)
    slope[~valid] = np.nan
    ctx = SimpleNamespace(dem=dem, valid=valid, grid=grid, step=step, cell_km2=cell_km2, recv=recv, acc=acc,
                          streams_flat=streams, zsurf=zsurf, slope=slope, core_flat=erode(valid, P["EDGE_BUFFER_CELLS"]).ravel(),
                          border_flat=border_cells(valid).ravel(), main_donor_all=build_main_donor(recv, acc, streams))
    ctx.starts, ctx.donors = build_donors(recv)
    ctx.order = order
    if method in ("paper", "both"):
        slope_s = disk_mean(slope, max(1, int(round(P["CONV_RADIUS_M"] / step))))
        thr1 = P["SLOPE_THR_1_DEG"] or otsu_threshold(slope_s[valid])
        log("  paper method: slope threshold 1 = %.1f deg (%s)" % (thr1, "given" if P["SLOPE_THR_1_DEG"] else "Otsu"))
        steep = np.nan_to_num(slope_s, nan=0.0) >= thr1
        curv = surface_curvature(zf, grid.dx, grid.dy, P["CURV_BASELINE_CELLS"])
        curv[~valid] = np.nan
        ctx.steep_flat, ctx.network_flat, ctx.curv_flat, ctx.slope_s = steep.ravel(), streams & steep.ravel(), curv.ravel(), slope_s
    else:
        steep, curv, slope_s = np.zeros((H, W), dtype=bool), np.full((H, W), np.nan), np.full((H, W), np.nan)

    man = snap_manual(manual or [], grid, streams, acc, P, log)
    known = [m for m in man if m["status"] == "candidate"]
    auto = []
    if mode in ("auto", "both"):
        if method in ("openness", "both"):
            a = detect_apexes_openness(ctx, P)
            log("Apex detection (valley mouths): %d found" % sum(1 for x in a if x["status"] == "candidate"))
            auto += a
        if method in ("paper", "both"):
            a = detect_apexes_paper(ctx, P)
            for x in a:
                x["method"] = "paper"
            log("Apex detection (paper steep-zone exits): %d exits, %d pass the curvature/length tests" %
                (len(a), sum(1 for x in a if x["status"] == "candidate")))
            auto += a
    dup_ref = known if mode in ("manual", "both") else []
    for a in auto:                                     # a detected apex next to a manual one is a duplicate
        for m in dup_ref:
            if a["status"] == "candidate" and math.hypot(a["row"] - m["row"], a["col"] - m["col"]) * step < P["MERGE_DIST_M"]:
                a["status"] = "rejected: duplicate of manual apex " + m["name"]
    for m in known:                                    # validation: how far is the detector from each known apex?
        near = [math.hypot(a["row"] - m["row"], a["col"] - m["col"]) * step for a in auto
                if a["status"] == "candidate" or a["status"].startswith("rejected: duplicate")]
        log("  known apex %-8s (catchment %.1f km2): nearest detected apex is %s" %
            (m["name"], acc[m["cell"]] * cell_km2, ("%.0f m away" % min(near)) if near else "none"))
    cands = (man if mode in ("manual", "both") else []) + auto

    # ---- iterative selection of fan apexes: manual first, then the biggest catchments --------------------
    dem_flat = dem.ravel()
    fans = []
    fan_id = np.zeros((H, W), dtype=np.int32)
    cone_r = np.full((H, W), np.nan, dtype=np.float32)
    inside_any = np.zeros((H, W), dtype=bool)
    in_catch = np.zeros(H * W, dtype=bool)                     # cells draining through an accepted apex
    todo = [c for c in cands if c["status"] == "candidate"]
    todo.sort(key=lambda c: (0 if c["source"] == "manual" else 1, -acc[c["cell"]], -dem_flat[c["cell"]]))
    log("Delineating fans for %d apexes ..." % len(todo))
    for n_done, cand in enumerate(todo):
        r, c = cand["row"], cand["col"]
        ax, ay = grid.rc_to_xy(r, c)
        cand["x"], cand["y"], cand["z"] = ax, ay, float(dem[r, c])
        cand["acc_km2"] = float(acc[cand["cell"]] * cell_km2)
        auto_c = cand["source"] == "auto"
        if auto_c and P["SUPPRESS_APEX_INSIDE_FAN"] and inside_any[r, c]:
            cand["status"] = "rejected: inside a larger fan"
            continue
        if auto_c and P["APEX_REFINE"]:
            cell2, ra, fan = refine_apex(ctx, cand["cell"], P)
            if cell2 != cand["cell"]:
                cand["cell"], r, c = cell2, cell2 // W, cell2 % W
                ax, ay = grid.rc_to_xy(r, c)
                cand["row"], cand["col"] = int(r), int(c)
                cand["x"], cand["y"], cand["z"] = ax, ay, float(dem[r, c])
                cand["acc_km2"] = float(acc[cell2] * cell_km2)
                if inside_any[r, c]:
                    cand["status"] = "rejected: inside a larger fan"
                    continue
            if fan is None:
                cand["status"] = "rejected: no fan found (radial profiles)"
                continue
        else:
            ra = radial_analysis(zsurf, grid, ax, ay, P)
            if ra is None:
                cand["status"] = "rejected: apex outside DEM"
                continue
            if ra["s_ref"] < P["MIN_FAN_SLOPE_DEG"]:
                cand["status"] = "rejected: terrain too flat for a fan (%.1f deg)" % ra["s_ref"]
                continue
            fan = build_fan(zsurf, slope, grid, ax, ay, ra, P) if ra["n_fan"] >= P["MIN_FAN_RAYS"] else None
            if fan is None:
                cand["status"] = "rejected: no fan found (radial profiles)"
                continue
        length = float(ra["D"].max())
        span = float((ra["D"] >= 0.25 * length).sum() * P["RADIAL_STEP_DEG"])
        if fan["area_km2"] < P["MIN_FAN_AREA_KM2"] or length < P["MIN_FAN_LENGTH_M"]:
            cand["status"] = "rejected: fan too small (%.3f km2, %.0f m)" % (fan["area_km2"], length)
            continue
        if auto_c and span < P["MIN_SPAN_DEG"]:
            cand["status"] = "rejected: corridor-shaped (span %.0f deg)" % span
            continue
        wr0, wr1, wc0, wc1 = fan["win"]
        ovl = float((fan["mask"] & inside_any[wr0:wr1, wc0:wc1]).sum()) / float(fan["mask"].sum())
        nested = bool(in_catch[cand["cell"]])                  # apex lies upstream of an accepted apex
        # ---- watershed analysis: catchment + contributing drainage ----
        ups = upstream_cells(cand["cell"], ctx.starts, ctx.donors)
        area_c = ups.size * cell_km2
        if auto_c and fan["area_km2"] > P["MAX_FAN_CATCHMENT_RATIO"] * area_c:
            cand["status"] = "rejected: fan (%.2f km2) larger than its catchment (%.2f km2)" % (fan["area_km2"], area_c)
            continue
        if auto_c and ovl >= (P["MAX_NESTED_OVERLAP"] if nested else P["MAX_FAN_OVERLAP"]):
            cand["status"] = "rejected: duplicate (%.0f%% inside earlier fans%s)" % (100 * ovl, ", nested apex" if nested else "")
            continue
        fid = len(fans) + 1
        cand["status"] = "accepted"
        cand["fan_id"] = fid
        in_catch[ups] = True
        inC = np.zeros(H * W, dtype=bool)
        inC[ups] = True
        sc = ups[streams[ups]]
        rc_ = recv[sc]
        k = (rc_ >= 0) & inC[np.maximum(rc_, 0)]
        sc, rc_ = sc[k], rc_[k]
        x1, y1 = grid.rc_to_xy(sc // W, sc % W)
        x2, y2 = grid.rc_to_xy(rc_ // W, rc_ % W)
        zmax = float(np.nanmax(dem_flat[ups]))
        complete = not bool(ctx.border_flat[ups].any())
        d = ra["D"] > 0
        flags = []
        if ra["n_trunc"] > 0:
            flags.append("max-length")
        if ra["n_edge"] > 0:
            flags.append("DEM-edge")
        if span < 60:
            flags.append("narrow")
        if not complete:
            flags.append("catchment-cut")
        fan.update(fid=fid, cand=cand, ra=ra, length_m=length, span_deg=span, catch_cells=ups,
                   catch_area_km2=area_c, catch_relief_m=zmax - cand["z"],
                   melton=(zmax - cand["z"]) / 1000.0 / math.sqrt(max(area_c, 1e-9)), catch_complete=complete,
                   drain_segments=(x1, y1, x2, y2), drain_len_km=float(np.hypot(x2 - x1, y2 - y1).sum() / 1000.0),
                   cone_slope_deg=float(np.degrees(np.arctan((ra["z_apex"] - ra["Zend"][d]) / ra["D"][d])).mean()),
                   toe_elev_min=float(ra["Zend"][d].min()), flags=",".join(flags), complex_id=0, iou_ref=None)
        sel = fan["mask"]
        inside_any[wr0:wr1, wc0:wc1] |= sel
        free = sel & (fan_id[wr0:wr1, wc0:wc1] == 0)
        fan_id[wr0:wr1, wc0:wc1][free] = fid
        cone_r[wr0:wr1, wc0:wc1][free] = fan["cone"][free]
        fans.append(fan)
        if (n_done + 1) % 20 == 0:
            log("  ... %d / %d apexes processed" % (n_done + 1, len(todo)))
    log("Done: %d fans accepted, %d candidates rejected" % (len(fans), len(cands) - len(fans)))

    complexes = merge_fans(fans, ctx, P)
    complex_id = np.zeros((H, W), dtype=np.int32)
    for cx in reversed(complexes):
        r0, r1, c0, c1 = cx["win"]
        complex_id[r0:r1, c0:c1][cx["mask"]] = cx["cid"]
    log("Coalescence: %d fans -> %d fan complexes" % (len(fans), len(complexes)))
    if reference:
        log("VALIDATION against %d reference polygon(s)" % len(reference))
        for ref in reference:
            ref["cell_km2"] = cell_km2
        score_reference(fans, complexes, reference, log)
    return dict(grid=grid, ctx=ctx, params=P, candidates=cands, fans=fans, complexes=complexes,
                acc_km2=(acc * cell_km2).reshape(H, W), streams=streams.reshape(H, W), slope=slope, slope_s=slope_s,
                steep=steep, feeder=(streams & steep.ravel()).reshape(H, W), curv=curv, zsurf=zsurf,
                fan_id=fan_id, cone=cone_r, complex_id=complex_id, filled=filled)



# =============================================================================
# 4. AHP RISK MODULE  (pure numpy - needs no QGIS, so it can be unit-tested)
# =============================================================================
RANDOM_INDEX = {1: 0.0, 2: 0.0, 3: 0.58, 4: 0.90, 5: 1.12, 6: 1.24, 7: 1.32, 8: 1.41, 9: 1.45, 10: 1.49}

# Default pairwise matrices (Saaty 1-9 scale; expert judgement). Replace them per site with [ahp] catchment_matrix / fan_matrix in site.toml.
CATCH_NAMES = ["slope", "local_relief", "upstream_area", "flow_accumulation", "hist_rainfall", "ruggedness", "bare_sparse"]
CATCH_MATRIX = np.array([[1, 2, 2, 3, 3, 4, 5],
                         [1 / 2, 1, 2, 2, 3, 3, 4],
                         [1 / 2, 1 / 2, 1, 2, 2, 3, 4],
                         [1 / 3, 1 / 2, 1 / 2, 1, 2, 2, 3],
                         [1 / 3, 1 / 3, 1 / 2, 1 / 2, 1, 2, 3],
                         [1 / 4, 1 / 3, 1 / 3, 1 / 2, 1 / 2, 1, 2],
                         [1 / 5, 1 / 4, 1 / 4, 1 / 3, 1 / 3, 1 / 2, 1]], dtype=float)
FAN_NAMES = ["hand_inv", "fan_slope", "flow_accumulation", "upstream_area", "local_relief", "event_rainfall"]
FAN_MATRIX = np.array([[1, 2, 2, 3, 3, 2],
                       [1 / 2, 1, 2, 2, 2, 2],
                       [1 / 2, 1 / 2, 1, 2, 2, 2],
                       [1 / 3, 1 / 2, 1 / 2, 1, 2, 2],
                       [1 / 3, 1 / 2, 1 / 2, 1 / 2, 1, 2],
                       [1 / 2, 1 / 2, 1 / 2, 1 / 2, 1 / 2, 1]], dtype=float)


def ahp_solve(A):
    """Principal-eigenvector weights, lambda_max, consistency index and consistency ratio (Saaty)."""
    A = np.asarray(A, dtype=float)
    n = len(A)
    if n < 2:
        return np.array([1.0]), 1.0, 0.0, 0.0
    w, v = np.linalg.eig(A)
    k = int(np.argmax(w.real))
    lam = float(w[k].real)
    wt = np.abs(v[:, k].real)
    wt /= wt.sum()
    ci = (lam - n) / (n - 1)
    ri = RANDOM_INDEX.get(n, 1.49)
    return wt, lam, ci, (ci / ri if ri > 0 else 0.0)


def ahp_weights(names, A, present):
    """Weights for the criteria in `present`: the sub-matrix of the full AHP matrix is re-solved (so a missing optional
    layer never leaves a hole), then flow accumulation and upstream area - the same quantity when both come from a
    DEM - are merged into one criterion that keeps their combined weight."""
    idx = [i for i, nme in enumerate(names) if nme in present]
    w, lam, ci, cr = ahp_solve(np.asarray(A)[np.ix_(idx, idx)])
    wd = {names[i]: float(w[j]) for j, i in enumerate(idx)}
    if "flow_accumulation" in wd and "upstream_area" in wd:
        wd["upstream_area"] += wd.pop("flow_accumulation")
    return wd, dict(cr=float(cr), n=len(idx), lam=float(lam), dropped=[n_ for n_ in names if n_ not in present])


# ---- DEM-derived criteria ------------------------------------------------------
def focal_range(z, r):
    """Local relief: max - min inside a (2r+1) square window (separable, NaN-aware)."""
    r = int(r)
    H, W = z.shape
    fin = np.isfinite(z)

    def sweep(a, fn, pad):
        out = a.copy()
        p = np.pad(a, ((0, 0), (r, r)), constant_values=pad)
        for d in range(-r, r + 1):
            out = fn(out, p[:, r + d:r + d + W])
        res = out.copy()
        p = np.pad(out, ((r, r), (0, 0)), constant_values=pad)
        for d in range(-r, r + 1):
            res = fn(res, p[r + d:r + d + H, :])
        return res

    hi = sweep(np.where(fin, z, -np.inf), np.maximum, -np.inf)
    lo = sweep(np.where(fin, z, np.inf), np.minimum, np.inf)
    out = hi - lo
    out[~fin | ~np.isfinite(out)] = np.nan
    return out


def riley_tri(z, valid):
    """Terrain Ruggedness Index (Riley et al. 1999): sqrt of the summed squared differences to the 8 neighbours (m)."""
    zf = fill_nan(z)
    zp = np.pad(zf, 1, mode="edge")
    H, W = zf.shape
    s = np.zeros((H, W))
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr or dc:
                s += (zp[1 + dr:1 + dr + H, 1 + dc:1 + dc + W] - zf) ** 2
    out = np.sqrt(s)
    out[~valid] = np.nan
    return out


def compute_hand(dem_flat, recv, stream_flat, valid_flat):
    """Height Above Nearest Drainage: elevation above the first channel cell reached by following the flow path
    (pointer jumping, so it is vectorised). NaN where the flow path never reaches a channel."""
    n = recv.size
    nxt = np.where(stream_flat | (recv < 0), np.arange(n), recv)
    for _ in range(64):
        n2 = nxt[nxt]
        if np.array_equal(n2, nxt):
            break
        nxt = n2
    ok = stream_flat[nxt] & valid_flat
    return np.where(ok, np.maximum(dem_flat - dem_flat[nxt], 0.0), np.nan)


def pct_range(v, lo=2.0, hi=98.0):
    """(low, high) percentile pair used for p2-p98 min-max normalisation (the paper's rule)."""
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0, 1.0
    a, b = np.percentile(v, [lo, hi])
    return float(a), float(b if b > a else a + 1e-9)


def norm01(v, rng):
    return np.clip((v - rng[0]) / (rng[1] - rng[0]), 0.0, 1.0)


def log_area_scale(acc_km2, amin, aref):
    """Upstream area on a FIXED log scale (portable between DEMs): 0 at amin, 1 at aref km2."""
    a = np.maximum(np.asarray(acc_km2, dtype=float), amin)
    return np.clip((np.log10(a) - math.log10(amin)) / (math.log10(aref) - math.log10(amin)), 0.0, 1.0)


# ---- Strahler order, perimeters (the Stage 1 columns that were "unavailable") -----------------------------------
def strahler(recv, order, stream_flat):
    """Strahler order of every channel cell (0 elsewhere)."""
    n = recv.size
    mxl, cntl, sol = [0] * n, [0] * n, [0] * n
    rl = recv.tolist()
    for c in order[stream_flat[order]][::-1].tolist():          # upstream first
        o = mxl[c]
        o = 1 if o == 0 else (o + 1 if cntl[c] >= 2 else o)
        sol[c] = o
        r = rl[c]
        if r >= 0 and stream_flat[r]:
            if o > mxl[r]:
                mxl[r], cntl[r] = o, 1
            elif o == mxl[r]:
                cntl[r] += 1
    return np.array(sol, dtype=np.int32)


def horton_ratios(so, recv, in_catch, stream_flat, step):
    """(max order, mean bifurcation ratio, mean stream-length ratio) inside one catchment."""
    cells = np.flatnonzero(in_catch & stream_flat)
    if cells.size == 0:
        return 0, float("nan"), float("nan")
    umax = int(so[cells].max())
    if umax < 2:
        return umax, float("nan"), float("nan")
    r = recv[cells]
    rr = np.maximum(r, 0)
    ends = (r < 0) | ~stream_flat[rr] | ~in_catch[rr] | (so[rr] > so[cells])
    N = np.array([int((ends & (so[cells] == u)).sum()) for u in range(1, umax + 1)], dtype=float)
    Ltot = np.array([float((so[cells] == u).sum()) * step for u in range(1, umax + 1)])
    rb = [N[u] / N[u + 1] for u in range(umax - 1) if N[u + 1] > 0 and N[u] > 0]
    lbar = Ltot / np.maximum(N, 1)
    rl = [lbar[u] / lbar[u - 1] for u in range(1, umax) if lbar[u - 1] > 0 and N[u] > 0]
    return umax, (float(np.mean(rb)) if rb else float("nan")), (float(np.mean(rl)) if rl else float("nan"))


def raster_perimeter_km(mask, step, smooth=False):
    """Perimeter of a raster region. smooth=False counts cell faces, which is what a polygon vectorised from a raster
    (a catchment) has; smooth=True applies the pi/4 correction for hand-traced / smoothed outlines (fans)."""
    m = np.pad(mask, 1)
    edges = (m[1:, :] != m[:-1, :]).sum() + (m[:, 1:] != m[:, :-1]).sum()
    return float(edges) * step * (math.pi / 4.0 if smooth else 1.0) / 1000.0


# ---- the AHP driver ---------------------------------------------------------------------------------
def build_crit(res, P, layers, log):
    """Region-wide criteria rasters (DEM-derived) and the p2-p98 ranges used by the 'region' normalisation."""
    ctx, grid = res["ctx"], res["grid"]
    H, W = ctx.dem.shape
    step, valid = ctx.step, ctx.valid
    c = SimpleNamespace()
    c.H, c.W, c.step, c.valid = H, W, step, valid
    c.slope = ctx.slope
    c.fan_slope = horn_slope(fill_nan(res["zsurf"]), grid.dx, grid.dy)           # slope of the channel-filled fan surface
    c.fan_slope[~valid] = np.nan
    rad = max(1, int(round(mscale(P["LOCAL_RELIEF_RADIUS_M"], 3, step) / step)))
    c.relief = focal_range(ctx.dem, rad)
    c.tri = riley_tri(ctx.dem, valid)
    c.acc_km2 = ctx.acc.reshape(H, W) * ctx.cell_km2
    c.rng = {k: pct_range(getattr(c, k)[valid]) for k in ("slope", "fan_slope", "relief", "tri")}
    c.hand_cache = {}
    c.layers = layers
    for k in ("hist", "bare", "built"):
        v = layers.get(k)
        if v is not None:
            c.rng[k] = pct_range(v[valid]) if np.nanmax(v) > 1.001 else (0.0, 1.0)
    rain = layers.get("rain")
    c.rain_factor = None if rain is None else np.clip(rain / P["RAINFALL_REF_MM"], 0.0, 1.0)
    log("Criteria rasters: slope, fan slope, local relief (+-%d cells), ruggedness (TRI), upstream area%s" % (
        rad, "".join(", " + n for n, k in (("event rainfall", "rain"), ("historical rainfall", "hist"), ("bare/sparse", "bare"),
                                            ("built density (exposure)", "built")) if layers.get(k) is not None)))
    return c


def get_hand(crit, ctx, thr_km2):
    """Local HAND (m) for a channel-initiation threshold, plus its region p2-p98 range (cached)."""
    key = round(float(thr_km2), 4)
    if key not in crit.hand_cache:
        st = (ctx.acc * ctx.cell_km2 >= thr_km2) & crit.valid.ravel()
        h = compute_hand(ctx.dem.ravel(), ctx.recv, st, crit.valid.ravel()).reshape(crit.H, crit.W)
        crit.hand_cache[key] = (h, pct_range(h[np.isfinite(h)]))
    return crit.hand_cache[key]


def catchment_metrics(ctx, crit, so, st_m, pos, ups, outlet, P, rain_ups=None):
    """Stage-1 master-table metrics of one catchment (cells `ups`, drained through `outlet`)."""
    H, W = crit.H, crit.W
    step, ck = ctx.step, ctx.cell_km2
    dem_f = ctx.dem.ravel()
    area = ups.size * ck
    zmin = float(dem_f[outlet])
    zmax = float(np.nanmax(dem_f[ups]))
    rows, cols = ups // W, ups % W
    r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
    m = np.zeros((r1 - r0, c1 - c0), dtype=bool)
    m[rows - r0, cols - c0] = True
    in_c = np.zeros(H * W, dtype=bool)
    in_c[ups] = True
    sc = ups[st_m[ups]]                                                # channel length inside the catchment
    rc = ctx.recv[sc]
    k = (rc >= 0) & in_c[np.maximum(rc, 0)]
    sc, rc = sc[k], rc[k]
    diag = (np.abs(sc // W - rc // W) + np.abs(sc % W - rc % W)) == 2
    stream_km = float((step * np.where(diag, math.sqrt(2.0), 1.0)).sum() / 1000.0)
    fd = {int(outlet): 0.0}                                            # flow distance to the outlet, downstream first
    rl = ctx.recv
    for c in ups[np.argsort(pos[ups])].tolist():
        if c == outlet:
            continue
        r = int(rl[c])
        if r in fd:
            fd[c] = fd[r] + step * (math.sqrt(2.0) if (abs(c // W - r // W) + abs(c % W - r % W)) == 2 else 1.0)
    far = max(fd, key=fd.get)
    Lm = fd[far] / 1000.0
    umax, rb, rlr = horton_ratios(so, ctx.recv, in_c, st_m, step)
    per = raster_perimeter_km(m, step)
    relief = zmax - zmin
    tc = float("nan")
    if Lm > 0 and dem_f[far] > zmin:
        S = (dem_f[far] - zmin) / (Lm * 1000.0)
        tc = 0.0195 * (Lm * 1000.0) ** 0.77 * S ** -0.385                # Kirpich (1940), minutes
    tc_act = tc
    if rain_ups is not None and np.isfinite(rain_ups).any() and np.nanmax(rain_ups) > 0:
        act = ups[np.nan_to_num(rain_ups, nan=0.0) >= P["RAIN_ACTIVE_FRAC"] * np.nanmax(rain_ups)]
        cand = [c for c in act.tolist() if c in fd]
        if cand:                                                      # response time from the farthest cell that actually gets rain
            fa = max(cand, key=fd.get)
            La = fd[fa]
            if La > 0 and dem_f[fa] > zmin:
                tc_act = 0.0195 * La ** 0.77 * ((dem_f[fa] - zmin) / La) ** -0.385
    return dict(area_km2=area, perimeter_km=per, min_elev_m=zmin, max_elev_m=zmax, relief_m=relief,
                mean_slope_deg=float(np.nanmean(crit.slope.ravel()[ups])), stream_length_km=stream_km,
                drainage_density=stream_km / area, melton=relief / 1000.0 / math.sqrt(area),
                melton_m_per_sqrt_km2=relief / math.sqrt(area), main_channel_km=Lm,
                main_channel_gradient_m_per_km=(float(dem_f[far]) - zmin) / Lm if Lm > 0 else float("nan"),
                form_factor=area / (Lm * Lm) if Lm > 0 else float("nan"), perimeter_form_factor_proxy=16.0 * area / (per * per),
                max_stream_order=umax, bifurcation_ratio=rb, stream_length_ratio=rlr, kirpich_tc_min=tc, tc_rain_active_min=tc_act)


def catchment_score(crit, ups, area_km2, wts, P):
    """C_i: catchment source susceptibility. Spatial criteria are the mean of the region-normalised rasters over the
    catchment; the upstream-area criterion is the catchment area on a fixed log scale."""
    v = {}
    for name, key, arr in (("slope", "slope", crit.slope), ("local_relief", "relief", crit.relief), ("ruggedness", "tri", crit.tri)):
        v[name] = float(np.nanmean(norm01(arr.ravel()[ups], crit.rng[key])))
    v["upstream_area"] = float(log_area_scale(area_km2, P["CATCH_AREA_MIN_KM2"], P["CATCH_AREA_REF_KM2"]))
    for name, key in (("hist_rainfall", "hist"), ("bare_sparse", "bare")):
        if name in wts:
            v[name] = float(np.nanmean(norm01(crit.layers[key].ravel()[ups], crit.rng[key])))
    return float(sum(wts[k] * v[k] for k in wts)), v


def feeder_inventory(ctx, full_mask, members, extra, P):
    """The systems that feed a fan complex (A1, A2, A3 ... in the paper). Anchored at the fan apexes (member fans, biggest
    catchment first), then any other stream that enters the complex (>= FEEDER_MIN_KM2), then user pour points. Nested
    catchments are skipped, so every feeder is a distinct source."""
    grid = ctx.grid
    H, W = ctx.dem.shape
    inside = full_mask.ravel()
    out, taken = [], np.zeros(H * W, dtype=bool)

    def add(cell, name, source, fan_id=0, ups=None):
        x, y = grid.rc_to_xy(cell // W, cell % W)
        ups = upstream_cells(cell, ctx.starts, ctx.donors) if ups is None else ups
        taken[ups] = True
        out.append(dict(name=name, cell=int(cell), x=x, y=y, acc_km2=float(ctx.acc[cell] * ctx.cell_km2), ups=ups,
                        source=source, fan_id=fan_id))

    for fan in sorted(members, key=lambda f: -f["catch_area_km2"]):             # (a) the fan apexes
        c = fan["cand"]["cell"]
        if len(out) >= P["MAX_FEEDERS"] or taken[c]:
            continue
        if fan["catch_area_km2"] >= P["FEEDER_MIN_KM2"] or not out:
            add(c, "F%d" % (len(out) + 1), "fan apex", fan["fid"], fan["catch_cells"])
    ent = ctx.streams_flat & ~inside & (ctx.recv >= 0)                           # (b) other streams entering the complex
    ent &= inside[np.maximum(ctx.recv, 0)]
    cells = np.flatnonzero(ent & (ctx.acc * ctx.cell_km2 >= P["FEEDER_MIN_KM2"]))
    for c in cells[np.argsort(-ctx.acc[cells])].tolist():
        if len(out) >= P["MAX_FEEDERS"]:
            break
        x, y = grid.rc_to_xy(c // W, c % W)
        if taken[c] or any(math.hypot(x - f["x"], y - f["y"]) < P["FEEDER_MIN_SEP_M"] for f in out):
            continue
        ups = upstream_cells(c, ctx.starts, ctx.donors)
        if np.isin(np.array([f["cell"] for f in out], dtype=np.int64), ups).any():   # would contain an existing feeder
            continue
        add(c, "F%d" % (len(out) + 1), "entering stream", 0, ups)
    for e in extra or []:                                                        # (c) user pour points
        add(e["cell"], e["name"], "user pour point")
    for f in out:                                                                # flow direction (HEC-RAS inflow line)
        r = int(ctx.recv[f["cell"]])
        if r >= 0:
            xr, yr = grid.rc_to_xy(r // W, r % W)
            f["dir"] = (xr - f["x"], yr - f["y"])
        else:
            f["dir"] = (0.0, -1.0)
    return out


def hazard_stack(crit, ctx, win, mask, wf, hd, scope, P, ua=None, rngs=None):
    """Fan hazard H(x) in one window: weighted sum of the normalised criteria.
    scope "region": ABSOLUTE physical scales (HAND_REF_M, FAN_SLOPE_REF_DEG, RELIEF_REF_M, FLOW_AREA_MIN/REF_KM2) - the same
                    number means the same thing on every fan and every DEM.
    scope "fan"   : relative, p2-p98 inside the fan (or the ranges `rngs` of a reference run) - shows structure inside one fan.
    hd: HAND window (m) from hand_window(). ua: contributing area (km2) of each cell; default = the single-path D8 area (v3)."""
    r0, r1, c0, c1 = win
    sl, rl = crit.fan_slope[r0:r1, c0:c1], crit.relief[r0:r1, c0:c1]
    ac = crit.acc_km2[r0:r1, c0:c1] if ua is None else ua
    got = {}
    if scope == "region":
        amin, aref = P["FLOW_AREA_MIN_KM2"], P["FLOW_AREA_REF_KM2"]
        n = {"hand_inv": np.where(np.isfinite(hd), np.clip(1.0 - hd / P["HAND_REF_M"], 0.0, 1.0), 0.0),
             "fan_slope": np.clip(sl / P["FAN_SLOPE_REF_DEG"], 0.0, 1.0),
             "local_relief": np.clip(rl / P["RELIEF_REF_M"], 0.0, 1.0),
             "upstream_area": np.clip(np.log10(np.maximum(ac, 1e-9) / amin) / math.log10(aref / amin), 0.0, 1.0)}
    else:
        R = rngs or {}
        la = np.log10(np.maximum(ac, 1e-6))
        rh = R.get("hand") or pct_range(hd[mask & np.isfinite(hd)])
        rs_ = R.get("slope") or pct_range(sl[mask & np.isfinite(sl)])
        rr_ = R.get("relief") or pct_range(rl[mask & np.isfinite(rl)])
        ra = R.get("area") or pct_range(la[mask & (ac > 1e-6)])
        got = dict(hand=rh, slope=rs_, relief=rr_, area=ra)
        n = {"hand_inv": np.where(np.isfinite(hd), 1.0 - norm01(hd, rh), 0.0), "fan_slope": norm01(sl, rs_),
             "local_relief": norm01(rl, rr_), "upstream_area": norm01(la, ra)}
    if "event_rainfall" in wf:
        n["event_rainfall"] = np.nan_to_num(crit.rain_factor[r0:r1, c0:c1], nan=0.0)
    Hh = np.where(mask, sum(wf[k] * n[k] for k in wf), np.nan)
    rs = int(round(mscale(P["RISK_SMOOTH_M"], 1, ctx.step) / ctx.step)) if P["RISK_SMOOTH_M"] > 0 else 0
    if rs >= 1:                                          # a flood path is wider than one 27 m cell: soften the channel lines
        Hh = np.where(mask, disk_mean(Hh, rs), np.nan)
    n["_rngs"] = got
    return Hh, n


def hand_window(crit, ctx, thr_km2, win, ua=None):
    """HAND (m) inside a window. Channels = the single-path D8 streams (>= thr_km2) plus, where the flow has been spread over a
    fan (ua), every cell that carries that much routed area (scaled with the cell size so the channel width does not change with
    the DEM resolution). A fan's flood paths are in its interior, not only on the one D8 line that happens to hug its margin."""
    r0, r1, c0, c1 = win
    st = ((ctx.acc * ctx.cell_km2) >= thr_km2) & crit.valid.ravel()
    if ua is not None:
        st = st.reshape(crit.H, crit.W).copy()
        st[r0:r1, c0:c1] |= ua >= thr_km2 * (ctx.step / 30.0)
        st = st.ravel()
    h = compute_hand(ctx.dem.ravel(), ctx.recv, st, crit.valid.ravel()).reshape(crit.H, crit.W)
    return h[r0:r1, c0:c1]


# ---- how a flood spreads over a fan: multiple-flow-direction routing -----------------------------
def mfd_prepare(z, valid, exponent, dx, dy):
    """Holmgren (1994) multiple-flow-direction weights, w ~ tan(slope)**exponent, on a pit-filled copy of z.
    A single D8 line (v3's 'upstream area') is one cell wide; a flood leaving a fan apex is not."""
    H, W = z.shape
    ok = np.isfinite(z) & valid
    filled, recv, order = priority_flood(np.where(ok, z, np.nan), ok)
    pos = np.zeros(H * W, dtype=np.int64)
    pos[order] = np.arange(order.size)
    f = np.where(ok, filled, -np.inf).ravel()
    proc = np.lexsort((-pos, -f))                        # high -> low; on flats the cells far from the outlet go first
    rank = np.empty(H * W, dtype=np.int64)
    rank[proc] = np.arange(H * W)
    rank2 = rank.reshape(H, W)
    idx = np.arange(H * W).reshape(H, W)
    Wt = np.zeros((H, W, 8))
    NB = np.full((H, W, 8), -1, dtype=np.int64)
    offs = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
    with np.errstate(invalid="ignore"):
        for k, (di, dj) in enumerate(offs):
            d = math.hypot(di * dy, dj * dx)
            a = (slice(max(0, -di), H - max(0, di)), slice(max(0, -dj), W - max(0, dj)))
            b = (slice(max(0, -di) + di, H - max(0, di) + di), slice(max(0, -dj) + dj, W - max(0, dj) + dj))
            down = ok[a] & ok[b] & (rank2[b] > rank2[a])
            s = np.maximum((filled[a] - filled[b]) / d, 1e-4)
            Wt[a[0], a[1], k] = np.where(down, s ** exponent, 0.0)
            NB[a[0], a[1], k] = np.where(down, idx[b], -1)
    tot = Wt.sum(2)
    tot[tot == 0] = 1.0
    Wt /= tot[:, :, None]
    return dict(shape=(H, W), proc=proc, nb=NB.reshape(-1, 8).astype(np.int32), w=Wt.reshape(-1, 8))


def mfd_accumulate(prep, source):
    """Route `source` (any quantity per cell) downhill; returns the amount passing through every cell (source included)."""
    H, W = prep["shape"]
    flux = np.asarray(source, dtype=float).ravel().tolist()
    nb, w = prep["nb"], prep["w"]
    for c in prep["proc"].tolist():
        v = flux[c]
        if v == 0.0:
            continue
        n, ww = nb[c].tolist(), w[c].tolist()
        for k in range(8):
            j = n[k]
            if j >= 0:
                flux[j] += v * ww[k]
    return np.array(flux).reshape(H, W)


def route_complex(res, ctx, cx, feeders, P):
    """Spread every feeder's catchment, and the rain falling on the fan's own hillsides, over the channel-filled surface.
    Returns dict(win = routing window (r0, r1, c0, c1), total = contributing area (km2) over that window,
    feed = the part that comes from the feeders only). Also stores f['flow_km2'] (that feeder alone, same window).
    The window reaches RUNOUT_MAX_M beyond the fan so the flood can be followed past the fan toe."""
    H, W = ctx.dem.shape
    r0, r1, c0, c1 = cx["win"]
    b_loc = int(math.ceil(P["ROUTE_MARGIN_M"] / ctx.step))
    b = int(math.ceil((P["ROUTE_MARGIN_M"] + max(0.0, P["RUNOUT_MAX_M"])) / ctx.step))
    R0, R1, C0, C1 = max(0, r0 - b), min(H, r1 + b), max(0, c0 - b), min(W, c1 + b)
    zf, zr = res["zsurf"][R0:R1, C0:C1], ctx.dem[R0:R1, C0:C1]
    bl = float(P.get("ROUTE_SURFACE_BLEND", 1.0))                      # 0 = route on the real DEM (keeps incised channels), 1 = channel-filled
    zroute = zf if bl >= 1.0 else np.where(np.isfinite(zr) & np.isfinite(zf), zr + bl * (zf - zr), zf)
    prep = mfd_prepare(zroute, ctx.valid[R0:R1, C0:C1], P["ROUTE_EXPONENT"], ctx.grid.dx, ctx.grid.dy)
    h, w = R1 - R0, C1 - C0
    own = np.zeros((h, w))                                            # local runoff: only inside the original margin
    L0, L1, K0, K1 = max(0, r0 - b_loc), min(H, r1 + b_loc), max(0, c0 - b_loc), min(W, c1 + b_loc)
    own[L0 - R0:L1 - R0, K0 - C0:K1 - C0] = np.where(ctx.valid[L0:L1, K0:K1], ctx.cell_km2, 0.0)
    for f in feeders:                                                 # ... except cells already inside a feeder catchment
        rr, cc = f["ups"] // W - R0, f["ups"] % W - C0
        k = (rr >= 0) & (rr < h) & (cc >= 0) & (cc < w)
        own[rr[k], cc[k]] = 0.0
    total = own.copy()
    feed = np.zeros((h, w))
    for f in feeders:
        r, c = f["cell"] // W - R0, f["cell"] % W - C0
        if 0 <= r < h and 0 <= c < w:
            src = np.zeros((h, w))
            src[r, c] = f["metrics"]["area_km2"]
            f["flow_km2"] = mfd_accumulate(prep, src)
            if not any(o is not f and o["ups"].size > f["ups"].size and bool((o["ups"] == f["cell"]).any()) for o in feeders):   # nested: counted once
                total[r, c] += f["metrics"]["area_km2"]
                feed += f["flow_km2"]
        else:
            f["flow_km2"] = None
    return dict(win=(R0, R1, C0, C1), total=mfd_accumulate(prep, total), feed=feed)


def _grow(seed, allowed, n_iter):
    """8-neighbour region growing of `seed` into `allowed`, at most n_iter cells outward."""
    g_ = seed.copy()
    allowed = allowed & ~seed
    for _ in range(int(n_iter)):
        d = g_.copy()
        d[1:, :] |= g_[:-1, :]
        d[:-1, :] |= g_[1:, :]
        d[:, 1:] |= g_[:, :-1]
        d[:, :-1] |= g_[:, 1:]
        d[1:, 1:] |= g_[:-1, :-1]
        d[:-1, :-1] |= g_[1:, 1:]
        d[1:, :-1] |= g_[:-1, 1:]
        d[:-1, 1:] |= g_[1:, :-1]
        new = d & allowed & ~g_
        if not new.any():
            break
        g_ |= new
    return g_


def _grow_oct(seed, allowed, n_iter):
    """Region growing that alternates 8- and 4-neighbour steps: an octagonal approximation of a circle of n_iter cells."""
    g_ = seed.copy()
    allowed = allowed & ~seed
    for i in range(int(n_iter)):
        d = g_.copy()
        d[1:, :] |= g_[:-1, :]
        d[:-1, :] |= g_[1:, :]
        d[:, 1:] |= g_[:, :-1]
        d[:, :-1] |= g_[:, 1:]
        if i % 2 == 0:
            d[1:, 1:] |= g_[:-1, :-1]
            d[:-1, :-1] |= g_[1:, 1:]
            d[1:, :-1] |= g_[:-1, 1:]
            d[:-1, 1:] |= g_[1:, :-1]
        new = d & allowed & ~g_
        if not new.any():
            break
        g_ |= new
    return g_


def runout_corridor(ctx, rw, fan_rw, feed, feeders, a_max, a_total, hand, P):
    """Cells beyond the fan that the feeders' flood still reaches. Two parts, both stopping at a river whose D8 area exceeds
    RUNOUT_STOP_RATIO x the feeders' total (e.g. the Indus):
      (a) the spread flow: connected to the fan, carrying >= RUNOUT_MIN_FRACTION of the largest feeder's area per cell;
      (b) the channel: the D8 flood path of every feeder, followed RUNOUT_MAX_M past the fan toe and widened by
          RUNOUT_HALF_WIDTH_M, keeping only cells within HAND_REF_M of the channel (a flood cannot climb a hillside).
    Part (b) matters because the incised channel can run along the fan edge, outside the fan polygon, while the spread flow (a) follows
    the fan surface. All arrays are in the routing window rw (hand = HAND of that window)."""
    R0, R1, C0, C1 = rw
    H, W = ctx.dem.shape
    acc_km2 = (ctx.acc * ctx.cell_km2).reshape(H, W)[R0:R1, C0:C1]
    free = ctx.valid[R0:R1, C0:C1] & (acc_km2 < P["RUNOUT_STOP_RATIO"] * a_total)
    thr = max(P["FLOW_AREA_MIN_KM2"], P["RUNOUT_MIN_FRACTION"] * a_max)
    cor = _grow(fan_rw, free & (feed >= thr), math.ceil(P["RUNOUT_MAX_M"] / ctx.step))
    path = np.zeros_like(fan_rw)
    h, w = fan_rw.shape
    for f in feeders:
        c, cells, dist, d_ = f["cell"], [], [], 0.0
        while c >= 0 and len(cells) < 50000:
            r, k = divmod(int(c), W)
            if not (R0 <= r < R1 and C0 <= k < C1):
                break
            if cells:
                d_ += ctx.step * (math.sqrt(2.0) if (r - R0 != cells[-1][0] and k - C0 != cells[-1][1]) else 1.0)
            cells.append((r - R0, k - C0))
            dist.append(d_)
            c = int(ctx.recv[c])
        in_fan = [i for i, (r, k) in enumerate(cells) if fan_rw[r, k]]
        if not in_fan:
            continue
        limit = dist[in_fan[-1]] + P["RUNOUT_MAX_M"]
        for (r, k), dd in zip(cells, dist):
            if dd > limit or not free[r, k]:
                break
            path[r, k] = True
    if path.any() and P["RUNOUT_HALF_WIDTH_M"] > 0:
        low = free & np.isfinite(hand) & (hand <= P["HAND_REF_M"])
        band = _grow_oct(path, low | path, math.ceil(P["RUNOUT_HALF_WIDTH_M"] / ctx.step))
        cor |= band
    return cor & ~fan_rw


# ---- event hydrology: SCS triangular hydrograph ---------------------------------------------------
def cn_for_amc(cn, amc):
    """Curve number for antecedent moisture condition 1 (dry) / 3 (wet) from the AMC II value (Chow, Maidment & Mays 1988)."""
    cn = np.asarray(cn, dtype=np.float64)
    if int(amc) == 1:
        return 4.2 * cn / (10.0 - 0.058 * cn)
    if int(amc) == 3:
        return 23.0 * cn / (10.0 + 0.13 * cn)
    return cn


def runoff_depth_mm(rain_mm, P, cn=None):
    """-> (runoff depth in mm, mean curve number used or None) for an event of `rain_mm`.
    RUNOFF_METHOD 'coefficient': RUNOFF_COEFF x rain.
    RUNOFF_METHOD 'scs_cn': Q = (P - Ia)^2 / (P - Ia + S), S = 25400 / CN - 254 mm, Ia = CN_IA_RATIO x S. `cn` may be an array (the curve
    numbers of the catchment's cells): the runoff is computed per cell and averaged, which is not the same as using the mean CN."""
    rain = max(0.0, float(rain_mm))
    if str(P.get("RUNOFF_METHOD", "coefficient")).lower() != "scs_cn":
        return P["RUNOFF_COEFF"] * rain, None
    if cn is None:
        cn = P.get("CURVE_NUMBER")
    if cn is None:
        raise ValueError("RUNOFF_METHOD is 'scs_cn' but no curve number was given (CURVE_NUMBER or a CN raster)")
    cn = np.atleast_1d(np.asarray(cn, dtype=np.float64))
    cn = cn[np.isfinite(cn) & (cn > 0)]
    if cn.size == 0:
        raise ValueError("the curve-number raster has no valid cell in this catchment")
    cn = np.clip(cn_for_amc(np.clip(cn, 1.0, 100.0), P.get("CN_AMC", 2)), 1.0, 100.0)
    S = 25400.0 / cn - 254.0
    ia = float(P.get("CN_IA_RATIO", 0.2)) * S
    q = np.where(rain > ia, (rain - ia) ** 2 / np.maximum(rain - ia + S, 1e-9), 0.0)
    return float(q.mean()), float(cn.mean())


def event_hydrology(area_km2, tc_min, rain_mm, P, cn=None):
    """Peak discharge, volume and time to peak of one feeder for an event of `rain_mm`.
    Qp = 0.208 A Q / Tp (SCS, peak-rate factor 484), Tp = D/2 + 0.6 Tc, Q = runoff_depth_mm(rain)."""
    q_mm, cn_used = runoff_depth_mm(rain_mm, P, cn)
    q_mm = max(0.0, q_mm)
    tc = tc_min if (tc_min == tc_min and tc_min > 0) else 0.0
    tp = P["STORM_DURATION_MIN"] / 2.0 + 0.6 * tc
    qp = 0.208 * area_km2 * q_mm / (tp / 60.0) if tp > 0 else 0.0
    return dict(runoff_mm=q_mm, tp_min=tp, peak_m3s=qp, volume_m3=q_mm / 1000.0 * area_km2 * 1.0e6, curve_number=cn_used)


def trigger_from_peak(qp, P):
    """Cloudburst trigger factor 0-1: peak discharge on a log scale from TRIGGER_Q_LOW to TRIGGER_Q_HIGH m3/s."""
    lo, hi = P["TRIGGER_Q_LOW_M3S"], P["TRIGGER_Q_HIGH_M3S"]
    return float(np.clip(math.log10(qp / lo) / math.log10(hi / lo), 0.0, 1.0)) if qp > lo else 0.0


# ---- validation against a flood footprint / flood points (success-rate curve after Sar et al. 2025) ----
def presence_auc(score, presence, domain, nbins=100):
    """Cells are ranked from high to low susceptibility and cut into nbins equal-area classes; x = cumulative share of the
    area, y = cumulative share of the flooded cells; AUC = area under y(x). 0.5 = no skill, 1 = perfect."""
    ok = domain & np.isfinite(score)
    s, p = score[ok], presence[ok].astype(float)
    if p.sum() == 0 or s.size < nbins:
        return float("nan")
    o = np.argsort(-s, kind="stable")
    cx = np.arange(1, s.size + 1) / s.size
    cy = np.cumsum(p[o]) / p.sum()
    xb = np.linspace(0.0, 1.0, nbins + 1)
    yb = np.interp(xb, np.concatenate([[0.0], cx]), np.concatenate([[0.0], cy]))
    return float(np.sum((yb[1:] + yb[:-1]) / 2.0 * np.diff(xb)))


def roc_auc(score, presence, domain):
    """Standard ROC AUC (flooded vs not-flooded cells, tie-aware rank-sum): 1 = perfect, 0.5 = no skill."""
    ok = domain & np.isfinite(score)
    s, p = score[ok], presence[ok]
    npos, nneg = int(p.sum()), int((~p).sum())
    if npos == 0 or nneg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ss = s[order]
    uniq, inv, cnt = np.unique(ss, return_inverse=True, return_counts=True)
    cum = np.cumsum(cnt)
    avg = (cum - cnt + 1 + cum) / 2.0
    ranks = np.empty(s.size)
    ranks[order] = avg[inv]
    return float((ranks[p].sum() - npos * (npos + 1) / 2.0) / (npos * nneg))


def validate_maps(maps, mask, presence, n_total, breaks, cell_km2):
    """maps: {name: array}. presence: bool array in the same window (flood footprint or flooded points).
    auc_roc: flooded vs not-flooded cells of the fan (needs a footprint). auc_success: the paper's cumulative-area curve
    (max 1 - flooded share/2, so compare maps with each other, not with 1)."""
    out = []
    pres = presence & mask
    for name, arr in maps.items():
        ok = mask & np.isfinite(arr)
        idx = np.digitize(arr[ok], list(breaks))
        pk = presence[ok]
        classes = [dict(area_pct=100.0 * float((idx == k).mean()), flooded_pct=100.0 * float(pk[idx == k].sum() / max(1, pk.sum())))
                   for k in range(len(breaks) + 1)]
        pct = None
        if pres.any():                                   # where the flooded cells sit in the ranking (0 = top, 100 = bottom)
            allv = np.sort(arr[ok])
            pv = arr[pres & np.isfinite(arr)]
            pct = float(np.mean(100.0 * (1.0 - np.searchsorted(allv, pv, side="right") / allv.size)))
        out.append(dict(map=name, flooded_cells=int(pres.sum()), covered_pct=100.0 * float(pres.sum()) / max(1, n_total),
                        flooded_share_pct=100.0 * float(pres.sum()) / max(1, int(mask.sum())),
                        auc_roc=roc_auc(arr, presence, mask), auc_success=presence_auc(arr, presence, mask),
                        mean_rank_pct=pct, classes=classes))
    return out


def exposure_window(crit, win, mask):
    r0, r1, c0, c1 = win
    b = crit.layers.get("built")
    if b is None:
        return np.where(mask, 1.0, np.nan)
    bn = norm01(np.nan_to_num(b[r0:r1, c0:c1], nan=0.0), crit.rng["built"])
    return np.where(mask, 0.5 + 0.5 * bn, np.nan)


def class_areas(R, mask, breaks, cell_km2):
    idx = np.digitize(R[mask], list(breaks))
    return [float((idx == k).sum() * cell_km2) for k in range(len(breaks) + 1)]


def run_ahp(res, P, layers=None, target_xy=None, extra=None, validation=None, log=log):
    """AHP flash-flood risk on the fan complexes found by run_pipeline.
    layers (all optional, arrays on the DEM grid): rain (event mm), hist (historical extreme rainfall), bare (0-1),
    built (0-1 built density), lakes ([dict(area_m2, volume_m3, cells)]), rain_design (True when rain is one design number).
    validation (optional): dict(presence=bool array on the DEM grid, name=str) - flood footprint or flooded points.
    Returns a dict of results."""
    layers = {k: v for k, v in (layers or {}).items() if v is not None and (k != "lakes" or len(v))}
    ctx, grid = res["ctx"], res["grid"]
    H, W = ctx.dem.shape
    step, ck = ctx.step, ctx.cell_km2
    if not res["complexes"]:
        log("AHP: no fan complex was found - nothing to analyse")
        return None
    crit = build_crit(res, P, layers, log)
    cn = {"slope", "local_relief", "upstream_area", "flow_accumulation", "ruggedness"}
    if "hist" in layers:
        cn.add("hist_rainfall")
    if "bare" in layers:
        cn.add("bare_sparse")
    wc, ic = ahp_weights(CATCH_NAMES, CATCH_MATRIX, cn)
    fn = {"hand_inv", "fan_slope", "flow_accumulation", "upstream_area", "local_relief"}
    if "rain" in layers and P["RAIN_IN_HAZARD"]:
        fn.add("event_rainfall")
    wf, icf = ahp_weights(FAN_NAMES, FAN_MATRIX, fn)
    log("AHP weights (consistency ratio must be <= 0.10):")
    log("  catchment source (%d criteria, CR %.4f): %s" % (ic["n"], ic["cr"], ", ".join("%s %.3f" % kv for kv in wc.items())))
    log("  fan hazard       (%d criteria, CR %.4f): %s" % (icf["n"], icf["cr"], ", ".join("%s %.3f" % kv for kv in wf.items())))
    for nm, info in (("catchment", ic), ("fan", icf)):
        if info["dropped"]:
            log("  (%s AHP: %s not supplied -> removed from the matrix and the weights re-solved)" % (nm, ", ".join(info["dropped"])))

    st_m = (ctx.acc * ck >= P["MASTER_STREAM_KM2"]) & ctx.valid.ravel()          # Stage-1 style network
    so = strahler(ctx.recv, ctx.order, st_m)
    pos = np.empty(ctx.recv.size, dtype=np.int64)
    pos[ctx.order] = np.arange(ctx.order.size)

    masters, cplx_out = [], []
    rain = layers.get("rain")
    for cx in res["complexes"]:
        r0, r1, c0, c1 = cx["win"]
        full = np.zeros((H, W), dtype=bool)
        full[r0:r1, c0:c1] = cx["mask"]
        members = [res["fans"][i] for i in cx["members"]]
        ex = [e for e in (extra or []) if full[max(0, e["cell"] // W - 60):e["cell"] // W + 60, max(0, e["cell"] % W - 60):e["cell"] % W + 60].any()]
        feeders = feeder_inventory(ctx, full, members, ex, P)
        for f in feeders:
            f["metrics"] = catchment_metrics(ctx, crit, so, st_m, pos, f["ups"], f["cell"], P,
                                             rain_ups=None if rain is None else rain.ravel()[f["ups"]])
            f["C"], f["parts"] = catchment_score(crit, f["ups"], f["metrics"]["area_km2"], wc, P)
            area = f["metrics"]["area_km2"]
            f["rain_mm"], f["hydro"], f["T"] = None, None, 1.0
            if rain is not None:
                mm = float(np.nanmean(rain.ravel()[f["ups"]]))
                if layers.get("rain_design") and P["STORM_ARF_EXPONENT"] > 0:       # one design number: a cloudburst covers a small area
                    mm *= min(1.0, (area / P["STORM_CORE_KM2"]) ** -P["STORM_ARF_EXPONENT"])
                f["rain_mm"] = mm
                f["hydro"] = event_hydrology(area, f["metrics"]["tc_rain_active_min"], mm, P,
                                             cn=None if layers.get("cn") is None else layers["cn"].ravel()[f["ups"]])
                f["T"] = trigger_from_peak(f["hydro"]["peak_m3s"], P)
            lk = [l for l in layers.get("lakes", []) if np.isin(l["cells"], f["ups"]).mean() >= 0.5]
            vol = sum(l["volume_m3"] if l.get("volume_m3") else 0.104 * l["area_m2"] ** 1.42 for l in lk)
            f["lakes"], f["lake_volume_m3"] = len(lk), float(vol)
            f["G"] = float(np.clip(vol / P["GLOF_REF_VOLUME_M3"], 0.0, 1.0)) if lk else 0.0
        E = exposure_window(crit, cx["win"], cx["mask"])
        cplx_out.append(dict(cx=cx, full=full, feeders=feeders, C_max=max((f["C"] for f in feeders), default=0.0),
                             E_mean=float(np.nanmean(E)), members=members))
        for f in feeders:
            row = dict(system="C%d-%s" % (cx["cid"], f["name"]), complex_id=cx["cid"], feeder=f["name"], source=f["source"],
                       feeder_x=f["x"], feeder_y=f["y"], stream_threshold_km2=P["MASTER_STREAM_KM2"], C_source=f["C"])
            row.update(f["metrics"])
            if f["hydro"]:
                row.update(event_rain_mm=f["rain_mm"], event_runoff_mm=f["hydro"]["runoff_mm"], event_tc_min=f["metrics"]["tc_rain_active_min"], event_tp_min=f["hydro"]["tp_min"],
                           event_peak_m3s=f["hydro"]["peak_m3s"], event_volume_m3=f["hydro"]["volume_m3"], trigger_T=f["T"])
            row.update(fan_area_km2=cx["area_km2"], fan_mean_slope_deg=float(np.nanmean(ctx.slope[r0:r1, c0:c1][cx["mask"]])),
                       fan_perimeter_km=raster_perimeter_km(cx["mask"], step, smooth=True),
                       fan_length_km=max(m["length_m"] for m in members) / 1000.0, fan_span_deg=max(m["span_deg"] for m in members),
                       apex_elev_m=min(m["cand"]["z"] for m in members), n_member_fans=len(members), flags=cx["flags"])
            masters.append(row)
    for co in cplx_out:
        co["priority"] = co["C_max"] * (co["E_mean"] if "built" in layers else 1.0)
    log("FAN COMPLEXES ranked by catchment-source susceptibility (C) and feeders")
    log("  %-3s %9s %5s %7s %8s  %s" % ("id", "area km2", "C max", "prior.", "feeders", "feeder catchments (km2) : C"))
    for co in sorted(cplx_out, key=lambda t: -t["priority"])[:12]:
        log("  %-3d %9.2f %5.2f %7.2f %8d  %s" % (co["cx"]["cid"], co["cx"]["area_km2"], co["C_max"], co["priority"], len(co["feeders"]),
                                                  "; ".join("%s %.1f : %.2f" % (f["name"], f["metrics"]["area_km2"], f["C"]) for f in co["feeders"])))
    tgt = P["AHP_TARGET"]
    if target_xy is not None:
        d = []
        rr, cc = grid.cell_of(*target_xy)
        for co in cplx_out:
            r0, r1, c0, c1 = co["cx"]["win"]
            inside = r0 <= rr < r1 and c0 <= cc < c1 and bool(co["cx"]["mask"][rr - r0, cc - c0])
            rows, cols = np.nonzero(co["cx"]["mask"])
            dist = float(np.hypot(rows + r0 - rr, cols + c0 - cc).min()) * step
            d.append((0.0 if inside else dist, co))
        chosen = [min(d, key=lambda t: t[0])[1]]
    elif tgt == "all":                                   # every major complex, largest fan first
        big = [co for co in cplx_out if co["cx"]["area_km2"] >= P["AHP_MIN_FAN_KM2"]] or cplx_out
        big.sort(key=lambda t: -t["cx"]["area_km2"])
        chosen = big[:P["AHP_MAX_TARGETS"]]
        if len(big) > len(chosen):
            log("  (%d complexes >= %.1f km2; the %d largest are analysed - raise AHP_MAX_TARGETS or name them in AHP_TARGET)" % (
                len(big), P["AHP_MIN_FAN_KM2"], len(chosen)))
    elif tgt == "top":
        chosen = [max(cplx_out, key=lambda t: t["priority"])]
    else:
        ids = {int(i) for i in (tgt if isinstance(tgt, (list, tuple, set)) else [tgt])}
        chosen = [co for co in cplx_out if co["cx"]["cid"] in ids]
    log("Target fan complex(es): %s" % ", ".join("C%d" % co["cx"]["cid"] for co in chosen))

    targets = []
    hand_d8 = get_hand(crit, ctx, P["HAND_CHANNEL_KM2"])[0]
    for co in chosen:
        cx = co["cx"]
        win, mask = cx["win"], cx["mask"]
        r0, r1, c0, c1 = win
        t = dict(cid=cx["cid"], win=win, mask=mask, feeders=co["feeders"], fan_mask=mask, corridor=np.zeros(mask.shape, dtype=bool))
        try:
            rt = route_complex(res, ctx, cx, co["feeders"], P)
            R0, R1, C0, C1 = rt["win"]
            fan_rw = np.zeros((R1 - R0, C1 - C0), dtype=bool)
            fan_rw[r0 - R0:r1 - R0, c0 - C0:c1 - C0] = mask
            cor = np.zeros_like(fan_rw)
            if P["RUNOUT_MAX_M"] > 0:
                cor = runout_corridor(ctx, rt["win"], fan_rw, rt["feed"], co["feeders"], max(f["metrics"]["area_km2"] for f in co["feeders"]),
                                      sum(f["metrics"]["area_km2"] for f in co["feeders"]), hand_d8[R0:R1, C0:C1], P)
            ext = fan_rw | cor                                           # fan + run-out corridor: where the hazard is evaluated
            rr_, cc_ = np.nonzero(ext)
            nr0, nr1, nc0, nc1 = rr_.min(), rr_.max() + 1, cc_.min(), cc_.max() + 1
            win = (R0 + nr0, R0 + nr1, C0 + nc0, C0 + nc1)
            r0, r1, c0, c1 = win
            mask = ext[nr0:nr1, nc0:nc1]
            t.update(win=win, mask=mask, fan_mask=fan_rw[nr0:nr1, nc0:nc1], corridor=cor[nr0:nr1, nc0:nc1])
            t["flow_km2"] = rt["total"][nr0:nr1, nc0:nc1]
            for f in co["feeders"]:
                if f.get("flow_km2") is not None:
                    f["flow_km2"] = f["flow_km2"][nr0:nr1, nc0:nc1]
            if cor.any():
                log("  C%d: run-out corridor adds %.2f km2 beyond the %.2f km2 fan (flood path followed up to %.0f m past the toe)" % (
                    cx["cid"], cor.sum() * ck, fan_rw.sum() * ck, P["RUNOUT_MAX_M"]))
        except MemoryError:
            t["flow_km2"] = None
            for f in co["feeders"]:
                f["flow_km2"] = None
            log("  (C%d: routing skipped - not enough memory; using single-path D8 area)" % cx["cid"])
        hw = hand_window(crit, ctx, P["HAND_CHANNEL_KM2"], win, t["flow_km2"])
        for scope in ("region", "fan"):
            Hh, n = hazard_stack(crit, ctx, win, mask, wf, hw, scope, P, ua=t["flow_km2"])
            t["H_" + scope], t["n_" + scope] = Hh, n
        t["E"] = exposure_window(crit, win, mask)
        for scope in ("region", "fan"):
            t["R_" + scope] = t["H_" + scope] * t["E"]
        t["hand"] = np.where(mask, hw, np.nan)
        t["scen"], t["worst"] = [], {}
        for f in co["feeders"]:                         # one scenario per feeder: its own flood spreading over the fan
            Hi = {}
            for scope in ("region", "fan"):
                ua = f["flow_km2"] if f.get("flow_km2") is not None else t["flow_km2"]
                rg = t["n_fan"]["_rngs"] if scope == "fan" else None
                Hi[scope] = hazard_stack(crit, ctx, win, mask, wf, hand_window(crit, ctx, P["HAND_CHANNEL_KM2"], win, ua), scope, P, ua=ua, rngs=rg)[0]
            for kind, fac in (("cloudburst", f["T"]), ("glof", f["G"])):
                if kind == "glof" and not f["lakes"]:
                    continue
                for scope in ("region", "fan"):
                    t["scen"].append(dict(feeder=f["name"], kind=kind, scope=scope, factor=fac,
                                          arr=Hi[scope] * t["E"] * (0.5 + 0.5 * f["C"]) * (0.5 + 0.5 * fac)))
        for kind in {s["kind"] for s in t["scen"]}:
            for scope in ("region", "fan"):
                t["worst"][(kind, scope)] = np.fmax.reduce([s["arr"] for s in t["scen"] if s["kind"] == kind and s["scope"] == scope])
        t["classes"] = {sc: class_areas(t["R_" + sc], mask, P["RISK_CLASS_BREAKS"], ck) for sc in ("region", "fan")}
        base = t["H_fan"][mask]                                # HAND sensitivity: effect of the channel-initiation threshold
        t["sens"] = []
        for thr in P["HAND_SENSITIVITY_KM2"]:
            hws = hand_window(crit, ctx, thr, win, t["flow_km2"])
            Hs, _ = hazard_stack(crit, ctx, win, mask, wf, hws, "fan", P, ua=t["flow_km2"], rngs=None)
            hs = Hs[mask]
            top_b, top_s = base >= np.percentile(base, 80), hs >= np.percentile(hs, 80)
            hd = hws[mask]
            t["sens"].append(dict(threshold_km2=thr, corr=float(np.corrcoef(base, hs)[0, 1]),
                                  top20_overlap=float((top_b & top_s).sum() / (top_b | top_s).sum()),
                                  hand_median_m=float(np.nanmedian(hd)), share_hand_lt5m=float(np.mean(hd < 5))))
        b = int(math.ceil(P["HANDOFF_BUFFER_M"] / step))
        t["clip_win"] = (max(0, r0 - b), min(H, r1 + b), max(0, c0 - b), min(W, c1 + b))
        t["inflow"] = []
        for f in co["feeders"]:
            dx, dy = f["dir"]
            nrm = math.hypot(dx, dy) or 1.0
            px, py = -dy / nrm, dx / nrm
            h = P["INFLOW_HALF_WIDTH_M"]
            t["inflow"].append((f["name"], f["x"] - px * h, f["y"] - py * h, f["x"] + px * h, f["y"] + py * h))
        t["validation"] = []
        if validation is not None and validation.get("presence") is not None:
            maps = {"risk_abs": t["R_region"], "risk_rel": t["R_fan"]}
            for (kind, scope), arr in t["worst"].items():
                maps["worst_%s_%s" % (kind, "abs" if scope == "region" else "rel")] = arr
            t["validation"] = validate_maps(maps, mask, validation["presence"][r0:r1, c0:c1], int(validation["presence"].sum()),
                                            P["RISK_CLASS_BREAKS"], ck)
        targets.append(t)
        log("Complex C%d risk R = H x E (%s); area per class [very low .. very high], km2:" % (
            cx["cid"], "exposure from the built layer" if "built" in layers else "no exposure layer, E = 1"))
        for sc, lab in (("region", "absolute"), ("fan", "relative")):
            log("  %-8s scale: %s" % (lab, "  ".join("%.2f" % a for a in t["classes"][sc])))
        log("  feeders (catchment km2, C source, trigger T%s):" % (", peak m3/s" if rain is not None else ""))
        for f in co["feeders"]:
            log("    %-4s %7.1f km2  C %.2f  T %.2f%s" % (f["name"], f["metrics"]["area_km2"], f["C"], f["T"],
                                                         "  Qp %.0f m3/s, %.0f mm runoff, Tp %.0f min" % (
                                                             f["hydro"]["peak_m3s"], f["hydro"]["runoff_mm"], f["hydro"]["tp_min"]) if f["hydro"] else ""))
        log("  HAND channel-threshold sensitivity (relative hazard vs the %.2f km2 default):" % P["HAND_CHANNEL_KM2"])
        for s in t["sens"]:
            log("    %.2f km2: r = %.2f, top-20%% overlap %.2f, median HAND %.0f m, %.0f%% of the fan within 5 m" % (
                s["threshold_km2"], s["corr"], s["top20_overlap"], s["hand_median_m"], 100 * s["share_hand_lt5m"]))
        for v in t["validation"]:
            log("  validation %-22s ROC-AUC %s | success-rate AUC %s | flooded cells inside this fan: %d (%.0f%% of all) | they sit in the top %s%% of the map" % (
                v["map"], "%.3f" % v["auc_roc"] if v["auc_roc"] == v["auc_roc"] else "n/a", "%.3f" % v["auc_success"] if v["auc_success"] == v["auc_success"] else "n/a",
                v["flooded_cells"], v["covered_pct"], "%.0f" % v["mean_rank_pct"] if v["mean_rank_pct"] is not None else "n/a"))
    return dict(weights=dict(catchment=(wc, ic), fan=(wf, icf)), crit=crit, complexes=cplx_out, targets=targets, master=masters)



# ======================================================================================================
# 2d. HEC-RAS HANDOFF: turns one flood path from the AHP stage into the files a HEC-RAS 2D model needs
#     (flow-area polygon, boundary lines, terrain, roughness classes, hydrographs, scenario table, README)
# ======================================================================================================
SCS_T = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6,
         2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0)                     # NRCS curvilinear unit hydrograph (peak-rate factor 484)
SCS_Q = (0.0, 0.030, 0.100, 0.190, 0.310, 0.470, 0.660, 0.820, 0.930, 0.990, 1.000, 0.990, 0.930, 0.860, 0.780, 0.680, 0.560, 0.460,
         0.390, 0.330, 0.280, 0.207, 0.147, 0.107, 0.077, 0.055, 0.040, 0.029, 0.021, 0.015, 0.011, 0.005, 0.0)


def scs_hydrograph(qp, tp_min, dt_min=1.0):
    """Flow hydrograph (minutes, m3/s) from the peak flow and time to peak. Its volume equals A x Q whenever qp = 0.208 A Q / Tp,
    because the area under this dimensionless curve is 1.336 qp Tp."""
    t = np.arange(0.0, SCS_T[-1] * tp_min + dt_min, dt_min)
    q = np.interp(t / tp_min, SCS_T, SCS_Q, right=0.0) * qp
    return t, q


def triangular_hydrograph(qp, tp_min, base_min, dt_min=1.0):
    """Triangular hydrograph of an observed event: linear rise to qp at tp, linear fall to zero at base_min."""
    t = np.arange(0.0, base_min + dt_min, dt_min)
    q = np.where(t <= tp_min, qp * t / tp_min, qp * np.clip((base_min - t) / (base_min - tp_min), 0.0, 1.0))
    return t, q


def hydro_volume_m3(t_min, q):
    return float(np.sum((q[1:] + q[:-1]) / 2.0 * np.diff(t_min * 60.0)))


def polyline_frame(px, py, X, Y):
    """For every point (X, Y): chainage s along the polyline (m), signed offset d (m, positive to the LEFT of travel) and `inside`
    (False where the nearest point lies beyond the first or last vertex)."""
    n = len(px)
    seg = np.hypot(np.diff(px), np.diff(py))
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    best = np.full(X.shape, np.inf)
    s_out = np.zeros(X.shape)
    d_out = np.zeros(X.shape)
    inside = np.ones(X.shape, dtype=bool)
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


def _smooth_resample(x, y, step, win_m, ds):
    """D8 paths are staircases: moving-average them (ends pinned) and resample every ds metres."""
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    k = max(1, int(round(win_m / step / 2.0)))
    xs, ys = x.copy(), y.copy()
    for i in range(len(x)):
        a, b = max(0, i - k), min(len(x), i + k + 1)
        m = min(i - a, b - 1 - i)                      # symmetric window: pins the ends
        xs[i], ys[i] = x[i - m:i + m + 1].mean(), y[i - m:i + m + 1].mean()
    s2 = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))])
    sn = np.arange(0.0, s2[-1], ds)
    sn = np.append(sn, s2[-1]) if s2[-1] - sn[-1] > 0.25 * ds else sn
    return np.interp(sn, s2, xs), np.interp(sn, s2, ys), sn


def _ring_simple(ring):
    """True when no two non-adjacent edges of the closed ring cross."""
    P = np.asarray(ring, dtype=float)
    n = len(P) - 1
    A, B = P[:-1], P[1:]

    def ccw(a, b, c):
        return (c[..., 1] - a[..., 1]) * (b[..., 0] - a[..., 0]) > (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])
    for i in range(n):
        j = np.arange(i + 2, n)
        j = j[~((i == 0) & (j == n - 1))]
        if j.size == 0:
            continue
        c, d = A[j], B[j]
        a, b = A[i][None, :], B[i][None, :]
        if np.any((ccw(a, c, d) != ccw(b, c, d)) & (ccw(a, b, c) != ccw(a, b, d))):
            return False
    return True


def _dp(P, tol):
    """Douglas-Peucker simplification of an open chain (n, 2)."""
    n = len(P)
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        A, B = P[i], P[j]
        v = B - A
        L2 = float(v @ v)
        Q = P[i + 1:j] - A
        t = np.clip((Q @ v) / L2, 0.0, 1.0) if L2 > 0 else np.zeros(len(Q))
        dist = np.hypot(*(Q - np.outer(t, v)).T)
        k = int(np.argmax(dist))
        if dist[k] > tol:
            keep[i + 1 + k] = True
            stack += [(i, i + 1 + k), (i + 1 + k, j)]
    return P[keep]


def _mask_ring(mask, x0, y0, res, tol):
    """Outer boundary of a boolean mask (its largest loop) as a closed (n, 2) ring in world coordinates, simplified to `tol` metres.
    Cell (r, c) covers x0 + c*res .. +res and y0 - r*res .. -res (x0, y0 = top-left corner). Built from the cell edges, so it cannot fold."""
    H, W = mask.shape
    m = np.pad(mask, 1)
    r, c = np.nonzero(mask)
    edges = {}

    def add(sel, a, b):
        for rr, cc in zip(r[sel], c[sel]):
            edges.setdefault((a[0] + cc, a[1] + rr), []).append((b[0] + cc, b[1] + rr))
    add(~m[r, c + 1], (0, 0), (1, 0))                  # top
    add(~m[r + 1, c + 2], (1, 0), (1, 1))              # right
    add(~m[r + 2, c + 1], (1, 1), (0, 1))              # bottom
    add(~m[r + 1, c], (0, 1), (0, 0))                  # left
    used, loops = set(), []
    for sk in list(edges.keys()):
        for ek in edges[sk]:
            if (sk, ek) in used:
                continue
            used.add((sk, ek))
            loop, cur, prev = [sk, ek], ek, sk
            while cur != loop[0]:
                nxt = [e for e in edges.get(cur, []) if (cur, e) not in used]
                if not nxt:
                    break
                if len(nxt) > 1:                            # two cells meeting at a corner: always turn towards the cell we are on
                    ix, iy = cur[0] - prev[0], cur[1] - prev[1]
                    nxt.sort(key=lambda e: -(ix * (e[1] - cur[1]) - iy * (e[0] - cur[0])))
                used.add((cur, nxt[0]))
                loop.append(nxt[0])
                prev, cur = cur, nxt[0]
            loops.append(loop)
    best, best_a = None, 0.0
    for lp in loops:
        xy = np.array([(x0 + cx * res, y0 - cy * res) for cx, cy in lp])
        ar = 0.5 * np.sum(xy[:-1, 0] * xy[1:, 1] - xy[1:, 0] * xy[:-1, 1])
        if abs(ar) > abs(best_a):
            best, best_a = xy, ar
    ring = best[:-1]
    k = int(np.argmax(np.hypot(ring[:, 0] - ring[0, 0], ring[:, 1] - ring[0, 1])))
    ring = np.vstack([_dp(ring[:k + 1], tol)[:-1], _dp(np.vstack([ring[k:], ring[:1]]), tol)])
    ring = np.vstack([ring, ring[:1]]) if not np.allclose(ring[0], ring[-1]) else ring
    if 0.5 * np.sum(ring[:-1, 0] * ring[1:, 1] - ring[1:, 0] * ring[:-1, 1]) < 0:
        ring = ring[::-1]
    return ring


def _snap_edge(ring, A, B):
    """The ring vertices between the ones nearest A and B, going the shorter way round (an end edge of the polygon)."""
    R = ring[:-1]
    n = len(R)
    iA = int(np.argmin(np.hypot(R[:, 0] - A[0], R[:, 1] - A[1])))
    iB = int(np.argmin(np.hypot(R[:, 0] - B[0], R[:, 1] - B[1])))
    if iA == iB:
        return np.array([A, B])

    def path(i, j):
        idx, k = [i], i
        while k != j:
            k = (k + 1) % n
            idx.append(k)
        return R[idx]
    p1, p2 = path(iA, iB), path(iB, iA)
    ln = lambda p: float(np.hypot(*np.diff(p, axis=0).T).sum())
    return p1 if ln(p1) <= ln(p2) else p2


def ras_reach(ctx, grid, f, P, line_xy=None):
    """The flood path of one feeder as a HEC-RAS 2D flow area: a smoothed centreline (the D8 path from the feeder entry, or the line
    you drew along the real stream: line_xy) and the area within a funnel-shaped distance of it (narrow at the inflow, RAS_HALF_WIDTH_M
    further down). The outline is traced from that area, so it cannot fold over itself however sinuous the stream is."""
    H, W = ctx.dem.shape
    why = ""
    if line_xy is not None:
        xy = np.asarray(line_xy, dtype=float)
        d = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xy[:, 0]), np.diff(xy[:, 1])))])
        i0 = int(np.searchsorted(d, P["RAS_START_M"]))
        i1 = int(min(len(xy) - 1, np.searchsorted(d, P["RAS_START_M"] + P["RAS_REACH_KM"] * 1000.0)))
        if d[-1] < P["RAS_START_M"] + P["RAS_REACH_KM"] * 1000.0:
            why = "your line is only %.1f km long" % (d[-1] / 1000.0)
    else:
        c, cells = int(f["cell"]), []
        while c >= 0 and len(cells) < 100000:
            cells.append(c)
            c = int(ctx.recv[c])
        xy = np.array([grid.rc_to_xy(p // W, p % W) for p in cells])
        d = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xy[:, 0]), np.diff(xy[:, 1])))])
        i0 = int(np.searchsorted(d, P["RAS_START_M"]))
        i1 = int(np.searchsorted(d, P["RAS_START_M"] + P["RAS_REACH_KM"] * 1000.0))
        if i1 >= len(cells):
            i1 = len(cells) - 1
            why = "the flow path leaves the DEM after %.1f km" % (d[-1] / 1000.0)
        if P["RAS_STOP_RATIO"] > 0:
            acc = ctx.acc[np.array(cells)] * ctx.cell_km2
            stop = np.nonzero((acc > P["RAS_STOP_RATIO"] * f["metrics"]["area_km2"]) & (np.arange(len(cells)) > i0))[0]
            if stop.size and stop[0] < i1:
                i1 = int(stop[0])
                why = "the flow path joins a stream with %.0f times the feeder's catchment area after %.1f km (RAS_STOP_RATIO; 0 = follow it on)" % (
                    P["RAS_STOP_RATIO"], d[i1] / 1000.0)
    i1 = max(i1, i0 + 6)
    x, y, s = _smooth_resample(xy[i0:i1 + 1, 0], xy[i0:i1 + 1, 1], ctx.step, P["RAS_SMOOTH_M"], P["RAS_VERTEX_M"])
    tx, ty = np.gradient(x), np.gradient(y)
    nrm = np.hypot(tx, ty)
    nx, ny = -ty / nrm, tx / nrm                                      # left normal
    hw = P["RAS_INFLOW_HALF_WIDTH_M"] + (P["RAS_HALF_WIDTH_M"] - P["RAS_INFLOW_HALF_WIDTH_M"]) * np.clip(s / max(P["RAS_FLARE_M"], 1.0), 0.0, 1.0)
    res_t = float(P["RAS_CELL_SIZE_M"])                               # the outline is traced on a grid of the mesh size
    pad = hw.max() + 3 * res_t
    xmin, xmax, ymin, ymax = x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad
    while ((xmax - xmin) / res_t) * ((ymax - ymin) / res_t) > 1.5e6:
        res_t *= 1.25
    x0, y0 = np.floor(xmin / res_t) * res_t, np.ceil(ymax / res_t) * res_t
    nc, nr = int(math.ceil((xmax - x0) / res_t)), int(math.ceil((y0 - ymin) / res_t))
    X = np.broadcast_to(x0 + (np.arange(nc) + 0.5) * res_t, (nr, nc))
    Y = np.broadcast_to((y0 - (np.arange(nr) + 0.5) * res_t)[:, None], (nr, nc))
    step_ = max(1, int(round(80.0 / max(P["RAS_VERTEX_M"], 1.0))))
    sel = np.unique(np.append(np.arange(0, len(x), step_), len(x) - 1))
    sf, df, insf = polyline_frame(x[sel], y[sel], X, Y)
    mask = insf & (np.abs(df) <= np.interp(sf, s[sel], hw[sel]))
    ux, uy = tx / nrm, ty / nrm
    mask &= ((X - x[-1]) * ux[-1] + (Y - y[-1]) * uy[-1] <= 0) & ((X - x[0]) * ux[0] + (Y - y[0]) * uy[0] >= 0)    # straight ends, square to the path
    for tol_f in (0.7, 0.35, 0.0):                                   # simplify less if that ever makes the outline touch itself
        ring = _mask_ring(mask, x0, y0, res_t, tol_f * res_t)
        if _ring_simple(ring):
            break
    area = 0.5 * abs(np.sum(ring[:-1, 0] * ring[1:, 1] - ring[1:, 0] * ring[:-1, 1])) / 1.0e6
    iA, iB = np.array([x[0] - nx[0] * hw[0], y[0] - ny[0] * hw[0]]), np.array([x[0] + nx[0] * hw[0], y[0] + ny[0] * hw[0]])
    oA, oB = np.array([x[-1] + nx[-1] * hw[-1], y[-1] + ny[-1] * hw[-1]]), np.array([x[-1] - nx[-1] * hw[-1], y[-1] - ny[-1] * hw[-1]])
    z = np.array([ctx.dem[grid.cell_of(px, py)] for px, py in zip(x, y)], dtype=float)
    return dict(x=x, y=y, s=s, z=z, hw=hw, hwl=hw, hwr=hw, nx=nx, ny=ny, ring=ring, inflow=_snap_edge(ring, iA, iB), outflow=_snap_edge(ring, oA, oB),
                area_km2=float(area), length_m=float(s[-1]), shrink=1.0, simple=_ring_simple(ring), smooth_m=float(P["RAS_SMOOTH_M"]),
                from_line=line_xy is not None, why_short=why, trace_res_m=res_t, ux=ux, uy=uy)


def ras_frame_masks(reach, X, Y, P):
    """Chainage / offset of every cell of a window and the masks derived from it."""
    s, d, inside = polyline_frame(reach["x"], reach["y"], X, Y)
    hwl, hwr = np.interp(s, reach["s"], reach["hwl"]), np.interp(s, reach["s"], reach["hwr"])
    area = inside & (np.where(d >= 0, d <= hwl, -d <= hwr))
    if "ux" in reach:
        area &= ((X - reach["x"][-1]) * reach["ux"][-1] + (Y - reach["y"][-1]) * reach["uy"][-1] <= 0) & ((X - reach["x"][0]) * reach["ux"][0] + (Y - reach["y"][0]) * reach["uy"][0] >= 0)
    channel = inside & (np.abs(d) <= P["RAS_CHANNEL_WIDTH_M"] / 2.0)
    return s, d, area, channel


def flow_type(melton_ratio, P):
    lo, hi = P["RAS_MELTON_BREAKS"]
    k = 0 if melton_ratio < lo else (1 if melton_ratio < hi else 2)
    return ("clear-water flood", "debris flood (hyperconcentrated)", "debris flow")[k], P["RAS_CV_BY_TYPE"][k]


def ras_hydrographs(f, P):
    """Hydrograph parameters for the feeder: observed override first, otherwise the pipeline's SCS numbers."""
    ov = P.get("RAS_OBSERVED_HYDROGRAPH")
    if ov:
        base = ov.get("base_min") or 2.67 * ov["tp_min"]
        t, q = triangular_hydrograph(ov["peak_m3s"], ov["tp_min"], base, P["RAS_DT_MIN"])
        return dict(source="observed (RAS_OBSERVED_HYDROGRAPH)", peak=float(ov["peak_m3s"]), tp_min=float(ov["tp_min"]), base_min=float(base),
                    shape="triangular", volume_m3=hydro_volume_m3(t, q), runoff_mm=None, t=t, q=q)
    h = f.get("hydro")
    if not h:
        return None
    t, q = scs_hydrograph(h["peak_m3s"], h["tp_min"], P["RAS_DT_MIN"])
    return dict(source="pipeline (SCS) for the rainfall given", peak=float(h["peak_m3s"]), tp_min=float(h["tp_min"]), base_min=float(t[-1]),
                shape="SCS curvilinear", volume_m3=hydro_volume_m3(t, q), runoff_mm=float(h["runoff_mm"]), t=t, q=q)


def ras_scenarios(hydro, ftype, cv_default, P):
    """Run table: every flow scale x every rheology option."""
    rows, k = [], 0
    rheo = [("clear water", 0.0, "Newtonian")]
    if cv_default > 0:
        rheo += [("bulked clear water", cv_default, "Newtonian + bulking"), ("non-Newtonian (Bingham)", cv_default, "Bingham")]
    for sc in P["RAS_Q_SCALES"]:
        for name, cv, model in rheo:
            k += 1
            rows.append(dict(run="R%02d" % k, flow_scale=sc, peak_m3s=hydro["peak"] * sc, volume_m3=hydro["volume_m3"] * sc, rheology=name, model=model,
                             Cv=cv, bulking_factor=1.0 / (1.0 - cv) if cv < 1 else float("nan"), n_channel=P["RAS_N_CHANNEL"],
                             hydrograph_file="hydrograph_Qx%.2f.csv" % sc))
    return rows


def ras_landcover(channel, built, P):
    """Roughness classes for RAS Mapper's land-cover layer: 1 channel, 2 fan / bare ground, 3 built-up (only with a built layer)."""
    lc = np.full(channel.shape, 2, dtype=np.uint8)
    table = [(1, "channel", P["RAS_N_CHANNEL"]), (2, "fan / bare ground", P["RAS_N_FAN"])]
    if built is not None:
        lc[np.nan_to_num(built, nan=0.0) >= P["RAS_BUILT_MIN"]] = 3
        table.append((3, "built-up", P["RAS_N_BUILT"]))
    lc[channel] = 1
    return lc, table


def ras_burn(z, s, channel, reach, depth_m):
    """Optional: lower the channel strip to depth_m below the (never rising) bed profile of the centreline. A 30 m cell is about as wide
    as many fan channels, so a DEM that coarse cannot hold a real 3-5 m deep channel; this puts a notch of that depth back."""
    bed = np.minimum.accumulate(reach["z"])
    ref = np.interp(s, reach["s"], bed)
    out = z.copy()
    out[channel] = np.minimum(z[channel], (ref - depth_m)[channel])
    return out


def _fmt_t(x):
    return "%.1f" % x if x == x else "n/a"


def ras_build(res, ahp, P, line_xy=None, log=print):
    """Everything a HEC-RAS 2D run of one flood path needs. Returns a dict (arrays + text-ready content)."""
    ctx, grid = res["ctx"], res["grid"]
    H, W = ctx.dem.shape
    tg = [t for t in ahp["targets"] if P["RAS_COMPLEX"] is None or t["cid"] == P["RAS_COMPLEX"]]
    if not tg:
        raise RuntimeError("RAS_COMPLEX = %r is not among the analysed fan complexes (%s)" % (
            P["RAS_COMPLEX"], ", ".join("C%d" % t["cid"] for t in ahp["targets"])))
    t = tg[0]
    fs = [f for f in t["feeders"] if not P["RAS_FEEDER"] or f["name"] == P["RAS_FEEDER"]]
    if not fs:
        raise RuntimeError("feeder '%s' not found in C%d (%s)" % (P["RAS_FEEDER"], t["cid"], ", ".join(f["name"] for f in t["feeders"])))
    f = max(fs, key=lambda f: f["metrics"]["area_km2"]) if not P["RAS_FEEDER"] else fs[0]
    reach = ras_reach(ctx, grid, f, P, line_xy)
    notes = []
    if not reach["simple"]:
        notes.append("WARNING: the polygon outline is not simple - fix it in QGIS (Vector > Geometry Tools > Check Validity) before using it")
    if reach["why_short"]:
        notes.append("the reach is %.1f km, not the %.1f km asked: %s" % (reach["length_m"] / 1000.0, P["RAS_REACH_KM"], reach["why_short"]))
    if reach["length_m"] < 500.0:
        notes.append("VERY SHORT reach (%.0f m): this feeder's flood joins a larger stream almost at once - model that stream instead (pick a bigger feeder, or set RAS_STOP_RATIO=0)" % reach["length_m"])
    ring = reach["ring"]
    b = 4
    ra, ca = zip(*[grid.cell_of(x, y) for x, y in ring])
    R0, R1, C0, C1 = max(0, min(ra) - b), min(H, max(ra) + b + 1), max(0, min(ca) - b), min(W, max(ca) + b + 1)
    if min(ra) < 0 or max(ra) >= H or min(ca) < 0 or max(ca) >= W:
        notes.append("the polygon reaches beyond the DEM: terrain is clipped at the DEM edge")
    cols, rows = np.arange(C0, C1), np.arange(R0, R1)
    X = np.broadcast_to(grid.x0 + (cols + 0.5) * grid.dx, (len(rows), len(cols)))
    Y = np.broadcast_to((grid.y0 - (rows + 0.5) * grid.dy)[:, None], (len(rows), len(cols)))
    s, d, area, channel = ras_frame_masks(reach, X, Y, P)
    z = ctx.dem[R0:R1, C0:C1].astype(np.float32)
    built = ahp["crit"].layers.get("built")
    lc, lc_table = ras_landcover(channel, None if built is None else built[R0:R1, C0:C1], P)
    z_burn = ras_burn(np.where(np.isfinite(z), z, np.nan), s, channel, reach, P["RAS_BURN_DEPTH_M"]).astype(np.float32) if P["RAS_BURN_DEPTH_M"] > 0 else None
    if ctx.step > P["RAS_CELL_SIZE_M"]:
        notes.append("the DEM cell (%.1f m) is coarser than the suggested mesh (%.0f m): finer mesh cells add no terrain detail - use a finer DEM if you have one" % (
            ctx.step, P["RAS_CELL_SIZE_M"]))
    mt = f["metrics"]
    ratio = (mt.get("melton_m_per_sqrt_km2", mt.get("melton", float("nan")))) / 1000.0
    ftype, cv = flow_type(ratio, P)
    hy = ras_hydrographs(f, P)
    if hy is None:
        notes.append("NO HYDROGRAPH: this feeder has no event rainfall. Give [storm] rainfall_mm (or rainfall_raster) in site.toml, or an observed hydrograph in a validation case.")
    zc = np.minimum.accumulate(reach["z"])
    k = int(np.searchsorted(reach["s"], max(0.0, reach["length_m"] - 500.0)))
    slope_out = float(max(0.001, (zc[k] - zc[-1]) / max(reach["s"][-1] - reach["s"][k], 1.0)))
    scen = ras_scenarios(hy, ftype, cv if cv > 0 else P["RAS_CV_BY_TYPE"][1], P) if hy else []
    mesh = {sz: reach["area_km2"] * 1e6 / (sz * sz) for sz in (10.0, 5.0, 3.0)}
    return dict(cid=t["cid"], feeder=f, reach=reach, win=(R0, R1, C0, C1), terrain=z, terrain_burned=z_burn, landcover=lc, lc_table=lc_table,
                area_mask=area, hydro=hy, flow_type=ftype, cv=cv, melton=ratio, scenarios=scen, slope_out=slope_out, mesh=mesh, notes=notes,
                step=ctx.step, gt_window=(grid.x0 + C0 * grid.dx, grid.dx, 0.0, grid.y0 - R0 * grid.dy, 0.0, -grid.dy))


def _csv(rows, cols):
    out = [",".join(cols)]
    for r in rows:
        out.append(",".join(("%.6g" % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols))
    return "\n".join(out) + "\n"


def _geojson(geom_type, coords, props, epsg):
    import json
    crs = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::%s" % epsg}} if epsg else None
    d = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": props, "geometry": {"type": geom_type, "coordinates": coords}}]}
    if crs:
        d["crs"] = crs
    return json.dumps(d)


def ras_texts(pkg, P, epsg=None):
    """File name -> text for everything that is not a raster or a shapefile."""
    import json
    r, hy = pkg["reach"], pkg["hydro"]
    out = {}
    ring = [[float(x), float(y)] for x, y in r["ring"]]
    out["flow_area_polygon.geojson"] = _geojson("Polygon", [ring], {"name": "2D flow area", "area_km2": round(r["area_km2"], 3)}, epsg)
    out["centerline.geojson"] = _geojson("LineString", [[float(a), float(b)] for a, b in zip(r["x"], r["y"])], {"name": "centreline / breakline"}, epsg)
    out["bc_inflow.geojson"] = _geojson("LineString", [[float(a), float(b)] for a, b in r["inflow"]], {"name": "inflow"}, epsg)
    out["bc_outflow.geojson"] = _geojson("LineString", [[float(a), float(b)] for a, b in r["outflow"]], {"name": "outflow"}, epsg)
    if hy:
        for sc in P["RAS_Q_SCALES"]:
            out["hydrograph_Qx%.2f.csv" % sc] = _csv([dict(time_min=float(a), time_hr=float(a) / 60.0, Q_m3s=float(b) * sc) for a, b in zip(hy["t"], hy["q"])],
                                                     ["time_min", "time_hr", "Q_m3s"])
        out["scenario_table.csv"] = _csv(pkg["scenarios"], ["run", "flow_scale", "peak_m3s", "volume_m3", "rheology", "model", "Cv", "bulking_factor",
                                                             "n_channel", "hydrograph_file"])
    out["landcover_manning_n.csv"] = _csv([dict(code=c, name=n, manning_n=v) for c, n, v in pkg["lc_table"]], ["code", "name", "manning_n"])
    summ = dict(complex=pkg["cid"], feeder=pkg["feeder"]["name"], epsg=epsg, centerline=dict(x=[float(v) for v in r["x"]], y=[float(v) for v in r["y"]],
                s=[float(v) for v in r["s"]], z=[float(v) for v in r["z"]]), flow_area_km2=r["area_km2"], length_m=r["length_m"],
                outflow_friction_slope=pkg["slope_out"], flow_type=pkg["flow_type"], melton_ratio=pkg["melton"], Cv_default=pkg["cv"],
                hydrograph=None if not hy else dict(source=hy["source"], peak_m3s=hy["peak"], tp_min=hy["tp_min"], base_min=hy["base_min"],
                volume_m3=hy["volume_m3"], shape=hy["shape"]), feeder_entry=[float(r["x"][0]), float(r["y"][0])], terrain_window=list(pkg["win"]),
                terrain_cell_m=pkg["step"], mesh_cells_estimate=pkg["mesh"], notes=pkg["notes"], settings={k: v for k, v in P.items() if k.startswith("RAS_")})
    out["handoff_summary.json"] = json.dumps(summ, indent=2, default=str)
    out["README.txt"] = ras_readme(pkg, P, epsg)
    return out


def ras_readme(pkg, P, epsg):
    r, hy = pkg["reach"], pkg["hydro"]
    L = []
    A = L.append
    A("HEC-RAS 2D HANDOFF - flood path of feeder %s, fan complex C%d" % (pkg["feeder"]["name"], pkg["cid"]))
    A("=" * 78)
    A("Reach %.2f km long, 2D flow area %.2f km2, bed %.0f m -> %.0f m. CRS: EPSG:%s (same as your DEM)." % (r["length_m"] / 1000, r["area_km2"], r["z"][0], r["z"][-1], epsg or "see terrain.tif"))
    A("")
    A("WHAT HEC-RAS IS DOING HERE: it solves water moving over the terrain, cell by cell and second by second, from the inflow hydrograph")
    A("at the top of the reach, and reports depth, speed and arrival time everywhere. The AHP maps only RANK places; this gives numbers.")
    A("")
    A("0. CHECK FIRST (5 minutes in QGIS)")
    A("   Load centerline.geojson and flow_area_polygon.geojson over satellite imagery. The centreline is the DEM's flow path, which can sit")
    A("   a few hundred metres from the real stream. If it does, draw a line along the real stream (from the feeder entry downstream),")
    A("   save it as a file, set [hecras] centerline to it in site.toml and re-run.")
    if pkg["notes"]:
        A("   Notes from the run:")
        for n in pkg["notes"]:
            A("     - " + n)
    A("")
    A("FILES")
    A("   terrain.tif%s          elevation (the DEM clipped to the polygon + a margin)" % (" / terrain_burned.tif" if pkg["terrain_burned"] is not None else ""))
    A("   flow_area_polygon.*       the 2D flow area; centerline.* the breakline (refine the mesh along it); bc_inflow.*, bc_outflow.* boundary lines")
    A("   landcover.tif + landcover_manning_n.csv   roughness classes (code -> Manning n)")
    A("   hydrograph_Qx*.csv        inflow hydrographs (time_min, time_hr, Q_m3s) at 1 flow scale each")
    A("   scenario_table.csv        the runs to make (flow scale x rheology)")
    A("   handoff_summary.json      everything above in one file ('fanrisk score' reads it)")
    A("")
    A("BUILD THE MODEL IN HEC-RAS (RAS Mapper). Menu wording differs between versions: HEC's '2D Modeling User's Manual' covers each step.")
    A("   1. New project. Set the projection to the terrain's CRS (EPSG:%s)." % (epsg or "of terrain.tif"))
    A("   2. RAS Mapper -> Terrains -> create a terrain from terrain.tif.")
    A("   3. Geometry -> 2D Flow Area: draw it, or import flow_area_polygon.shp. Mesh cell size about %.0f m; add centerline as a breakline with" % P["RAS_CELL_SIZE_M"])
    A("      cells about half that size along it. Estimated cells: %.0f at 10 m, %.0f at 5 m, %.0f at 3 m." % (pkg["mesh"][10.0], pkg["mesh"][5.0], pkg["mesh"][3.0]))
    A("   4. Boundary conditions: external BC lines from bc_inflow and bc_outflow (they sit on the polygon edge).")
    A("   5. Land cover: create a land-cover layer from landcover.tif, fill Manning's n from landcover_manning_n.csv.")
    A("   6. Roads, bridges, culverts: floods are often made worse where roads cross the channel. Add road embankments as")
    A("      terrain modifications or breaklines, and culverts / bridges as 2D connections, if you have them.")
    A("   7. Unsteady flow: Flow hydrograph at the inflow line (paste Q_m3s from hydrograph_Qx1.00.csv, interval 1 minute).")
    A("      Outflow line: normal depth, friction slope %.4f." % pkg["slope_out"])
    A("   8. Plan: use the full shallow-water (momentum) equation set - the flow is fast and supercritical, diffusion wave is not suitable.")
    A("      Adaptive time step with Courant about 1, dry start, run for at least %.0f minutes (peak arrives at %.0f min)." % (
        (hy["t"][-1] + 30) if hy else 120, hy["tp_min"] if hy else 0))
    A("   9. Export the MAXIMUM maps - Depth (Max), Velocity (Max) and, if offered, Arrival Time - as GeoTIFF into one folder (e.g. ras_export).")
    A("      Set the map profile to Max first. A map for one time step (its name carries a time like 00:17:00) shows only the flood at that moment,")
    A("      not its extent; 'fanrisk score' refuses such maps. The plan must run long enough for the WHOLE flood to pass, or the Max maps are cut short.")
    A("")
    A("DEBRIS FLOW (optional second pass)")
    A("   Flow type from the feeder's Melton ratio %.2f: %s." % (pkg["melton"], pkg["flow_type"]))
    A("   The flow type of a past or future event is rarely certain (flash flood, mud flow, debris flow), so scenario_table.csv runs clear water, bulked")
    A("   clear water and a non-Newtonian (Bingham) option. In HEC-RAS: Unsteady flow -> non-Newtonian / mud and debris options, enter the")
    A("   volumetric concentration Cv from the table and a bulking method that adds solids to this clear-water hydrograph. Yield stress and")
    A("   viscosity are rarely known for a site: treat them as calibration parameters. Cv values in the table are starting values only.")
    A("")
    if hy:
        A("INFLOW: %s. Peak %.0f m3/s at %.0f min, volume %.2f million m3, shape %s." % (hy["source"], hy["peak"], hy["tp_min"], hy["volume_m3"] / 1e6, hy["shape"]))
        A("   The SCS hydrograph assumes the whole catchment responds together: for a large catchment hit by a localised storm it can be well off.")
        A("   For a hindcast of a documented flood give the observed hydrograph in a validation case; for design storms / forecasts use the pipeline's.")
    A("")
    A("THEN: 'fanrisk hecras run' makes the runs and exports the maps; 'fanrisk score' compares them with observations and with the AHP map.")
    return "\n".join(L) + "\n"
