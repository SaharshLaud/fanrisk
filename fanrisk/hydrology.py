# -*- coding: utf-8 -*-
"""Curve numbers from land cover and soil group (so runoff does not have to be guessed for a new place).

Inputs are two free global rasters:
    land cover   ESA WorldCover (10 m; class codes 10 ... 100)
    soil group   HYSOGs250m hydrologic soil groups (Ross et al. 2018; codes 1-4 = A-D, 11-14 = dual groups A/D ... D/D)
The table gives the SCS curve number (antecedent moisture condition II) for each pair. Values follow USDA TR-55 (1986), Table 2-2, for the
closest cover type. They are STARTING VALUES: replace them with [runoff.cn_table] in site.toml where local tables exist.
"""
import numpy as np

# WorldCover class -> (name, CN for soil group A, B, C, D)
CN_TABLE = {
    10: ("tree cover", (36, 60, 73, 79)),                   # woods, fair condition
    20: ("shrubland", (35, 56, 70, 77)),                    # brush, fair
    30: ("grassland", (49, 69, 79, 84)),                    # pasture / range, fair
    40: ("cropland", (67, 78, 85, 89)),                     # row crops, straight rows, good
    50: ("built-up", (77, 85, 90, 92)),                     # residential, about 65 % impervious
    60: ("bare / sparse vegetation", (77, 86, 91, 94)),     # fallow, bare soil
    70: ("snow and ice", (98, 98, 98, 98)),
    80: ("permanent water", (98, 98, 98, 98)),
    90: ("herbaceous wetland", (98, 98, 98, 98)),
    95: ("mangroves", (98, 98, 98, 98)),
    100: ("moss and lichen", (63, 77, 85, 88)),             # desert shrub, poor (sparse cover on rock)
}
HSG_INDEX = {1: 0, 2: 1, 3: 2, 4: 3, 11: 3, 12: 3, 13: 3, 14: 3}     # dual groups are taken as D (undrained): the cautious choice


def cn_from_landcover(landcover, soil_group, table=None, default_group=None):
    """Arrays of WorldCover classes and HYSOGs soil groups (same grid) -> curve-number array (NaN where either is unknown).
    default_group: 'A'..'D' used where the soil raster has no data (None = leave NaN)."""
    table = {int(k): v for k, v in (table or {}).items()} or {k: v[1] for k, v in CN_TABLE.items()}
    lc = np.where(np.isfinite(landcover), landcover, -1).astype(np.int64)
    sg = np.where(np.isfinite(soil_group), soil_group, -1).astype(np.int64)
    idx = np.full(sg.shape, -1, dtype=np.int64)
    for code, k in HSG_INDEX.items():
        idx[sg == code] = k
    if default_group:
        idx[idx < 0] = "ABCD".index(str(default_group).upper())
    cn = np.full(lc.shape, np.nan)
    for code, vals in table.items():
        sel = (lc == code) & (idx >= 0)
        cn[sel] = np.asarray(vals, dtype=float)[idx[sel]]
    return cn


def describe_cn(cn, valid):
    v = cn[valid & np.isfinite(cn)]
    if v.size == 0:
        return "no valid curve number on the DEM"
    return "curve number: mean %.0f, range %.0f-%.0f, covers %.0f %% of the DEM" % (v.mean(), v.min(), v.max(), 100.0 * v.size / max(1, int(valid.sum())))
