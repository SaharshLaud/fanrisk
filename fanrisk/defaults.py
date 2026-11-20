# -*- coding: utf-8 -*-
"""Method defaults of the fanrisk engine.

These are METHOD parameters: they describe how fans are found and how hazard is scored, and are meant to work unchanged on a new DEM.
Nothing here is calibrated to one place. Numbers that depend on climate / terrain type live in a region profile (fanrisk/profiles/*.toml),
numbers of one site in its site.toml, and observations of one flood in a validation case (fanrisk/cases/*.toml).
Precedence: defaults < region profile < site.toml sections < [advanced] in site.toml.
"""

ENGINE_DEFAULTS = dict(
    # ---- input / output -----------------------------------------------------

    # ---- DEM preparation ----------------------------------------------------
    NODATA_VALUES=[],                # extra nodata values (NaN, the declared nodata, -9999, -32768 always count)
    ZERO_IS_NODATA="auto",           # zero-filled margins next to much higher terrain: "auto" | True | False
    REPROJECT_GEOGRAPHIC=True,       # a DEM in degrees is warped to UTM automatically

    # ---- which apexes? ------------------------------------------------------
    APEX_MODE="auto",                # "auto"   = detect apexes (MANUAL_APEXES only validate the detector)
                                     # "manual" = only MANUAL_APEXES become fans
                                     # "both"   = MANUAL_APEXES + detected apexes
    MANUAL_APEXES=[],                # [(id, system, lon, lat), ...] WGS84
    SNAP_RADIUS_M=250.0,             # manual apex -> nearest sizeable stream cell within this radius
    MERGE_DIST_M=150.0,              # detected apex closer than this to a manual apex counts as a duplicate
    APEX_METHOD="openness",          # "openness" (valley mouth) | "paper" (steep-zone exit) | "both"

    # ---- hydrological analysis ----------------------------------------------
    FLOW_ACC_THRESHOLD_KM2=None,     # None = auto: max(0.5 km2, 400 cells)

    # ---- apex detector "openness": channel leaves a confined valley ---------
    OPEN_RADIUS_M=1000.0,            # terrain is "walled in" if it rises above OPEN_BLOCK_DEG within this distance
    OPEN_BLOCK_DEG=8.0,
    OPEN_SMOOTH_M=130.0,             # smoothing of the openness signal along the channel
    OPEN_PLATEAU_M=2000.0,           # apex = first point reaching OPEN_FRACTION of the openness found within this
    OPEN_FRACTION=0.90,              #   distance downstream ...
    OPEN_MIN_RISE=5.0,               #   ... after rising by at least this many of 16 directions

    # ---- apex detector "paper": steep zone + feeder-channel curvature --------
    CONV_RADIUS_M=150.0,             # [paper: convolution filter radius]
    SLOPE_THR_1_DEG=None,            # [paper: slope threshold 1]  None = Otsu threshold of the smoothed slope map
    CURV_BASELINE_CELLS=3,
    CURV_THRESHOLD=0.0,              # [paper: curvature threshold] mean curvature <= this (negative = concave)
    TRUNK_LENGTH_M=600.0,
    MIN_FEEDER_LENGTH_M=300.0,
    EDGE_BUFFER_CELLS=5,             # ignore apexes this close to the DEM edge / nodata

    # ---- fan surface and radial profiles --------------------------------------
    CHANNEL_FILL_M=80.0,             # channels/gullies narrower than ~2x this are filled before profiling
    RADIAL_STEP_DEG=5.0,
    MAX_FAN_LENGTH_M=8000.0,         # raise it if the fan layer reports truncated profiles (n_trunc)
    REF_LEN_M=1200.0,                # fan reference slope = steep end of the mean slope over this first stretch
    TOE_FRACTION=0.6,                # [paper: slope threshold 2, variable] toe = slope stays below this fraction
    TOE_MIN_DEG=1.0,                 #   of the fan's reference slope (clipped to this range) ...
    TOE_MAX_DEG=3.0,
    MIN_RUN_M=300.0,                 #   ... for at least this length
    SLOPE_WINDOW_M=300.0,            # window of the along-profile slope
    PROFILE_SMOOTH_M=90.0,
    GRACE_M=300.0,                   # the toe test starts after this distance from the apex
    RISE_TOLERANCE_M=10.0,           # profile climbing this much above its lowest point ends the profile
    MIN_RADIAL_RELIEF_M=3.0,
    MIN_RADIAL_LENGTH_M=90.0,
    RADIAL_MEDIAN_WINDOW=3,          # median filter over neighbouring radials (odd; 1 = off)

    # ---- conical surface x DEM -----------------------------------------------
    CONE_TOLERANCE_M=3.0,            # surface may exceed the cone by tol_m + tol_frac x local relief
    CONE_TOLERANCE_FRAC=0.12,
    MASK_SMOOTH_CELLS=1,
    HEAD_RECOVERY_M=0.0,             # optional: fan may extend this far UP above the apex (0 = off; APEX_REFINE does it better)

    # ---- apex refinement along the channel (detected apexes only) --------------
    APEX_REFINE=True,
    APEX_REFINE_UP_M=2000.0,         # the detector's apex may be this far below the fan head
    APEX_REFINE_STEP_M=150.0,
    APEX_REFINE_GAIN=0.05,           # keep moving up while a step gains at least this fraction of fan area

    # ---- acceptance filters ---------------------------------------------------
    MIN_FAN_AREA_KM2=0.1,
    MIN_FAN_LENGTH_M=300.0,
    MIN_FAN_SLOPE_DEG=1.0,           # a fan must be at least this steep near its apex (rejects plains and bowls)
    MIN_SPAN_DEG=30.0,               # rejects corridor-shaped valley floors
    MIN_FAN_RAYS=4,
    SUPPRESS_APEX_INSIDE_FAN=True,   # detected apex inside an already accepted fan is dropped
    MAX_FAN_OVERLAP=0.8,             # detected fan with >= this fraction inside earlier fans = duplicate (>1 = off)
    MAX_NESTED_OVERLAP=0.5,          # same, but for an apex UPSTREAM of an accepted apex (its fan is mostly the same fan)
    MAX_FAN_CATCHMENT_RATIO=1.5,     # fan larger than this x its own catchment = apex sits on another fan (rejected)

    # ---- coalescence into fan complexes --------------------------------------
    MERGE_FANS=True,
    MERGE_GAP_M=60.0,                # fans closer than this are merged and the gap is closed
    MERGE_MIN_OVERLAP=0.02,          # or overlapping by this fraction of the smaller fan

    # ---- AHP flash-flood risk mapping (runs on the fan complexes found above) -----------------------
    RUN_AHP=True,
    AHP_TARGET="all",                # "all" = every fan complex >= AHP_MIN_FAN_KM2 (largest first) | "top" = highest source susceptibility | [ids]
    AHP_TARGET_LONLAT=None,          # (lon, lat): analyse the fan complex containing / nearest to this point (overrides AHP_TARGET)
    AHP_MIN_FAN_KM2=1.0,             # AHP_TARGET = "all" skips complexes smaller than this ...
    AHP_MAX_TARGETS=6,               # ... and analyses at most this many
    EXTRA_FEEDERS=[],                # [(name, lon, lat)]: extra pour points ON a tributary
    FEEDER_MIN_KM2=1.0,              # a stream entering a fan counts as a feeder from this catchment area ...
    FEEDER_MIN_SEP_M=300.0,          # ... and distinct entries must be this far apart
    MAX_FEEDERS=12,                  # small (2-3 km2) streams can be the lethal ones: do not cap them out
    MASTER_STREAM_KM2=0.05,          # channel threshold of the master morphometric table
    HAND_CHANNEL_KM2=0.25,           # channel-initiation threshold for the local HAND criterion
    HAND_SENSITIVITY_KM2=(0.05, 0.25, 0.5, 1.0, 2.0),
    LOCAL_RELIEF_RADIUS_M=500.0,     # local relief = max - min inside +- this distance
    CATCH_AREA_MIN_KM2=0.5,          # catchment source criterion (area): fixed log scale, 0 at this area ...
    CATCH_AREA_REF_KM2=100.0,        # ... and 1 at this area (portable between DEMs, unlike min-max over a few systems)
    # fan-hazard criteria on ABSOLUTE physical scales (the same score means the same thing on every fan and every DEM)
    HAND_REF_M=10.0,                 # height above the channel that scores 0 (about twice the largest credible flow depth; region profile)
    FAN_SLOPE_REF_DEG=10.0,          # fan slope that scores 1 (upper end of alluvial-fan slopes)
    RELIEF_REF_M=300.0,              # local relief that scores 1
    FLOW_AREA_MIN_KM2=0.05,          # routed contributing area: 0 at this area ...
    FLOW_AREA_REF_KM2=100.0,         # ... and 1 at this area
    # how a flood spreads over the fan (multiple-flow-direction routing of every feeder + the fan's own hillsides)
    ROUTE_EXPONENT=2.0,              # Holmgren exponent: 1 = spreads widely, 4+ = nearly one path
    ROUTE_SURFACE_BLEND=0.0,         # surface the flood is routed on: 0 = real DEM (incised channels capture the flow, as in HEC-RAS), 1 = channel-filled fan surface (v0.1 behaviour)
    ROUTE_MARGIN_M=1000.0,           # routing window around the fan complex
    # event hydrology of each feeder (SCS triangular hydrograph)
    STORM_DURATION_MIN=15.0,         # storm duration in minutes (site input: [storm] duration_min)
    RUNOFF_COEFF=1.0,                # RUNOFF_METHOD 'coefficient': share of the rain that runs off (1 = bare rock; lower for soil / vegetation)
    STORM_ARF_EXPONENT=0.0,          # design storm only: depth x (area/STORM_CORE_KM2)**-exponent for catchments larger than the core (0 = off)
    STORM_CORE_KM2=1.0,              # area over which a design cloudburst depth applies in full
    RUNOUT_MAX_M=3000.0,             # follow the flood this far beyond the fan toe (0 = fans only): settlements below the toe are often the ones hit
    RUNOUT_MIN_FRACTION=0.01,        # ... where a cell still carries at least this share of the largest feeder's area
    RUNOUT_HALF_WIDTH_M=400.0,       # ... the channel path is widened by this much each side (region profile)
    RUNOUT_STOP_RATIO=4.0,           # ... and stop at a river whose D8 area exceeds this multiple of the feeders' total (the trunk river)
    RAIN_ACTIVE_FRAC=0.5,            # response time is taken from the farthest cell that receives >= this share of the peak rain (localised storms)
    TRIGGER_Q_LOW_M3S=5.0,           # cloudburst trigger factor: log scale of the feeder peak flow, 0 at this discharge ...
    TRIGGER_Q_HIGH_M3S=1000.0,       # ... and 1 at this one (a very large flash flood; region profile)
    RAIN_IN_HAZARD=True,             # event rainfall is also a fan-hazard criterion when rainfall is given
    RAINFALL_REF_MM=100.0,           # event rainfall that counts as factor 1.0
    RISK_CLASS_BREAKS=(0.2, 0.4, 0.55, 0.7),  # very low | low | moderate | high | very high. Equal 0.2 steps leave 'very high' empty (no fan cell scores high on every
                                            # criterion at once); 0.7 is about the top 1 % of fan cells on the test DEMs (apex mouths, main flow paths)
    GLOF_REF_VOLUME_M3=1.0e6,        # GLOF screening: lake volume that counts as factor 1.0
    RISK_SMOOTH_M=60.0,              # light smoothing of the hazard surface (0 = off)
    HANDOFF_BUFFER_M=500.0,          # HEC-RAS terrain clip margin around the fan complex
    INFLOW_HALF_WIDTH_M=150.0,       # half-length of the inflow boundary line drawn across each feeder
    # optional inputs (names of QGIS layers; "" = not used). Rainfall/bare/built are rasters, lakes are polygons.
    RAINFALL_DESIGN_MM=None,         # or one number (mm) applied everywhere
    # validation against an observed flood (the success-rate / AUC method of Sar et al. 2025)
    SCENARIO_MAPS="auto",            # per-feeder scenario maps: True | False | "auto" (written when <= 2 complexes are analysed)
    # HEC-RAS 2D handoff: files for ONE flood path (see "HEC-RAS HANDOFF" in the header)
    RUN_RAS_HANDOFF=False,           # True = write the package into <OUTPUT_DIR>/ras_C<id>_<feeder>
    RAS_COMPLEX=None,                # fan complex id (None = the first analysed complex)
    RAS_FEEDER="",                   # feeder name such as "F1" ("" = the largest feeder of that complex)
    RAS_START_M=0.0,                 # start the reach this far below the feeder entry
    RAS_REACH_KM=4.5,                # reach length (site input: [hecras] reach_km)
    RAS_HALF_WIDTH_M=700.0,          # half-width of the 2D flow area (site input: wider than the widest credible flood)
    RAS_INFLOW_HALF_WIDTH_M=200.0,   # half-width at the inflow end; the area funnels out from there over RAS_FLARE_M
    RAS_FLARE_M=800.0,
    RAS_SMOOTH_M=300.0,              # smoothing of the DEM flow path (it is a staircase)
    RAS_STOP_RATIO=6.0,              # end the reach where the flow path joins a stream with this many times the feeder's catchment area (0 = never)
    RAS_VERTEX_M=40.0,               # vertex spacing of the polygon and centreline
    RAS_CHANNEL_WIDTH_M=40.0,        # width of the 'channel' roughness class
    RAS_BURN_DEPTH_M=0.0,            # > 0 also writes terrain_burned.tif with the channel strip lowered by this much
    RAS_CELL_SIZE_M=10.0,            # suggested HEC-RAS mesh cell size
    RAS_DT_MIN=1.0,                  # hydrograph time step, minutes
    RAS_Q_SCALES=(0.65, 1.0, 1.35),  # flow scales to run (+-35 % is the paper's uncertainty on its peak flows)
    RAS_OBSERVED_HYDROGRAPH=None,    # dict(peak_m3s=, tp_min=, base_min=) replaces the pipeline hydrograph (hindcast of an observed event: validation case)
    RAS_N_CHANNEL=0.05,              # Manning n: channel, fan / bare ground, built-up (region profile)
    RAS_N_FAN=0.035,
    RAS_N_BUILT=0.10,
    RAS_BUILT_MIN=0.3,               # built density above which a cell is 'built-up'
    RAS_MELTON_BREAKS=(0.3, 0.6),    # Melton ratio: below = clear-water flood, between = debris flood, above = debris flow (after Wilford et al. 2004)
    RAS_CV_BY_TYPE=(0.0, 0.25, 0.45),  # starting sediment concentration for each type (calibrate!)

    # ---- output ---------------------------------------------------------------
    SMOOTH_OUTLINES=True,
    WRITE_DEBUG_RASTERS=True,

    # ---- runoff (added in fanrisk): how event rain becomes runoff ---------------------------------------
    RUNOFF_METHOD="coefficient",     # "coefficient" = RUNOFF_COEFF x rain | "scs_cn" = SCS curve number (CURVE_NUMBER or a CN raster)
    CURVE_NUMBER=None,               # SCS curve number (AMC II) used when no CN raster is given
    CN_AMC=2,                        # antecedent moisture: 1 dry | 2 average | 3 wet (CN converted with the Chow et al. 1988 formulas)
    CN_IA_RATIO=0.2,                 # initial abstraction Ia = ratio x S (0.2 classic; 0.05 after Woodward et al. 2003)
)
