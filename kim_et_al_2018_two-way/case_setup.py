"""Simple case setup for the ANUGA Kim et al. (2018) ATM example.

This file contains only the parameters and schematic polygons needed to run the
case.  It is intentionally compact: edit this file when you want to change the
floor plan, the roughness, the boundary condition, or the staircase rating.

Important limitation: the original Kim et al. paper publishes floor-plan images,
but not the CAD files, mesh files, or exact stair dimensions.  The polygons below
are therefore a schematic teaching geometry, not the original KICT CAD model.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

# -----------------------------------------------------------------------------
# Basic run parameters
# -----------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

G = 9.81                    # gravitational acceleration, m/s^2
ROUGHNESS_N = 0.011         # Kim et al. smooth-concrete Manning roughness
GROUND_INLET_STAGE = 1.0    # prescribed stage/depth at the main entrance, m
RUN_DURATION = 100.0        # s
OUTPUT_INTERVAL = 1.0       # s between Python outputs and SWW frames
FLOW_ALGORITHM = "DE0"      # robust ANUGA discontinuous-elevation algorithm
MINIMUM_ALLOWED_DEPTH = 1.0e-4
MINIMUM_STORABLE_DEPTH = 1.0e-3

# Schematic floor dimensions.  Kim et al.'s Figure 5 labels the ground-floor
# length as approximately 75.3 m.  The remaining dimensions are schematic.
GROUND_LENGTH = 75.3
GROUND_WIDTH = 14.0
BASEMENT_LENGTH = 75.3
BASEMENT_WIDTH = 22.0

# Mesh controls.  Smaller values create finer triangles and longer runtimes.
GROUND_MAX_TRIANGLE_AREA = 0.35
BASEMENT_MAX_TRIANGLE_AREA = 0.45
SMALL_FEATURE_TRIANGLE_AREA = 0.05

# -----------------------------------------------------------------------------
# Staircase transfer parameters
# -----------------------------------------------------------------------------
# Kim et al. use Q = 0.544 A sqrt(2 g h_s).
# In this simplified code:
#   - h_s is computed from local stage in the stair sampling polygon.
#   - A is a vertical/normal flow cross-sectional area, not the horizontal
#     staircase polygon area.
#   - By default A = width_m * h_s, a rectangular-section closure.
#   - To prescribe a constant flow area instead, set area_m2 to a positive value.
# The widths are effective calibrated widths chosen to produce Kim-like
# staircase discharges when the local head is about 0.30 m.
STAIRS = [
    {
        "label": "STAIR_L",
        "target_q_m3s": 0.74,
        "width_m": 1.87,
        "area_m2": None,      # None means use A = width_m * h_s
        "crest_m": 0.001,      # effective hydraulic control elevation, m
    },
    {
        "label": "STAIR_R",
        "target_q_m3s": 0.81,
        "width_m": 2.05,
        "area_m2": None,
        "crest_m": 0.001,
    },
]

# -----------------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------------
Point = List[float]
Polygon = List[Point]


def rect(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    """Return a counter-clockwise rectangular polygon."""
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def wall(x0: float, y0: float, x1: float, y1: float, thickness: float = 0.25) -> Polygon:
    """Return a thin wall polygon around an axis-aligned line segment."""
    if abs(y1 - y0) < 1.0e-12:  # horizontal wall
        return rect(min(x0, x1), y0 - 0.5 * thickness, max(x0, x1), y0 + 0.5 * thickness)
    if abs(x1 - x0) < 1.0e-12:  # vertical wall
        return rect(x0 - 0.5 * thickness, min(y0, y1), x0 + 0.5 * thickness, max(y0, y1))
    raise ValueError("Only horizontal and vertical wall segments are supported in this simple geometry.")


# -----------------------------------------------------------------------------
# Ground floor schematic geometry
# -----------------------------------------------------------------------------
GROUND_OUTER = [
    [0.0, 0.0],
    [32.0, 0.0],      # main entrance starts
    [43.0, 0.0],      # main entrance ends
    [GROUND_LENGTH, 0.0],
    [GROUND_LENGTH, GROUND_WIDTH],
    [0.0, GROUND_WIDTH],
]
GROUND_BOUNDARY_TAGS = {"wall": [0, 2, 3, 4, 5], "inlet": [1]}
GROUND_INLET_POLYGON = rect(32.0, 0.05, 43.0, 1.30)
GROUND_STAIR_L_POLYGON = rect(23.0, 9.7, 29.0, 12.8)
GROUND_STAIR_R_POLYGON = rect(47.0, 7.1, 53.0, 10.2)

GROUND_WALLS = [
    wall(31.8, 0.25, 31.8, 7.8, 0.30),
    wall(43.2, 0.25, 43.2, 7.8, 0.30),
    wall(5.5, 3.9, 26.0, 3.9),
    wall(5.5, 9.1, 21.0, 9.1),
    wall(6.5, 4.0, 6.5, 8.7),
    wall(14.0, 4.0, 14.0, 8.7),
    wall(21.0, 4.0, 21.0, 8.7),
    wall(22.0, 9.3, 22.0, 13.6),
    wall(54.0, 4.0, 70.0, 4.0),
    wall(55.0, 9.2, 70.0, 9.2),
    wall(55.0, 4.1, 55.0, 8.7),
    wall(62.5, 4.1, 62.5, 8.7),
    wall(70.0, 4.1, 70.0, 8.7),
    wall(54.0, 9.4, 54.0, 13.6),
    rect(35.7, 7.2, 39.1, 8.7),
    rect(48.0, 11.7, 52.0, 12.8),
]

GROUND_INTERIOR_REGIONS = [
    [GROUND_INLET_POLYGON, SMALL_FEATURE_TRIANGLE_AREA],
    [GROUND_STAIR_L_POLYGON, SMALL_FEATURE_TRIANGLE_AREA],
    [GROUND_STAIR_R_POLYGON, SMALL_FEATURE_TRIANGLE_AREA],
]

# -----------------------------------------------------------------------------
# Basement schematic geometry
# -----------------------------------------------------------------------------
BASEMENT_OUTER = rect(0.0, 0.0, BASEMENT_LENGTH, BASEMENT_WIDTH)
BASEMENT_BOUNDARY_TAGS = {"wall": [0, 1, 2, 3]}
BASEMENT_STAIR_L_POLYGON = rect(23.0, 11.7, 28.5, 15.0)
BASEMENT_STAIR_R_POLYGON = rect(47.0, 11.4, 52.0, 14.8)

BASEMENT_WALLS = [
    wall(1.0, 8.0, 20.0, 8.0),
    wall(30.0, 8.0, 74.0, 8.0),
    wall(1.0, 15.8, 20.0, 15.8),
    wall(31.0, 15.8, 74.0, 15.8),
    wall(8.0, 0.8, 8.0, 7.6),
    wall(16.0, 0.8, 16.0, 7.6),
    wall(25.0, 0.8, 25.0, 7.6),
    wall(33.0, 0.8, 33.0, 7.6),
    wall(8.0, 16.1, 8.0, 21.2),
    wall(16.0, 16.1, 16.0, 21.2),
    wall(31.0, 16.1, 31.0, 21.2),
    rect(20.0, 15.8, 22.0, 18.8),
    rect(20.0, 8.0, 22.0, 11.0),
    wall(34.0, 8.2, 34.0, 13.2),
    wall(38.0, 13.0, 43.5, 13.0),
    wall(43.5, 13.0, 43.5, 15.5),
    wall(56.5, 8.2, 56.5, 15.5),
    wall(64.0, 8.2, 64.0, 15.5),
    wall(71.0, 8.2, 71.0, 15.5),
    wall(56.5, 16.0, 56.5, 21.2),
    wall(64.0, 16.0, 64.0, 21.2),
    wall(71.0, 16.0, 71.0, 21.2),
    rect(66.5, 16.2, 72.6, 21.0),
    rect(33.8, 16.2, 37.0, 21.0),
]

BASEMENT_INTERIOR_REGIONS = [
    [BASEMENT_STAIR_L_POLYGON, SMALL_FEATURE_TRIANGLE_AREA],
    [BASEMENT_STAIR_R_POLYGON, SMALL_FEATURE_TRIANGLE_AREA],
]

# Attach polygons to the stair parameter dictionaries.
STAIRS[0]["ground_polygon"] = GROUND_STAIR_L_POLYGON
STAIRS[0]["basement_polygon"] = BASEMENT_STAIR_L_POLYGON
STAIRS[1]["ground_polygon"] = GROUND_STAIR_R_POLYGON
STAIRS[1]["basement_polygon"] = BASEMENT_STAIR_R_POLYGON


def _apply_common_domain_settings(domain, name: str):
    """Apply common settings to an ANUGA domain."""
    domain.set_name(name)
    domain.set_datadir(str(RESULTS_DIR))
    if hasattr(domain, "set_flow_algorithm"):
        domain.set_flow_algorithm(FLOW_ALGORITHM)
    if hasattr(domain, "set_minimum_allowed_height"):
        domain.set_minimum_allowed_height(MINIMUM_ALLOWED_DEPTH)
    if hasattr(domain, "set_minimum_storable_height"):
        domain.set_minimum_storable_height(MINIMUM_STORABLE_DEPTH)
    # Discontinuous storage is clearer for finite-volume diagnostics.
    if hasattr(domain, "set_store_vertices_uniquely"):
        domain.set_store_vertices_uniquely(True)
    return domain


def build_ground_domain():
    """Build and initialise the ground-floor ANUGA domain."""
    import anuga

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    domain = anuga.create_domain_from_regions(
        GROUND_OUTER,
        boundary_tags=GROUND_BOUNDARY_TAGS,
        maximum_triangle_area=GROUND_MAX_TRIANGLE_AREA,
        interior_holes=GROUND_WALLS,
        interior_regions=GROUND_INTERIOR_REGIONS,
    )
    _apply_common_domain_settings(domain, "kim2018_ground_simple_v4")
    domain.set_quantity("elevation", 0.0, location="centroids")
    domain.set_quantity("friction", ROUGHNESS_N, location="centroids")
    domain.set_quantity("stage", expression="elevation", location="centroids")

    wall_bc = anuga.Reflective_boundary(domain)
    inlet_bc = anuga.Dirichlet_boundary([GROUND_INLET_STAGE, 0.0, 0.0])
    domain.set_boundary({"wall": wall_bc, "inlet": inlet_bc, "interior": wall_bc})
    return domain


def build_basement_domain():
    """Build and initialise the basement ANUGA domain."""
    import anuga

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    domain = anuga.create_domain_from_regions(
        BASEMENT_OUTER,
        boundary_tags=BASEMENT_BOUNDARY_TAGS,
        maximum_triangle_area=BASEMENT_MAX_TRIANGLE_AREA,
        interior_holes=BASEMENT_WALLS,
        interior_regions=BASEMENT_INTERIOR_REGIONS,
    )
    _apply_common_domain_settings(domain, "kim2018_basement_simple_v4")
    domain.set_quantity("elevation", 0.0, location="centroids")
    domain.set_quantity("friction", ROUGHNESS_N, location="centroids")
    domain.set_quantity("stage", expression="elevation", location="centroids")

    wall_bc = anuga.Reflective_boundary(domain)
    domain.set_boundary({"wall": wall_bc, "interior": wall_bc})
    return domain


def geometry_catalog() -> Dict[str, object]:
    """Return polygons used for plotting and report figures."""
    return {
        "ground_outer": GROUND_OUTER,
        "ground_walls": GROUND_WALLS,
        "ground_inlet": GROUND_INLET_POLYGON,
        "ground_stairs": [GROUND_STAIR_L_POLYGON, GROUND_STAIR_R_POLYGON],
        "basement_outer": BASEMENT_OUTER,
        "basement_walls": BASEMENT_WALLS,
        "basement_stairs": [BASEMENT_STAIR_L_POLYGON, BASEMENT_STAIR_R_POLYGON],
    }
