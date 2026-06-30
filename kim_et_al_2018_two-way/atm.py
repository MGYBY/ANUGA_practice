"""Minimal Adaptive Transfer Method operators for the simplified ANUGA case.

The implementation follows a conservative sequential workflow:

1. The ground-floor stair operator computes Q = 0.544 A sqrt(2 g h_s).
2. The same volume Q dt is removed from the ground-floor stair region.
3. The actually removed volume is written to a CSV interval hydrograph.
4. The basement source operator reads the CSV and adds the same interval volume
   to the basement stair landing region.

This file intentionally avoids advanced options from earlier versions.  It is
therefore easier to read, but it remains a semi-3D approximation rather than a
full reproduction of Kim et al.'s in-house solver.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np

from case_setup import G

try:
    from anuga.operators.base_operator import Operator
except Exception:  # Allows syntax checking without ANUGA installed.
    class Operator(object):  # type: ignore
        def __init__(self, *args, **kwargs):
            pass


# -----------------------------------------------------------------------------
# Geometry and numerical helper functions
# -----------------------------------------------------------------------------
def polygon_indices(domain, polygon: Iterable[Iterable[float]], label: str) -> np.ndarray:
    """Return ANUGA triangle indices whose centroids are inside a polygon."""
    import anuga

    region = anuga.Region(domain, polygon=list(polygon))
    indices = np.asarray(region.get_indices() if hasattr(region, "get_indices") else region.indices)
    if indices.size == 0:
        raise RuntimeError(f"No ANUGA triangles were found inside the polygon for {label}.")
    return indices.astype(np.int64)


def sync_height_quantity(domain, indices: np.ndarray) -> None:
    """Keep ANUGA's optional height quantity consistent after direct edits."""
    if hasattr(domain, "quantities") and "height" in domain.quantities:
        stage = domain.quantities["stage"].centroid_values[indices]
        elev = domain.quantities["elevation"].centroid_values[indices]
        domain.quantities["height"].centroid_values[indices] = np.maximum(stage - elev, 0.0)


def effective_head(stage: np.ndarray, cell_area: np.ndarray, crest_m: float) -> float:
    """Return a single h_s for a non-uniform staircase region.

    Local head is h_i = max(stage_i - crest_m, 0).  The returned h_s preserves
    the area-weighted average of h_i^(3/2), which is the correct weighting for
    the default rectangular-section closure A = B h_s, because then
    Q is proportional to h_s^(3/2).
    """
    local_h = np.maximum(np.asarray(stage, dtype=float) - float(crest_m), 0.0)
    area = np.asarray(cell_area, dtype=float)
    total_area = float(np.sum(area))
    if total_area <= 0.0 or np.all(local_h <= 0.0):
        return 0.0
    mean_h_32 = float(np.sum(area * local_h ** 1.5) / total_area)
    return mean_h_32 ** (2.0 / 3.0)


def flow_area(hs: float, stair: dict) -> float:
    """Return A in Kim et al.'s formula.

    If stair['area_m2'] is supplied, it is treated as a fixed vertical/normal
    cross-sectional flow area.  Otherwise the default rectangular closure is
    used: A = width_m * h_s.
    """
    hs = max(float(hs), 0.0)
    fixed_area = stair.get("area_m2")
    if fixed_area is not None:
        return max(float(fixed_area), 0.0)
    return max(float(stair["width_m"]), 0.0) * hs


def kim_stair_discharge(hs: float, stair: dict) -> tuple[float, float]:
    """Return (Q, A) from Q = 0.544 A sqrt(2 g h_s)."""
    hs = max(float(hs), 0.0)
    if hs <= 0.0:
        return 0.0, 0.0
    area = flow_area(hs, stair)
    if area <= 0.0:
        return 0.0, area
    q = 0.544 * area * math.sqrt(2.0 * G * hs)
    return q, area


# -----------------------------------------------------------------------------
# Ground-floor sink operator
# -----------------------------------------------------------------------------
class GroundStairSink(Operator):
    """Compute staircase discharge and remove the same volume upstairs."""

    def __init__(self, domain, stair: dict):
        super().__init__(domain, description=f"Ground stair sink {stair['label']}", label=stair["label"])
        self.stair = stair
        self.label = stair["label"]
        self.indices = polygon_indices(domain, stair["ground_polygon"], self.label)
        self.cell_area = np.asarray(domain.areas[self.indices], dtype=float)
        self.cumulative_volume = 0.0
        self.rows: List[dict] = []

    def __call__(self) -> None:
        t0 = float(self.domain.get_time())
        dt = float(self.domain.get_timestep())
        if dt <= 0.0:
            return
        t1 = t0 + dt

        stage = self.stage_c[self.indices]
        elev = self.elev_c[self.indices]
        depth = np.maximum(stage - elev, 0.0)
        available_volume = float(np.sum(depth * self.cell_area))

        hs = effective_head(stage, self.cell_area, self.stair["crest_m"])
        q_formula, area = kim_stair_discharge(hs, self.stair)
        requested_volume = q_formula * dt
        removed_volume = min(requested_volume, available_volume) if requested_volume > 0.0 else 0.0
        q_applied = removed_volume / dt if dt > 0.0 else 0.0

        if removed_volume > 0.0 and available_volume > 0.0:
            fraction = removed_volume / available_volume
            self.stage_c[self.indices] -= fraction * depth
            self.xmom_c[self.indices] *= 1.0 - fraction
            self.ymom_c[self.indices] *= 1.0 - fraction
            sync_height_quantity(self.domain, self.indices)
            if hasattr(self.domain, "fractional_step_volume_integral"):
                self.domain.fractional_step_volume_integral -= removed_volume

        self.cumulative_volume += removed_volume
        self.rows.append(
            {
                "t_start_s": t0,
                "t_end_s": t1,
                "dt_s": dt,
                "hs_m": hs,
                "flow_area_A_m2": area,
                "q_formula_m3s": q_formula,
                "q_applied_m3s": q_applied,
                "available_volume_m3": available_volume,
                "removed_volume_m3": removed_volume,
                "cumulative_removed_m3": self.cumulative_volume,
            }
        )

    def status(self) -> str:
        if not self.rows:
            return f"{self.label}: Q=0.000 m3/s"
        row = self.rows[-1]
        return (
            f"{self.label}: h_s={row['hs_m']:.3f} m, "
            f"A={row['flow_area_A_m2']:.3f} m2, "
            f"Q={row['q_applied_m3s']:.3f} m3/s"
        )

    def write_csv(self, path: Path) -> None:
        write_rows(path, self.rows)


# -----------------------------------------------------------------------------
# Transfer hydrograph and basement source operator
# -----------------------------------------------------------------------------
class TransferHydrograph:
    """Piecewise-constant hydrograph read from a ground-floor transfer CSV."""

    def __init__(self, csv_path: Path):
        self.csv_path = Path(csv_path)
        if not self.csv_path.exists():
            raise FileNotFoundError(self.csv_path)
        rows = read_rows(self.csv_path)
        if not rows:
            raise RuntimeError(f"Transfer file is empty: {self.csv_path}")
        self.t0 = np.asarray([r["t_start_s"] for r in rows], dtype=float)
        self.t1 = np.asarray([r["t_end_s"] for r in rows], dtype=float)
        self.q = np.asarray([max(0.0, r["q_applied_m3s"]) for r in rows], dtype=float)
        self.vol = np.asarray([max(0.0, r["removed_volume_m3"]) for r in rows], dtype=float)

    @property
    def total_volume(self) -> float:
        return float(np.sum(self.vol))

    @property
    def final_time(self) -> float:
        return float(np.max(self.t1))

    def volume_between(self, a: float, b: float) -> float:
        """Integrate the interval hydrograph over [a, b]."""
        left = np.maximum(self.t0, float(a))
        right = np.minimum(self.t1, float(b))
        overlap = np.maximum(right - left, 0.0)
        return float(np.sum(overlap * self.q))

    def mean_q(self, start: float, end: Optional[float] = None) -> float:
        end = self.final_time if end is None else float(end)
        if end <= start:
            return float("nan")
        return self.volume_between(start, end) / (end - start)


class BasementStairSource(Operator):
    """Add the same staircase volume to the basement landing region."""

    def __init__(self, domain, stair: dict, hydrograph: TransferHydrograph):
        super().__init__(domain, description=f"Basement stair source {stair['label']}", label=stair["label"])
        self.stair = stair
        self.label = stair["label"]
        self.hydrograph = hydrograph
        self.indices = polygon_indices(domain, stair["basement_polygon"], self.label)
        self.cell_area = np.asarray(domain.areas[self.indices], dtype=float)
        self.total_area = float(np.sum(self.cell_area))
        self.cumulative_volume = 0.0
        self.rows: List[dict] = []

    def __call__(self) -> None:
        t0 = float(self.domain.get_time())
        dt = float(self.domain.get_timestep())
        if dt <= 0.0:
            return
        t1 = t0 + dt

        added_volume = self.hydrograph.volume_between(t0, t1)
        q_applied = added_volume / dt if dt > 0.0 else 0.0
        if added_volume > 0.0 and self.total_area > 0.0:
            self.stage_c[self.indices] += added_volume / self.total_area
            # The simplified source adds vertical/down-stair inflow with zero
            # horizontal momentum.  This is volume-conservative but does not
            # reproduce Kim et al.'s full momentum-transfer term.
            sync_height_quantity(self.domain, self.indices)
            if hasattr(self.domain, "fractional_step_volume_integral"):
                self.domain.fractional_step_volume_integral += added_volume

        self.cumulative_volume += added_volume
        self.rows.append(
            {
                "t_start_s": t0,
                "t_end_s": t1,
                "dt_s": dt,
                "q_applied_m3s": q_applied,
                "added_volume_m3": added_volume,
                "cumulative_added_m3": self.cumulative_volume,
            }
        )

    def status(self) -> str:
        q = self.rows[-1]["q_applied_m3s"] if self.rows else 0.0
        return f"{self.label}: Q={q:.3f} m3/s, added={self.cumulative_volume:.3f} m3"

    def write_csv(self, path: Path) -> None:
        write_rows(path, self.rows)


# -----------------------------------------------------------------------------
# CSV helpers
# -----------------------------------------------------------------------------
def write_rows(path: Path, rows: List[dict]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: f"{v:.10f}" if isinstance(v, float) else v for k, v in row.items()})


def read_rows(path: Path) -> List[dict]:
    with Path(path).open("r", newline="") as f:
        reader = csv.DictReader(f)
        return [{k: float(v) for k, v in row.items()} for row in reader]


def write_transfer_summary(path: Path, transfer_csvs: Iterable[Path], steady_start_s: float = 50.0) -> None:
    rows = []
    for csv_path in transfer_csvs:
        h = TransferHydrograph(csv_path)
        rows.append(
            {
                "label": Path(csv_path).stem.replace("transfer_", ""),
                "total_removed_volume_m3": h.total_volume,
                "mean_q_after_50s_m3s": h.mean_q(steady_start_s),
                "max_q_m3s": float(np.max(h.q)),
            }
        )
    write_rows(path, rows)
