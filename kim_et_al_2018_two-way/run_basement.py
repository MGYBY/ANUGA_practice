"""Run the schematic basement model using transfer CSVs from run_ground.py."""
from __future__ import annotations

from pathlib import Path

from atm import BasementStairSource, TransferHydrograph
from case_setup import OUTPUT_INTERVAL, RESULTS_DIR, RUN_DURATION, STAIRS, build_basement_domain


def run_basement() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    domain = build_basement_domain()
    sources = []

    for stair in STAIRS:
        transfer_csv = RESULTS_DIR / f"transfer_{stair['label']}.csv"
        hydrograph = TransferHydrograph(transfer_csv)
        sources.append(BasementStairSource(domain, stair, hydrograph))
        print(f"Attached {stair['label']} source from {transfer_csv}; volume={hydrograph.total_volume:.3f} m3")

    print("Running basement model")
    print(domain.statistics())
    for t in domain.evolve(yieldstep=OUTPUT_INTERVAL, finaltime=RUN_DURATION):
        domain.print_timestepping_statistics()
        print("  " + " | ".join(s.status() for s in sources))

    for source in sources:
        path = RESULTS_DIR / f"source_{source.label}.csv"
        source.write_csv(path)
        print(f"Wrote {path}")


if __name__ == "__main__":
    run_basement()
