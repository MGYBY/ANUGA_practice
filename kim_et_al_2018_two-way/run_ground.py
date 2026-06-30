"""Run the schematic ground-floor model and write staircase transfer CSVs."""
from __future__ import annotations

from pathlib import Path

from atm import GroundStairSink, write_transfer_summary
from case_setup import OUTPUT_INTERVAL, RESULTS_DIR, RUN_DURATION, STAIRS, build_ground_domain


def run_ground() -> list[Path]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    domain = build_ground_domain()
    sinks = [GroundStairSink(domain, stair) for stair in STAIRS]

    print("Running ground-floor model")
    print(domain.statistics())
    for t in domain.evolve(yieldstep=OUTPUT_INTERVAL, finaltime=RUN_DURATION):
        domain.print_timestepping_statistics()
        print("  " + " | ".join(s.status() for s in sinks))

    csv_paths: list[Path] = []
    for sink in sinks:
        path = RESULTS_DIR / f"transfer_{sink.label}.csv"
        sink.write_csv(path)
        csv_paths.append(path)
        print(f"Wrote {path}")

    summary = RESULTS_DIR / "transfer_summary.csv"
    write_transfer_summary(summary, csv_paths)
    print(f"Wrote {summary}")
    return csv_paths


if __name__ == "__main__":
    run_ground()
