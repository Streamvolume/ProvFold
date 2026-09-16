#!/usr/bin/env python3
"""Replay the first and final replicate of every independent benchmark scenario and compare streams."""

from __future__ import annotations

from datetime import datetime, timezone

from collections import defaultdict
import csv
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import os
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from provfold.compare import compare_factorial_methods  # noqa: E402
from provfold.config import AggregationConfig, SimulationConfig  # noqa: E402
from provfold.simulate import simulate_records  # noqa: E402


GRID_PATH = Path(__file__).resolve().parent / "benchmark_grid.json"
OUTPUT = Path(os.environ.get("PROVFOLD_BENCHMARK_OUTPUT", Path(__file__).resolve().parent / "generated_run")).resolve()
GENERATED = OUTPUT / "generated"
REPORT_JSON = OUTPUT / "replay_validation.json"
REPORT_MD = OUTPUT / "replay_validation.md"


def update_digest(digest: hashlib._Hash, row: dict[str, object]) -> None:
    digest.update((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8"))


def stored_hashes(path: Path, selected: set[tuple[str, int]]) -> dict[tuple[str, int], str]:
    digests: dict[tuple[str, int], hashlib._Hash] = defaultdict(hashlib.sha256)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            key = (str(row["scenario_id"]), int(row["replicate_index"]))
            if key in selected:
                update_digest(digests[key], row)
    return {key: digest.hexdigest() for key, digest in digests.items()}


def main() -> None:
    grid = json.loads(GRID_PATH.read_text(encoding="utf-8"))
    aggregation_config = AggregationConfig.from_dict(grid["aggregation_config"])
    selected: set[tuple[str, int]] = set()
    scenario_by_id = {}
    for scenario in grid["scenarios"]:
        scenario_id = str(scenario["scenario_id"])
        scenario_by_id[scenario_id] = scenario
        selected.add((scenario_id, 1))
        selected.add((scenario_id, int(scenario["replicates"])))

    stored = {
        "truth": stored_hashes(GENERATED / "truth.jsonl.gz", selected),
        "cohort_effects": stored_hashes(GENERATED / "cohort_effects.jsonl.gz", selected),
        "records": stored_hashes(GENERATED / "records.jsonl.gz", selected),
        "predictions": stored_hashes(GENERATED / "predictions.jsonl.gz", selected),
    }
    audit = {}
    with (OUTPUT / "replicate_audit.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["scenario_id"], int(row["replicate_index"]))
            if key in selected:
                audit[key] = row["simulation_hash"]

    mismatches: list[dict[str, object]] = []
    for scenario_index, scenario in enumerate(grid["scenarios"], start=1):
        scenario_id = str(scenario["scenario_id"])
        parameters = dict(grid["base_simulation"])
        parameters.update(scenario["overrides"])
        for replicate_index in (1, int(scenario["replicates"])):
            key = (scenario_id, replicate_index)
            seed = int(grid["seed_schedule"]["base_seed"]) + scenario_index * int(grid["seed_schedule"]["scenario_stride"]) + replicate_index
            simulation_parameters = dict(parameters)
            simulation_parameters.update({"simulation_id": f"{scenario_id}_r{replicate_index:03d}", "seed": seed})
            simulation_config = SimulationConfig.from_dict(simulation_parameters)
            simulation = simulate_records(simulation_config)
            scenario_config = aggregation_config
            if simulation_config.primary_comparisons_per_family > 1:
                scenario_config = replace(
                    aggregation_config,
                    primary_comparisons=tuple(
                        f"case_vs_control_panel_{index:02d}"
                        for index in range(1, simulation_config.primary_comparisons_per_family + 1)
                    ),
                )
            predictions = compare_factorial_methods(simulation.records, scenario_config)
            prefix = {"scenario_id": scenario_id, "scenario_role": scenario["role"], "replicate_index": replicate_index, "seed": seed}
            collections = {
                "truth": [{**prefix, **row} for row in simulation.truth],
                "cohort_effects": [{**prefix, **row} for row in simulation.cohort_effects],
                "records": [{**prefix, **row} for row in simulation.records],
                "predictions": [
                    {
                        **prefix,
                        "method": row["method"],
                        "taxon_key": row["taxon_key"],
                        "selected_direction": row["selected_direction"],
                        "signal_classification": row["signal_classification"],
                        "voting_unit": row["voting_unit"],
                        "taxon_key_policy": row["taxon_key_policy"],
                        "comparison_policy": row["comparison_policy"],
                        "scope_policy": row["scope_policy"],
                    }
                    for row in predictions
                ],
            }
            for stream, rows in collections.items():
                digest = hashlib.sha256()
                for row in rows:
                    update_digest(digest, row)
                observed = stored[stream].get(key, "")
                if digest.hexdigest() != observed:
                    mismatches.append({"scenario_id": scenario_id, "replicate_index": replicate_index, "stream": stream, "expected": observed, "replayed": digest.hexdigest()})
            if simulation.simulation_hash != audit.get(key):
                mismatches.append({"scenario_id": scenario_id, "replicate_index": replicate_index, "stream": "simulation_hash", "expected": audit.get(key, ""), "replayed": simulation.simulation_hash})

    report = {
        "schema_version": "1.0",
        "executed_at": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if not mismatches else "FAIL",
        "scenario_count": len(grid["scenarios"]),
        "replayed_replicate_count": len(selected),
        "stream_count": 4,
        "simulation_hashes_checked": len(selected),
        "mismatch_count": len(mismatches),
        "mismatch_examples": mismatches[:20],
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    REPORT_MD.write_text(
        "# independent benchmark deterministic replay validation\n\n"
        f"Status: `{report['status']}`\n\n"
        f"Replayed the first and final replicate of all {report['scenario_count']} scenarios ({report['replayed_replicate_count']} replicates). Four stored streams and the simulation hash were compared; mismatches: {report['mismatch_count']}.\n",
        encoding="utf-8",
    )
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
