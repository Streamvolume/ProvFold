#!/usr/bin/env python3
"""Run the frozen independent hierarchical benchmark without parameter adaptation."""

from __future__ import annotations

from datetime import datetime, timezone

from collections import defaultdict
import csv
from dataclasses import replace
import gzip
import hashlib
import io
import os
import json
from pathlib import Path
import statistics
import sys
import time
import tracemalloc


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from provfold.benchmark import benchmark_factorial_comparisons  # noqa: E402
from provfold.compare import PROVFOLD_FACTORIAL_METHOD, compare_factorial_methods  # noqa: E402
from provfold.config import AggregationConfig, SimulationConfig  # noqa: E402
from provfold.simulate import simulate_records  # noqa: E402


GRID_PATH = Path(__file__).resolve().parent / "benchmark_grid.json"
OUTPUT = Path(os.environ.get("PROVFOLD_BENCHMARK_OUTPUT", Path(__file__).resolve().parent / "generated_run")).resolve()
GENERATED = OUTPUT / "generated"
METRICS = (
    "precision",
    "false_discovery_proportion",
    "recall",
    "specificity",
    "false_positive_rate",
    "f1",
    "false_recurrent_signal_rate",
    "heterogeneity_preservation_rate",
    "abstention_rate",
    "out_of_truth_selected_pair_count",
)


class DeterministicGzipJsonl:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.binary = path.open("wb")
        self.gzip = gzip.GzipFile(filename="", mode="wb", fileobj=self.binary, mtime=0)
        self.text = io.TextIOWrapper(self.gzip, encoding="utf-8", newline="\n")
        self.count = 0

    def write(self, row: dict[str, object]) -> None:
        self.text.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        self.count += 1

    def close(self) -> None:
        self.text.flush()
        self.text.detach()
        self.gzip.close()
        self.binary.close()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        if fields:
            writer.writeheader()
            writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def summarize(rows: list[dict[str, object]], metrics: tuple[str, ...]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["scenario_id"]), str(row["scenario_role"]), str(row["method"]))].append(row)
    output: list[dict[str, object]] = []
    for (scenario_id, role, method), group in sorted(grouped.items()):
        summary: dict[str, object] = {
            "scenario_id": scenario_id,
            "scenario_role": role,
            "method": method,
            "replicate_count": len(group),
        }
        for field in (
            "method_display_name",
            "voting_unit",
            "taxon_key_policy",
            "comparison_policy",
            "scope_policy",
        ):
            if field in group[0]:
                summary[field] = group[0][field]
        for metric in metrics:
            values = [float(row[metric]) for row in group if row.get(metric) is not None and row.get(metric) != ""]
            summary[f"{metric}_defined_replicates"] = len(values)
            summary[f"{metric}_undefined_replicates"] = len(group) - len(values)
            summary[f"{metric}_mean"] = round(statistics.fmean(values), 8) if values else None
            summary[f"{metric}_standard_deviation"] = round(statistics.stdev(values), 8) if len(values) > 1 else (0.0 if values else None)
            summary[f"{metric}_median"] = round(statistics.median(values), 8) if values else None
            summary[f"{metric}_p05"] = round(quantile(values, 0.05), 8) if values else None
            summary[f"{metric}_p95"] = round(quantile(values, 0.95), 8) if values else None
        output.append(summary)
    return output


def failure_diagnostics(performance_summary: list[dict[str, object]], grid: dict[str, object]) -> dict[str, object]:
    full_method = PROVFOLD_FACTORIAL_METHOD
    lookup = {(str(row["scenario_id"]), str(row["method"])): row for row in performance_summary}
    reference = lookup[("reference", full_method)]
    mandatory = {str(item["scenario_id"]) for item in grid["scenarios"] if item["role"] == "mandatory_failure_region"}
    rows = []
    failure_present = False
    for scenario_id in sorted(mandatory):
        current = lookup[(scenario_id, full_method)]
        recall_change = float(current["recall_mean"]) - float(reference["recall_mean"])
        abstention_change = float(current["abstention_rate_mean"]) - float(reference["abstention_rate_mean"])
        qualifies = recall_change <= -0.20 or abstention_change >= 0.20
        failure_present = failure_present or qualifies
        rows.append({
            "scenario_id": scenario_id,
            "recall_change_from_reference": round(recall_change, 8),
            "abstention_change_from_reference": round(abstention_change, 8),
            "mandatory_failure_criterion_met": qualifies,
        })

    uniform = True
    for scenario in grid["scenarios"]:
        scenario_id = str(scenario["scenario_id"])
        full = lookup[(scenario_id, full_method)]
        others = [row for (name, method), row in lookup.items() if name == scenario_id and method != full_method]
        if any(
            float(full["f1_mean"]) < float(other["f1_mean"])
            or float(full["false_recurrent_signal_rate_mean"]) > float(other["false_recurrent_signal_rate_mean"])
            for other in others
        ):
            uniform = False
            break
    return {
        "mandatory_failure_region_rows": rows,
        "mandatory_failure_region_present": failure_present,
        "uniform_optimality_inspection_triggered": uniform,
        "status": "PASS" if failure_present and not uniform else "INSPECTION_REQUIRED",
    }


def main() -> None:
    grid = json.loads(GRID_PATH.read_text(encoding="utf-8"))
    if grid.get("status") not in {"FROZEN_BEFORE_PRODUCTION_RESULTS", "AMENDED_BEFORE_PUBLIC_RELEASE", "RELEASED_CONFIGURATION"}:
        raise SystemExit("grid is neither frozen, transparently amended nor a released configuration")
    config = AggregationConfig.from_dict(grid["aggregation_config"])
    OUTPUT.mkdir(parents=True, exist_ok=True)
    GENERATED.mkdir(parents=True, exist_ok=True)

    truth_stream = DeterministicGzipJsonl(GENERATED / "truth.jsonl.gz")
    effect_stream = DeterministicGzipJsonl(GENERATED / "cohort_effects.jsonl.gz")
    record_stream = DeterministicGzipJsonl(GENERATED / "records.jsonl.gz")
    prediction_stream = DeterministicGzipJsonl(GENERATED / "predictions.jsonl.gz")
    performance_rows: list[dict[str, object]] = []
    runtime_rows: list[dict[str, object]] = []
    replicate_audit: list[dict[str, object]] = []

    try:
        for scenario_index, scenario in enumerate(grid["scenarios"], start=1):
            scenario_id = str(scenario["scenario_id"])
            role = str(scenario["role"])
            parameters = dict(grid["base_simulation"])
            parameters.update(scenario["overrides"])
            replicates = int(scenario["replicates"])
            print(f"independent benchmark {scenario_index}/{len(grid['scenarios'])}: {scenario_id} ({replicates} replicates)", flush=True)
            for replicate_index in range(1, replicates + 1):
                seed = int(grid["seed_schedule"]["base_seed"]) + scenario_index * int(grid["seed_schedule"]["scenario_stride"]) + replicate_index
                simulation_parameters = dict(parameters)
                simulation_parameters.update({
                    "simulation_id": f"{scenario_id}_r{replicate_index:03d}",
                    "seed": seed,
                })
                simulation_config = SimulationConfig.from_dict(simulation_parameters)
                tracemalloc.start()
                start = time.perf_counter()
                simulation = simulate_records(simulation_config)
                scenario_config = config
                if simulation_config.primary_comparisons_per_family > 1:
                    scenario_config = replace(
                        config,
                        primary_comparisons=tuple(
                            f"case_vs_control_panel_{index:02d}"
                            for index in range(1, simulation_config.primary_comparisons_per_family + 1)
                        ),
                    )
                predictions = compare_factorial_methods(simulation.records, scenario_config)
                metrics = benchmark_factorial_comparisons(predictions, simulation.truth, scenario_config)
                elapsed = time.perf_counter() - start
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()

                prefix = {
                    "scenario_id": scenario_id,
                    "scenario_role": role,
                    "replicate_index": replicate_index,
                    "seed": seed,
                }
                for row in simulation.truth:
                    truth_stream.write({**prefix, **row})
                for row in simulation.cohort_effects:
                    effect_stream.write({**prefix, **row})
                for row in simulation.records:
                    record_stream.write({**prefix, **row})
                for row in predictions:
                    prediction_stream.write({
                        **prefix,
                        "method": row["method"],
                        "taxon_key": row["taxon_key"],
                        "selected_direction": row["selected_direction"],
                        "signal_classification": row["signal_classification"],
                        "voting_unit": row["voting_unit"],
                        "taxon_key_policy": row["taxon_key_policy"],
                        "comparison_policy": row["comparison_policy"],
                        "scope_policy": row["scope_policy"],
                    })
                for row in metrics:
                    performance_rows.append({**prefix, **row})
                runtime_rows.append({
                    **prefix,
                    "record_count": len(simulation.records),
                    "truth_taxon_count": len(simulation.truth),
                    "runtime_seconds": round(elapsed, 8),
                    "peak_memory_mib": round(peak / (1024 * 1024), 8),
                })
                replicate_audit.append({
                    **prefix,
                    "simulation_hash": simulation.simulation_hash,
                    "truth_rows": len(simulation.truth),
                    "cohort_effect_rows": len(simulation.cohort_effects),
                    "record_rows": len(simulation.records),
                    "prediction_rows": len(predictions),
                    "metric_rows": len(metrics),
                    "excluded": False,
                    "exclusion_reason": "",
                })
    finally:
        truth_stream.close()
        effect_stream.close()
        record_stream.close()
        prediction_stream.close()

    performance_rows.sort(key=lambda row: (str(row["scenario_id"]), int(row["replicate_index"]), str(row["method"])))
    runtime_rows.sort(key=lambda row: (str(row["scenario_id"]), int(row["replicate_index"])))
    replicate_audit.sort(key=lambda row: (str(row["scenario_id"]), int(row["replicate_index"])))
    performance_summary = summarize(performance_rows, METRICS)
    runtime_for_summary = [
        {**row, "method": "all_methods_combined"}
        for row in runtime_rows
    ]
    runtime_summary = summarize(runtime_for_summary, ("runtime_seconds", "peak_memory_mib", "record_count"))
    diagnostics = failure_diagnostics(performance_summary, grid)

    write_csv(OUTPUT / "replicate_metrics.csv", performance_rows)
    write_csv(OUTPUT / "replicate_audit.csv", replicate_audit)
    write_csv(OUTPUT / "runtime_memory.csv", runtime_rows)
    write_csv(OUTPUT / "scenario_performance_summary.csv", performance_summary)
    write_csv(OUTPUT / "runtime_memory_summary.csv", runtime_summary)
    write_csv(OUTPUT / "failure_region_summary.csv", diagnostics["mandatory_failure_region_rows"])

    generated_files = sorted(GENERATED.glob("*.jsonl.gz"))
    generated_manifest_rows = [
        {
            "file": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
            "redistribution_status": "deterministic_synthetic_stream_regenerated_on_demand",
        }
        for path in generated_files
    ]
    write_csv(OUTPUT / "generated_stream_manifest.csv", generated_manifest_rows)
    result = {
        "schema_version": "1.0",
        "executed_at": datetime.now(timezone.utc).isoformat(),
        "status": diagnostics["status"],
        "grid_status": grid["status"],
        "grid_sha256": sha256(GRID_PATH),
        "policy_hash": config.policy_hash,
        "scenario_count": len(grid["scenarios"]),
        "replicate_count": len(replicate_audit),
        "excluded_replicate_count": sum(bool(row["excluded"]) for row in replicate_audit),
        "metric_row_count": len(performance_rows),
        "configured_factorial_cell_count": 18,
        "distinct_estimator_count": 10,
        "generated_stream_rows": {
            "truth": truth_stream.count,
            "cohort_effects": effect_stream.count,
            "records": record_stream.count,
            "predictions": prediction_stream.count,
        },
        "diagnostics": diagnostics,
        "generated_stream_manifest": generated_manifest_rows,
    }
    (OUTPUT / "benchmark_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if diagnostics["status"] != "PASS":
        raise SystemExit(3)


if __name__ == "__main__":
    main()
