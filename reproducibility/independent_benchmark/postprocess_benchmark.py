#!/usr/bin/env python3
"""Rebuild threshold calibration and estimator-equivalence tables from generated streams."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import replace
from datetime import datetime, timezone
import csv, gzip, hashlib, itertools, json, math, os, statistics, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from provfold.config import AggregationConfig
from provfold.core import aggregate, evaluate_direction_counts
from provfold.compare import factorial_method_key
BENCHMARK = Path(os.environ.get('PROVFOLD_BENCHMARK_OUTPUT', Path(__file__).resolve().parent/'generated_run')).resolve()
GRID = Path(__file__).resolve().parent/'benchmark_grid.json'

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))

def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        if fields:
            writer.writeheader()
            writer.writerows(rows)

def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

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
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction

def summarise(values: list[float | None]) -> dict[str, object]:
    defined = [float(value) for value in values if value is not None]
    if not defined:
        return {
            "defined_replicates": 0,
            "undefined_replicates": len(values),
            "mean": None,
            "monte_carlo_se": None,
            "p05": None,
            "p95": None,
        }
    sd = statistics.stdev(defined) if len(defined) > 1 else 0.0
    return {
        "defined_replicates": len(defined),
        "undefined_replicates": len(values) - len(defined),
        "mean": round(statistics.fmean(defined), 8),
        "monte_carlo_se": round(sd / math.sqrt(len(defined)), 8),
        "p05": round(quantile(defined, 0.05), 8),
        "p95": round(quantile(defined, 0.95), 8),
    }

def score(selected: dict[str, str], truth: list[dict[str, object]], directions: tuple[str, ...]) -> dict[str, float | None]:
    truth_map = {str(row["taxon"]): str(row["latent_state"]) for row in truth}
    direction_set = set(directions)
    tp = fp = fn = 0
    for taxon, state in truth_map.items():
        predicted = selected.get(taxon, "")
        true_states = {state} if state in direction_set else set()
        predicted_states = {predicted} if predicted in direction_set else set()
        tp += len(true_states & predicted_states)
        fp += len(predicted_states - true_states)
        fn += len(true_states - predicted_states)
    fp += sum(taxon not in truth_map and direction in direction_set for taxon, direction in selected.items())
    selected_count = tp + fp
    precision = None if selected_count == 0 else tp / selected_count
    recall = tp / (tp + fn) if tp + fn else None
    f1 = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None
    return {
        "precision": precision,
        "false_discovery_proportion": None if selected_count == 0 else fp / selected_count,
        "recall": recall,
        "f1": f1,
        "selected_pair_count": float(selected_count),
    }

def truth_groups() -> dict[tuple[str, int], list[dict[str, object]]]:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    with gzip.open(BENCHMARK / "generated/truth.jsonl.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            grouped[(str(row["scenario_id"]), int(row["replicate_index"]))].append(row)
    return grouped

def benchmark_threshold_calibration() -> list[dict[str, object]]:
    grid = json.loads(GRID.read_text(encoding="utf-8"))
    base = AggregationConfig.from_dict(grid["aggregation_config"])
    truths = truth_groups()
    values: dict[tuple[str, str, int], list[float | None]] = defaultdict(list)
    rows_path = BENCHMARK / "generated/records.jsonl.gz"
    with gzip.open(rows_path, "rt", encoding="utf-8") as handle:
        iterator = (json.loads(line) for line in handle)
        key_function = lambda row: (str(row["scenario_id"]), int(row["replicate_index"]))
        for key, group in itertools.groupby(iterator, key_function):
            records = list(group)
            primary_comparisons = sorted({
                str(row["comparison_id"])
                for row in records
                if str(row.get("evidence_tier", "")) == "primary"
            })
            reference_config = replace(
                base,
                analysis_id=f"calibration_{key[0]}_{key[1]}",
                primary_comparisons=tuple(primary_comparisons),
                minimum_support=1,
            )
            result = aggregate(records, reference_config)
            state_inputs = []
            for row in result.taxon_summaries:
                state_inputs.append((
                    str(row["target_taxon"]),
                    json.loads(str(row["direction_counts_json"])),
                    int(row["conflict_family_count"]),
                ))
            max_support = min(8, max(1, len({str(row["family_id"]) for row in records if str(row.get("family_id", ""))})))
            for support in range(1, max_support + 1):
                config = replace(reference_config, minimum_support=support)
                selected: dict[str, str] = {}
                for taxon, counts, conflicts in state_inputs:
                    direction, _, _ = evaluate_direction_counts(counts, conflicts, config)
                    if direction:
                        selected[taxon] = direction
                metrics = score(selected, truths[key], config.direction_states)
                for metric, value in metrics.items():
                    values[(key[0], metric, support)].append(value)

    output: list[dict[str, object]] = []
    for (scenario, metric, support), observed in sorted(values.items()):
        summary = summarise(observed)
        output.append({
            "scenario_id": scenario,
            "estimator": "M8_provenance_group_resolved_taxon_collapsed",
            "minimum_support": support,
            "maximum_opposition": 0,
            "metric": metric,
            **summary,
        })
    return output

def equivalence_audit() -> tuple[list[dict[str, object]], dict[str, object]]:
    methods = {
        method
        for method in (row["method"] for row in read_csv(BENCHMARK / "scenario_performance_summary.csv"))
    }
    classes: dict[str, list[str]] = {}
    for taxon_key in ("source_label", "resolved_target"):
        prefix = "S" if taxon_key == "source_label" else "R"
        classes[f"M{prefix}1_row_unstratified"] = [
            factorial_method_key("row", taxon_key, "ignore"),
            factorial_method_key("row", taxon_key, "separate_votes"),
        ]
        classes[f"M{prefix}2_report_unstratified"] = [factorial_method_key("report", taxon_key, "ignore")]
        classes[f"M{prefix}3_report_comparison_votes"] = [factorial_method_key("report", taxon_key, "separate_votes")]
        classes[f"M{prefix}4_provenance_group_collapsed"] = [
            factorial_method_key("family", taxon_key, "ignore"),
            factorial_method_key("family", taxon_key, "family_state"),
            factorial_method_key("report", taxon_key, "family_state"),
            factorial_method_key("row", taxon_key, "family_state"),
        ]
        classes[f"M{prefix}5_provenance_group_comparison_votes"] = [factorial_method_key("family", taxon_key, "separate_votes")]
    listed = {method for members in classes.values() for method in members}
    if listed != methods:
        raise RuntimeError(f"method-class coverage mismatch: missing={sorted(methods - listed)} extra={sorted(listed - methods)}")

    state: dict[tuple[str, int, str], set[tuple[str, str]]] = defaultdict(set)
    with gzip.open(BENCHMARK / "generated/predictions.jsonl.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["selected_direction"]:
                state[(str(row["scenario_id"]), int(row["replicate_index"]), str(row["method"]))].add(
                    (str(row["taxon_key"]), str(row["selected_direction"]))
                )
    violations = []
    replicate_keys = sorted({(scenario, replicate) for scenario, replicate, _ in state} | {
        (str(row["scenario_id"]), int(row["replicate_index"]))
        for row in read_csv(BENCHMARK / "replicate_audit.csv")
    })
    rows = []
    for class_id, members in classes.items():
        differing = 0
        for scenario, replicate in replicate_keys:
            sets = [state[(scenario, replicate, method)] for method in members]
            if any(item != sets[0] for item in sets[1:]):
                differing += 1
                violations.append((class_id, scenario, replicate))
        rows.append({
            "distinct_estimator_id": class_id,
            "configured_cell_count": len(members),
            "configured_cells": ";".join(members),
            "replicates_checked": len(replicate_keys),
            "within_class_selected_set_differences": differing,
        })
    return rows, {
        "configured_cell_count": len(methods),
        "distinct_estimator_count": len(classes),
        "replicates_checked": len(replicate_keys),
        "violation_count": len(violations),
        "status": "PASS" if not violations else "FAIL",
    }


def main():
    started = datetime.now(timezone.utc).isoformat()
    thresholds = benchmark_threshold_calibration()
    classes, validation = equivalence_audit()
    write_csv(BENCHMARK/'support_threshold_calibration.csv', thresholds)
    write_csv(BENCHMARK/'distinct_estimator_equivalence_classes.csv', classes)
    validation.update({'started_at': started, 'completed_at': datetime.now(timezone.utc).isoformat(), 'threshold_rows': len(thresholds), 'threshold_scenarios': len({r['scenario_id'] for r in thresholds}), 'threshold_replicates': len(read_csv(BENCHMARK/'replicate_audit.csv')), 'scope': 'Deterministic post-processing of generated streams using the package aggregation rules; not an independent implementation of aggregation.'})
    write_json(BENCHMARK/'postprocessing_validation.json', validation)
    if validation['violation_count']: raise SystemExit(1)
    print(json.dumps(validation, indent=2))

if __name__ == '__main__':
    main()
