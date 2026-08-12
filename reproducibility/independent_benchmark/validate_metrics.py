#!/usr/bin/env python3
"""Independently recompute benchmark metrics from stored truth and predictions."""

from __future__ import annotations

from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import os


ROOT = Path(__file__).resolve().parents[2]
GRID_PATH = Path(__file__).resolve().parent / "benchmark_grid.json"
OUTPUT = Path(os.environ.get("PROVFOLD_BENCHMARK_OUTPUT", Path(__file__).resolve().parent / "generated_run")).resolve()
GENERATED = OUTPUT / "generated"
REPORT_JSON = OUTPUT / "independent_validation.json"
REPORT_MD = OUTPUT / "independent_validation.md"
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


def safe_divide(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def read_jsonl_gz(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def calculate(truth: dict[str, str], predicted: dict[str, tuple[str, str]], directions: tuple[str, str]) -> dict[str, float | int | None]:
    tp = within_truth_fp = tn = fn = 0
    heterogeneity_preserved = heterogeneity_total = abstained = 0
    direction_states = set(directions)
    for taxon, latent_state in truth.items():
        selected, classification = predicted.get(taxon, ("", "not_evaluable"))
        true_states = {latent_state} if latent_state in direction_states else set()
        predicted_states = {selected} if selected in direction_states else set()
        tp_i = len(true_states & predicted_states)
        fp_i = len(predicted_states - true_states)
        fn_i = len(true_states - predicted_states)
        tp += tp_i
        within_truth_fp += fp_i
        fn += fn_i
        tn += len(direction_states) - tp_i - fp_i - fn_i
        predicted_directional = bool(predicted_states)
        if latent_state == "heterogeneous":
            heterogeneity_total += 1
            heterogeneity_preserved += not predicted_directional
        abstained += classification in {"not_evaluable", "conflict_only"}
    out_of_truth_selected = sum(
        taxon not in truth and selected in direction_states
        for taxon, (selected, _) in predicted.items()
    )
    false_positive = within_truth_fp + out_of_truth_selected
    selected_count = tp + false_positive
    precision = None if selected_count == 0 else tp / selected_count
    false_discovery_proportion = None if selected_count == 0 else false_positive / selected_count
    recall = safe_divide(tp, tp + fn)
    return {
        "precision": None if precision is None else round(precision, 8),
        "false_discovery_proportion": None if false_discovery_proportion is None else round(false_discovery_proportion, 8),
        "recall": round(recall, 8),
        "specificity": round(safe_divide(tn, tn + within_truth_fp), 8),
        "f1": round(safe_divide(2 * tp, 2 * tp + false_positive + fn), 8),
        "false_positive_rate": round(safe_divide(within_truth_fp, within_truth_fp + tn), 8),
        "false_recurrent_signal_rate": round(safe_divide(within_truth_fp, within_truth_fp + tn), 8),
        "heterogeneity_preservation_rate": round(safe_divide(heterogeneity_preserved, heterogeneity_total), 8),
        "abstention_rate": round(safe_divide(abstained, len(truth)), 8),
        "out_of_truth_selected_pair_count": out_of_truth_selected,
        "true_positive": tp,
        "false_positive": false_positive,
        "within_truth_false_positive": within_truth_fp,
        "true_negative": tn,
        "false_negative": fn,
    }


def main() -> None:
    grid = json.loads(GRID_PATH.read_text(encoding="utf-8"))
    directions = tuple(grid["aggregation_config"]["direction_states"])
    methods = tuple(sorted(grid["comparators"]))
    expected_replicates = sum(int(item["replicates"]) for item in grid["scenarios"])
    truth: dict[tuple[str, int], dict[str, str]] = defaultdict(dict)
    for row in read_jsonl_gz(GENERATED / "truth.jsonl.gz"):
        truth[(str(row["scenario_id"]), int(row["replicate_index"]))][str(row["taxon"])] = str(row["latent_state"])
    predictions: dict[tuple[str, int, str], dict[str, tuple[str, str]]] = defaultdict(dict)
    for row in read_jsonl_gz(GENERATED / "predictions.jsonl.gz"):
        predictions[(str(row["scenario_id"]), int(row["replicate_index"]), str(row["method"]))][str(row["taxon_key"])] = (
            str(row["selected_direction"]),
            str(row["signal_classification"]),
        )

    observed: dict[tuple[str, int, str], dict[str, str]] = {}
    with (OUTPUT / "replicate_metrics.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            observed[(row["scenario_id"], int(row["replicate_index"]), row["method"])] = row

    mismatches: list[dict[str, object]] = []
    compared_rows = 0
    for (scenario_id, replicate_index), truth_map in sorted(truth.items()):
        for method in methods:
            key = (scenario_id, replicate_index, method)
            recalculated = calculate(truth_map, predictions.get(key, {}), directions)
            source = observed.get(key)
            if source is None:
                mismatches.append({"scenario_id": scenario_id, "replicate_index": replicate_index, "method": method, "field": "row", "expected": "present", "observed": "missing"})
                continue
            compared_rows += 1
            for field in METRICS + (
                "true_positive",
                "false_positive",
                "within_truth_false_positive",
                "true_negative",
                "false_negative",
            ):
                expected = recalculated[field]
                if field in METRICS and field != "out_of_truth_selected_pair_count":
                    actual = None if source[field] == "" else float(source[field])
                else:
                    actual = int(source[field])
                if actual != expected:
                    mismatches.append({
                        "scenario_id": scenario_id,
                        "replicate_index": replicate_index,
                        "method": method,
                        "field": field,
                        "expected": expected,
                        "observed": actual,
                    })

    report = {
        "schema_version": "1.0",
        "date": "2026-08-12",
        "status": "PASS" if not mismatches and len(truth) == expected_replicates and compared_rows == expected_replicates * len(methods) else "FAIL",
        "expected_replicates": expected_replicates,
        "truth_replicates": len(truth),
        "expected_metric_rows": expected_replicates * len(methods),
        "compared_metric_rows": compared_rows,
        "mismatch_count": len(mismatches),
        "mismatch_examples": mismatches[:20],
        "validator_uses_stored_truth_and_predictions": True,
        "metric_unit": "taxon-by-direction state",
        "runtime_memory_not_revalidated": True,
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    REPORT_MD.write_text(
        "# independent benchmark independent metric validation\n\n"
        f"Status: `{report['status']}`\n\n"
        f"Recomputed {compared_rows} method–replicate rows from stored truth and predictions; mismatches: {len(mismatches)}. Runtime and memory are descriptive run measurements and are not expected to be byte-reproducible.\n",
        encoding="utf-8",
    )
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
