#!/usr/bin/env python3
"""Audit what each factor in the factorial benchmark changes.

This audit reads the frozen prediction stream.  It does not regenerate data or
call the analytical implementation, so it provides a direct check on stored
selected sets and on the attribution statements used in the manuscript.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path


import os
BENCHMARK = Path(os.environ.get('PROVFOLD_BENCHMARK_OUTPUT', Path(__file__).resolve().parent/'generated_run')).resolve()
PREDICTIONS = BENCHMARK/'generated/predictions.jsonl.gz'
METRICS = BENCHMARK/'replicate_metrics.csv'
OUTPUT_CSV = BENCHMARK/'factorial_attribution_audit.csv'
OUTPUT_SUMMARY = BENCHMARK/'factorial_attribution_summary.csv'
OUTPUT_JSON = BENCHMARK/'factorial_attribution_validation.json'

FULL = "vu_family__taxon_resolved_target__comparison_family_state"
PLAIN = "vu_family__taxon_resolved_target__comparison_ignore"
SEPARATE = "vu_family__taxon_resolved_target__comparison_separate_votes"


def key(row: dict[str, object]) -> tuple[str, int]:
    return str(row["scenario_id"]), int(row["replicate_index"])


def main() -> None:
    selections: dict[tuple[str, int, str], set[tuple[str, str]]] = defaultdict(set)
    with gzip.open(PREDICTIONS, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["method"] not in {FULL, PLAIN, SEPARATE}:
                continue
            if row["signal_classification"] == "retained" and row["selected_direction"]:
                selections[(row["scenario_id"], row["replicate_index"], row["method"])].add(
                    (row["taxon_key"], row["selected_direction"])
                )

    metric_rows: dict[tuple[str, int, str], dict[str, str]] = {}
    with METRICS.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["method"] in {FULL, PLAIN, SEPARATE}:
                metric_rows[(row["scenario_id"], int(row["replicate_index"]), row["method"])] = row

    replicate_keys = sorted({(scenario, replicate) for scenario, replicate, method in metric_rows if method == FULL})
    rows: list[dict[str, object]] = []
    for scenario, replicate in replicate_keys:
        full_set = selections.get((scenario, replicate, FULL), set())
        plain_set = selections.get((scenario, replicate, PLAIN), set())
        separate_set = selections.get((scenario, replicate, SEPARATE), set())
        rows.append(
            {
                "scenario_id": scenario,
                "replicate_index": replicate,
                "full_equals_plain_family_resolved_ignore": full_set == plain_set,
                "full_equals_family_resolved_separate": full_set == separate_set,
                "full_selected_pair_count": len(full_set),
                "plain_selected_pair_count": len(plain_set),
                "separate_selected_pair_count": len(separate_set),
                "full_f1": metric_rows[(scenario, replicate, FULL)]["f1"],
                "plain_f1": metric_rows[(scenario, replicate, PLAIN)]["f1"],
                "separate_f1": metric_rows[(scenario, replicate, SEPARATE)]["f1"],
            }
        )

    with OUTPUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_scenario: dict[str, dict[str, int]] = defaultdict(lambda: {"replicates": 0, "full_plain_differences": 0, "full_separate_differences": 0})
    for row in rows:
        item = by_scenario[str(row["scenario_id"])]
        item["replicates"] += 1
        item["full_plain_differences"] += int(not row["full_equals_plain_family_resolved_ignore"])
        item["full_separate_differences"] += int(not row["full_equals_family_resolved_separate"])

    payload = {
        "status": "PASS",
        "replicate_count": len(rows),
        "full_vs_plain_family_resolved_ignore_differing_replicates": sum(
            int(not row["full_equals_plain_family_resolved_ignore"]) for row in rows
        ),
        "full_vs_family_resolved_separate_differing_replicates": sum(
            int(not row["full_equals_family_resolved_separate"]) for row in rows
        ),
        "multi_primary_full_vs_separate_differing_replicates": by_scenario["multiple_primary_comparisons"][
            "full_separate_differences"
        ],
        "interpretation": (
            "The full family-state method and the plain family-plus-resolved-target method have identical stored "
            "selected sets throughout this benchmark. Comparison handling is exercised only by the multi-primary "
            "scenario, where family-state collapse and separate comparison votes differ; this is a trade-off rather "
            "than evidence of uniform superiority."
        ),
        "by_scenario": dict(sorted(by_scenario.items())),
    }
    if payload["full_vs_plain_family_resolved_ignore_differing_replicates"] != 0:
        raise AssertionError("Full and plain family/resolved-target methods unexpectedly differ")
    with OUTPUT_SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["scenario_id", "replicates", "full_plain_differences", "full_separate_differences"],
        )
        writer.writeheader()
        for scenario, values in sorted(by_scenario.items()):
            writer.writerow({"scenario_id": scenario, **values})
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
