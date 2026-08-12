from __future__ import annotations

import csv
from dataclasses import replace
import json
from pathlib import Path
import unittest

from provfold.benchmark import benchmark_factorial_predictions, benchmark_predictions
from provfold.compare import (
    FACTORIAL_METHODS,
    METHOD_DEFINITIONS,
    PROVFOLD_FACTORIAL_METHOD,
    compare_factorial_methods,
    compare_methods,
    factorial_method_key,
)
from provfold.core import aggregate
from provfold.config import AggregationConfig, SimulationConfig
from provfold.sensitivity import leave_one_family_out, threshold_surface
from provfold.simulate import simulate_records


ROOT = Path(__file__).resolve().parents[1]


def load_example():
    with (ROOT / "examples/minimal/input.csv").open(encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    config = AggregationConfig.from_dict(json.loads((ROOT / "examples/minimal/config.json").read_text()))
    return records, config


class AnalysisTests(unittest.TestCase):
    def test_all_comparators_are_present(self):
        records, config = load_example()
        result = compare_methods(records, config)
        self.assertEqual({row["method"] for row in result}, set(METHOD_DEFINITIONS))

    def test_all_eligible_comparator_scope_does_not_widen_full_method(self):
        records, config = load_example()
        primary = compare_methods(records, config)
        widened = compare_methods(records, replace(config, comparator_scope="all_eligible"))
        primary_full = [row for row in primary if row["method"] == "provenance_rank_folded"]
        widened_full = [row for row in widened if row["method"] == "provenance_rank_folded"]
        self.assertEqual(
            [{key: value for key, value in row.items() if key != "policy_hash"} for row in primary_full],
            [{key: value for key, value in row.items() if key != "policy_hash"} for row in widened_full],
        )
        primary_naive = [row for row in primary if row["method"] == "source_era_taxon"]
        widened_naive = [row for row in widened if row["method"] == "source_era_taxon"]
        self.assertNotEqual(
            [{key: value for key, value in row.items() if key != "policy_hash"} for row in primary_naive],
            [{key: value for key, value in row.items() if key != "policy_hash"} for row in widened_naive],
        )

    def test_complete_threshold_surface(self):
        records, config = load_example()
        result = threshold_surface(records, config)
        settings = {(row["minimum_support"], row["maximum_opposition"]) for row in result}
        self.assertEqual(settings, {(support, opposition) for support in range(1, 4) for opposition in range(0, 3)})

    def test_leave_every_primary_family_out(self):
        records, config = load_example()
        result = leave_one_family_out(records, config)
        self.assertEqual({row["omitted_family_id"] for row in result}, {"F001", "F002", "F003"})

    def test_simulation_is_deterministic_and_benchmarkable(self):
        simulation_config = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        first = simulate_records(simulation_config)
        second = simulate_records(simulation_config)
        self.assertEqual(first.simulation_hash, second.simulation_hash)
        self.assertEqual(first.truth, second.truth)
        self.assertEqual(first.records, second.records)
        aggregation_config = AggregationConfig.from_dict({
            "analysis_id": "simulation_test",
            "primary_comparisons": ["case_vs_control"],
            "primary_evidence_tiers": ["primary"],
            "direction_states": ["increased", "decreased"],
            "eligible_family_states": ["cohort_family"],
            "eligible_mapping_states": ["exact", "rank_folded"],
            "eligible_record_states": ["eligible"],
            "target_rank": "genus",
            "minimum_support": 2,
            "maximum_opposition": 0
        })
        metrics = benchmark_predictions(first.records, first.truth, aggregation_config)
        self.assertEqual({row["method"] for row in metrics}, set(METHOD_DEFINITIONS))
        for row in metrics:
            if row["precision"] is not None:
                self.assertGreaterEqual(row["precision"], 0.0)
                self.assertLessEqual(row["precision"], 1.0)

    def test_empty_selection_reports_undefined_precision(self):
        simulation_config = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        simulation = simulate_records(simulation_config)
        config = AggregationConfig.from_dict({
            "analysis_id": "empty_selection_test",
            "primary_comparisons": ["case_vs_control"],
            "primary_evidence_tiers": ["primary"],
            "direction_states": ["increased", "decreased"],
            "eligible_family_states": ["cohort_family"],
            "eligible_mapping_states": ["exact", "rank_folded"],
            "eligible_record_states": ["eligible"],
            "target_rank": "genus",
            "minimum_support": 999,
            "maximum_opposition": 0,
        })
        metrics = benchmark_factorial_predictions(simulation.records, simulation.truth, config)
        self.assertTrue(all(row["precision"] is None for row in metrics))
        self.assertTrue(all(not row["precision_defined"] for row in metrics))

    def test_truth_precedes_reporting_perturbations(self):
        base = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        stressed = replace(
            base,
            reports_per_family=4,
            selective_reporting_probability=0.2,
            direction_missing_probability=0.8,
            rank_inflation_probability=1.0,
            ambiguous_taxonomy_probability=0.8,
            comparator_leakage_probability=1.0,
            missing_family_probability=0.8,
            family_merge_probability=0.8,
            family_split_probability=0.8,
        )
        reference_result = simulate_records(base)
        stressed_result = simulate_records(stressed)
        self.assertEqual(reference_result.truth, stressed_result.truth)
        self.assertEqual(reference_result.cohort_effects, stressed_result.cohort_effects)
        self.assertNotEqual(reference_result.records, stressed_result.records)

    def test_missing_family_state_is_explicit(self):
        base = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        result = simulate_records(replace(base, missing_family_probability=1.0))
        self.assertTrue(result.records)
        self.assertTrue(all(not row["family_id"] and row["family_state"] == "indeterminate" for row in result.records))

    def test_factorial_ablation_contains_all_declared_cells(self):
        records, config = load_example()
        rows = compare_factorial_methods(records, config)
        self.assertEqual({row["method"] for row in rows}, set(FACTORIAL_METHODS))

    def test_full_factorial_cell_matches_aggregate(self):
        records, config = load_example()
        factorial = [
            row for row in compare_factorial_methods(records, config)
            if row["method"] == PROVFOLD_FACTORIAL_METHOD
        ]
        full = aggregate(records, config).taxon_summaries
        self.assertEqual(
            {(row["taxon_key"], row["selected_direction"], row["signal_classification"]) for row in factorial},
            {(row["target_taxon"], row["selected_direction"], row["signal_classification"]) for row in full},
        )

    def test_multi_primary_comparisons_exercise_family_state_collapse(self):
        base = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        simulation_config = replace(
            base,
            primary_comparisons_per_family=3,
            primary_comparison_fraction=0.7,
        )
        simulation = simulate_records(simulation_config)
        config = AggregationConfig.from_dict({
            "analysis_id": "multi_primary_test",
            "primary_comparisons": ["case_vs_control_panel_01", "case_vs_control_panel_02", "case_vs_control_panel_03"],
            "primary_evidence_tiers": ["primary"],
            "direction_states": ["increased", "decreased"],
            "eligible_family_states": ["cohort_family"],
            "eligible_mapping_states": ["exact", "rank_folded"],
            "eligible_record_states": ["eligible"],
            "target_rank": "genus",
            "minimum_support": 2,
            "maximum_opposition": 0,
        })
        rows = compare_factorial_methods(simulation.records, config)
        separate = {
            (row["taxon_key"], row["selected_direction"])
            for row in rows
            if row["method"] == factorial_method_key("family", "resolved_target", "separate_votes")
            and row["selected_direction"]
        }
        collapsed = {
            (row["taxon_key"], row["selected_direction"])
            for row in rows
            if row["method"] == PROVFOLD_FACTORIAL_METHOD
            and row["selected_direction"]
        }
        self.assertNotEqual(separate, collapsed)

    def test_factorial_benchmark_reports_shared_space_and_out_of_space_errors(self):
        base = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        simulation = simulate_records(replace(base, rank_inflation_probability=1.0))
        config = AggregationConfig.from_dict({
            "analysis_id": "factorial_metric_test",
            "primary_comparisons": ["case_vs_control"],
            "primary_evidence_tiers": ["primary"],
            "direction_states": ["increased", "decreased"],
            "eligible_family_states": ["cohort_family"],
            "eligible_mapping_states": ["exact", "rank_folded"],
            "eligible_record_states": ["eligible"],
            "target_rank": "genus",
            "minimum_support": 2,
            "maximum_opposition": 0,
        })
        metrics = benchmark_factorial_predictions(simulation.records, simulation.truth, config)
        self.assertEqual({row["method"] for row in metrics}, set(FACTORIAL_METHODS))
        self.assertTrue(all("false_discovery_proportion" in row for row in metrics))
        self.assertTrue(all("out_of_truth_selected_pair_count" in row for row in metrics))

    def test_rank_replacement_with_lineage_changes_source_key_not_resolved_key(self):
        base = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        reference = simulate_records(replace(base, rank_inflation_probability=0.0))
        replacement = simulate_records(replace(
            base,
            rank_inflation_probability=0.0,
            rank_replacement_probability=1.0,
            rank_lineage_missing_probability=0.0,
            ambiguous_taxonomy_probability=0.0,
        ))
        self.assertEqual(reference.truth, replacement.truth)
        self.assertEqual(reference.cohort_effects, replacement.cohort_effects)
        self.assertTrue(replacement.records)
        self.assertTrue(all(row["source_rank"] == "species" for row in replacement.records))
        self.assertTrue(all(row["target_taxon"] for row in replacement.records))

    def test_rank_replacement_without_lineage_is_non_voting(self):
        base = SimulationConfig.from_dict(json.loads((ROOT / "examples/minimal/simulation_config.json").read_text()))
        replacement = simulate_records(replace(
            base,
            rank_inflation_probability=0.0,
            rank_replacement_probability=1.0,
            rank_lineage_missing_probability=1.0,
            ambiguous_taxonomy_probability=0.0,
        ))
        self.assertTrue(replacement.records)
        self.assertTrue(all(not row["target_taxon"] and row["mapping_state"] == "ambiguous" for row in replacement.records))


if __name__ == "__main__":
    unittest.main()
