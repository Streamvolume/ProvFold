from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
import random
import unittest

from provfold.config import AggregationConfig
from provfold.core import aggregate, validate_records


ROOT = Path(__file__).resolve().parents[1]


def load_example():
    with (ROOT / "examples/minimal/input.csv").open(encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    import json
    config = AggregationConfig.from_dict(json.loads((ROOT / "examples/minimal/config.json").read_text()))
    return records, config


class CoreTests(unittest.TestCase):
    def test_validation_and_aggregation(self):
        records, config = load_example()
        report = validate_records(records, config)
        self.assertTrue(report.valid)
        result = aggregate(records, config)
        summaries = {row["target_taxon"]: row for row in result.taxon_summaries}
        self.assertEqual(summaries["Taxon_A"]["selected_direction"], "increased")
        self.assertEqual(summaries["Taxon_A"]["signal_classification"], "retained")
        self.assertEqual(summaries["Taxon_B"]["signal_classification"], "not_retained")
        self.assertEqual(summaries["Taxon_C"]["signal_classification"], "conflict_only")
        self.assertEqual(len(result.abstentions), 1)
        self.assertEqual(result.abstentions[0]["abstention_reason"], "family_state=indeterminate")
        self.assertEqual(len(result.record_map) + len(result.abstentions), len(records))

    def test_duplicate_rows_do_not_inflate_family_vote(self):
        records, config = load_example()
        result = aggregate(records, config)
        taxon_a_units = [row for row in result.evidence_units if row["target_taxon"] == "Taxon_A" and row["comparison_id"] == "case_vs_control"]
        self.assertEqual(len(taxon_a_units), 2)
        self.assertEqual(max(row["record_count"] for row in taxon_a_units), 2)

    def test_empty_lineage_component_value_is_permitted(self):
        records, config = load_example()
        records[0]["lineage_component"] = ""
        report = validate_records(records, config)
        self.assertTrue(report.valid)

    def test_lineage_component_field_is_optional(self):
        records, config = load_example()
        records[0].pop("lineage_component")
        report = validate_records(records, config)
        self.assertTrue(report.valid)

    def test_one_report_cannot_map_to_multiple_families(self):
        records, config = load_example()
        row = dict(records[0])
        row.update(record_id="R995", family_id="F999")
        report = validate_records(records + [row], config)
        self.assertFalse(report.valid)
        self.assertIn("report_maps_to_multiple_families", {item["code"] for item in report.issues})

    def test_row_order_invariance(self):
        records, config = load_example()
        original = aggregate(records, config)
        shuffled = list(records)
        random.Random(17).shuffle(shuffled)
        repeated = aggregate(shuffled, config)
        self.assertEqual(original.taxon_summaries, repeated.taxon_summaries)
        self.assertEqual(original.evidence_units, repeated.evidence_units)

    def test_duplicate_identifier_is_an_error(self):
        records, config = load_example()
        records.append(dict(records[0]))
        report = validate_records(records, config)
        self.assertFalse(report.valid)
        self.assertIn("duplicate_identifier", {row["code"] for row in report.issues})

    def test_higher_rank_abstains_without_downward_inference(self):
        records, config = load_example()
        row = dict(records[0])
        row.update(record_id="R999", source_rank="family", harmonized_rank="family", harmonized_taxon="Family_A", target_taxon="Taxon_A")
        result = aggregate([row], config)
        self.assertEqual(result.taxon_summaries, [])
        self.assertEqual(result.abstentions[0]["abstention_reason"], "higher_rank_cannot_map_down")

    def test_multiple_comparisons_do_not_create_multiple_family_votes(self):
        records, config = load_example()
        row = dict(records[0])
        row.update(record_id="R998", comparison_id="case_vs_control_panel2", report_id="P009", direction="increased")
        result = aggregate(records + [row], replace(config, primary_comparisons=("case_vs_control", "case_vs_control_panel2")))
        summary = next(item for item in result.taxon_summaries if item["target_taxon"] == "Taxon_A")
        self.assertEqual(summary["family_count"], 2)

    def test_opposite_panels_within_one_family_remain_conflict(self):
        records, config = load_example()
        row = dict(records[0])
        row.update(record_id="R997", comparison_id="case_vs_control_panel2", report_id="P009", direction="decreased")
        result = aggregate(records + [row], replace(config, primary_comparisons=("case_vs_control", "case_vs_control_panel2")))
        family_state = next(item for item in result.family_taxon_states if item["family_id"] == "F001" and item["target_taxon"] == "Taxon_A")
        self.assertEqual(family_state["family_direction_state"], "conflict")

    def test_one_family_id_cannot_mix_semantic_states(self):
        records, config = load_example()
        row = dict(records[0])
        row.update(record_id="R996", family_state="project_family")
        report = validate_records(records + [row], config)
        self.assertFalse(report.valid)
        self.assertIn("inconsistent_family_state", {item["code"] for item in report.issues})

    def test_non_primary_tier_cannot_veto_primary_evidence_unit(self):
        records, config = load_example()
        primary = [dict(records[0]), dict(records[2])]
        primary[0].update(record_id="MIX-001", direction="increased")
        primary[1].update(record_id="MIX-002", direction="increased")
        corroborative = dict(primary[1])
        corroborative.update(
            record_id="MIX-003",
            evidence_tier="corroborative",
            direction="decreased",
        )
        result = aggregate(primary + [corroborative], config)
        summary = next(row for row in result.taxon_summaries if row["target_taxon"] == "Taxon_A")
        self.assertEqual(summary["selected_direction"], "increased")
        excluded = next(row for row in result.record_map if row["record_id"] == "MIX-003")
        self.assertEqual(excluded["aggregation_status"], "excluded_from_primary_scope")
        self.assertEqual(excluded["non_voting_reason"], "evidence_tier_not_primary")


if __name__ == "__main__":
    unittest.main()
