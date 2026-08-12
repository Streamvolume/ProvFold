from __future__ import annotations

from pathlib import Path
import unittest

from provfold.config import AggregationConfig
from provfold.core import aggregate
from provfold.io import read_csv, read_json


ROOT = Path(__file__).resolve().parents[1]


class AHCaseRegressionTests(unittest.TestCase):
    def test_ah_case_reference_signals_are_preserved(self):
        case = ROOT / "reproducibility/ah_case"
        records = read_csv(case / "adapted_input.csv")
        config = AggregationConfig.from_dict(read_json(case / "aggregation_config.json"))
        result = aggregate(records, config)
        retained = {
            (row["target_taxon"], row["selected_direction"])
            for row in result.taxon_summaries
            if row["signal_classification"] == "retained"
        }
        self.assertEqual(
            retained,
            {("Dialister", "decreased_in_AH"), ("Veillonella", "increased_in_AH")},
        )


if __name__ == "__main__":
    unittest.main()
