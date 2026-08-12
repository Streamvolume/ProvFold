#!/usr/bin/env python3
"""Run candidate-contained integrity, unit, example and empirical-case checks."""

from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from provfold.compare import compare_factorial_methods
from provfold.io import load_aggregation_config, read_csv as read_provfold_csv


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def compare_expected_csv(observed: Path, expected: Path) -> list[str]:
    """Return field-level differences; an empty list means exact CSV equality."""
    if not observed.is_file() or not expected.is_file():
        return [f"missing comparison file: {observed.name} or {expected.name}"]
    observed_rows = rows(observed)
    expected_rows = rows(expected)
    differences: list[str] = []
    if observed_rows and expected_rows and list(observed_rows[0]) != list(expected_rows[0]):
        differences.append(f"header differs: {observed.name}")
    if len(observed_rows) != len(expected_rows):
        differences.append(f"row count differs for {observed.name}: {len(observed_rows)} != {len(expected_rows)}")
    for index, (left, right) in enumerate(zip(observed_rows, expected_rows), start=2):
        if left != right:
            fields = sorted(field for field in set(left) | set(right) if left.get(field, "") != right.get(field, ""))
            differences.append(f"{observed.name} row {index} differs in: {','.join(fields)}")
            if len(differences) >= 20:
                differences.append("difference report truncated")
                break
    return differences


def compare_expected_rows(observed: list[dict[str, object]], expected: Path) -> list[str]:
    """Compare in-memory rows with a CSV using the CSV writer's string semantics."""
    expected_rows = rows(expected)
    observed_rows = [
        {key: "" if value is None else str(value) for key, value in row.items()}
        for row in observed
    ]
    differences: list[str] = []
    if len(observed_rows) != len(expected_rows):
        differences.append(
            f"row count differs for {expected.name}: {len(observed_rows)} != {len(expected_rows)}"
        )
    for index, (left, right) in enumerate(zip(observed_rows, expected_rows), start=2):
        if left != right:
            fields = sorted(field for field in set(left) | set(right) if left.get(field, "") != right.get(field, ""))
            differences.append(f"{expected.name} row {index} differs in: {','.join(fields)}")
            if len(differences) >= 20:
                differences.append("difference report truncated")
                break
    return differences


def selected(items: list[dict[str, str]], taxon_field: str) -> set[str]:
    return {
        f"{row[taxon_field]}|{row['selected_direction']}"
        for row in items if row.get("selected_direction")
    }


def main() -> None:
    errors: list[str] = []
    for row in rows(ROOT / "SHA256SUMS.csv"):
        path = ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
            errors.append(f"manifest mismatch: {row['path']}")

    loader = unittest.TestLoader()
    suite = loader.discover(str(ROOT / "tests"))
    test_count = suite.countTestCases()
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    if test_count != 28 or not result.wasSuccessful():
        errors.append(f"unit-test failure or count mismatch: {test_count}")

    with tempfile.TemporaryDirectory(prefix="provfold-release-check-") as temp:
        temporary = Path(temp)
        commands = [
            [sys.executable, "-m", "provfold.cli", "reproduce", "--recipe", str(ROOT / "examples/minimal/recipe.json"), "--output-dir", str(temporary / "minimal")],
            [sys.executable, "-m", "provfold.cli", "reproduce", "--recipe", str(ROOT / "reproducibility/ah_case/recipe.json"), "--output-dir", str(temporary / "ah")],
            [sys.executable, "-m", "provfold.cli", "reproduce", "--recipe", str(ROOT / "reproducibility/hcc_repeated_report_case/recipe.json"), "--output-dir", str(temporary / "hcc")],
        ]
        for command in commands:
            completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
            if completed.returncode:
                errors.append(f"command failed: {' '.join(command)}: {completed.stderr.strip()}")

        alpha_commands = [
            [sys.executable, str(ROOT / "reproducibility/ah_alpha_diversity/recalculate_alpha_diversity.py")],
            [sys.executable, str(ROOT / "reproducibility/ah_alpha_diversity/validate_alpha_diversity_independently.py")],
        ]
        for command in alpha_commands:
            completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
            if completed.returncode:
                errors.append(f"AH quantitative command failed: {' '.join(command)}: {completed.stderr.strip()}")

        ah = temporary / "ah"
        required_ah = [
            ah / "aggregate/taxon_summaries.csv",
            ah / "sensitivity/threshold_surface.csv",
            ah / "omit_family/leave_one_family_out.csv",
        ]
        if not all(path.is_file() for path in required_ah):
            errors.append("AH-case reproduction did not create all required outputs")
        else:
            expected_pairs = [
                (ah / "compare/method_comparison.csv", ROOT / "reproducibility/ah_case/expected/method_comparison.csv"),
                (ah / "sensitivity/threshold_surface.csv", ROOT / "reproducibility/ah_case/expected/threshold_surface.csv"),
                (ah / "omit_family/leave_one_family_out.csv", ROOT / "reproducibility/ah_case/expected/leave_one_family_out.csv"),
            ]
            exact_differences = [
                difference
                for observed, expected in expected_pairs
                for difference in compare_expected_csv(observed, expected)
            ]
            ah_factorial = compare_factorial_methods(
                read_provfold_csv(ROOT / "reproducibility/ah_case/adapted_input.csv"),
                load_aggregation_config(ROOT / "reproducibility/ah_case/aggregation_config.json"),
            )
            exact_differences.extend(
                compare_expected_rows(
                    ah_factorial,
                    ROOT / "reproducibility/ah_case/expected/factorial_method_comparison.csv",
                )
            )
            errors.extend(f"AH-case expected mismatch: {item}" for item in exact_differences)
            primary = selected(rows(required_ah[0]), "target_taxon")
            expected_primary = {"Dialister|decreased_in_AH", "Veillonella|increased_in_AH"}
            if primary != expected_primary:
                errors.append(f"AH-case primary signals differ: {sorted(primary)}")
            by_setting: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
            for row in rows(required_ah[1]):
                by_setting[(int(row["minimum_support"]), int(row["maximum_opposition"]))].append(row)
            if len(selected(by_setting[(1, 0)], "target_taxon")) != 18:
                errors.append("AH-case support-1 count differs from 18")
            if len(selected(by_setting[(2, 0)], "target_taxon")) != 2:
                errors.append("AH-case support-2 count differs from 2")
            by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in rows(required_ah[2]):
                by_family[row["omitted_family_id"]].append(row)
            if sorted(by_family) != ["COH0003", "COH0004"] or any(selected(group, "target_taxon") for group in by_family.values()):
                errors.append("AH-case provenance-group omission no longer removes both signals")

        hcc = temporary / "hcc"
        required_hcc = [
            hcc / "aggregate/evidence_units.csv",
            hcc / "aggregate/taxon_summaries.csv",
            hcc / "compare/method_comparison.csv",
            hcc / "sensitivity/threshold_surface.csv",
            hcc / "omit_family/leave_one_family_out.csv",
        ]
        if not all(path.is_file() for path in required_hcc):
            errors.append("HCC stress-test reproduction did not create all required outputs")
        else:
            expected_pairs = [
                (required_hcc[2], ROOT / "reproducibility/hcc_repeated_report_case/expected/method_comparison.csv"),
                (required_hcc[3], ROOT / "reproducibility/hcc_repeated_report_case/expected/threshold_surface.csv"),
                (required_hcc[4], ROOT / "reproducibility/hcc_repeated_report_case/expected/leave_one_family_out.csv"),
            ]
            exact_differences = [
                difference
                for observed, expected in expected_pairs
                for difference in compare_expected_csv(observed, expected)
            ]
            hcc_factorial = compare_factorial_methods(
                read_provfold_csv(ROOT / "reproducibility/hcc_repeated_report_case/input.csv"),
                load_aggregation_config(ROOT / "reproducibility/hcc_repeated_report_case/aggregation_config.json"),
            )
            exact_differences.extend(
                compare_expected_rows(
                    hcc_factorial,
                    ROOT / "reproducibility/hcc_repeated_report_case/expected/factorial_method_comparison.csv",
                )
            )
            errors.extend(f"HCC expected mismatch: {item}" for item in exact_differences)
            method_rows = rows(required_hcc[2])
            report_selected = selected(
                [row for row in method_rows if row["method"] == "unique_publication"],
                "taxon_key",
            )
            provenance_selected = selected(
                [row for row in method_rows if row["method"] == "provenance_rank_folded"],
                "taxon_key",
            )
            expected_report = {
                "Blautia|increased_in_hcc",
                "Megamonas|decreased_in_hcc",
                "Prevotella|increased_in_hcc",
            }
            expected_provenance = {
                "Blautia|increased_in_hcc",
                "Prevotella|increased_in_hcc",
            }
            if report_selected != expected_report or provenance_selected != expected_provenance:
                errors.append("HCC report-versus-family selected sets differ from the frozen result")
            megamonas_units = [
                row for row in rows(required_hcc[0])
                if row["target_taxon"] == "Megamonas"
                and row["comparison_id"] == "non_hbv_hcc_vs_healthy"
            ]
            if len(megamonas_units) != 1 or megamonas_units[0]["report_count"] != "2":
                errors.append("HCC Megamonas reports no longer collapse into one evidence unit")

    payload = {
        "status": "PASS" if not errors else "FAIL",
        "manifest_entries": len(rows(ROOT / "SHA256SUMS.csv")),
        "unit_test_count": test_count,
        "ah_case_reference_signal_count": 2,
        "ah_case_expected_csv_exact_count": 4,
        "hcc_report_counting_signal_count": 3,
        "hcc_provenance_signal_count": 2,
        "hcc_expected_csv_exact_count": 4,
        "errors": errors,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
