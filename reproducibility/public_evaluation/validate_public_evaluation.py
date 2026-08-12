#!/usr/bin/env python3
"""Independently validate public-resource from stored row-level derivatives."""

from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path


OUTPUT = Path(os.environ.get("PROVFOLD_PUBLIC_OUTPUT", Path(__file__).resolve().parent / "generated_run")).resolve()
GENERATED = OUTPUT / "generated"
REPORT_JSON = OUTPUT / "independent_validation.json"
REPORT_MD = OUTPUT / "independent_validation.md"
COMPLETION_JSON = OUTPUT / "completion_status.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def selected(rows: list[dict[str, str]], taxon_field: str) -> set[str]:
    return {f"{row[taxon_field]}|{row['selected_direction']}" for row in rows if row["selected_direction"]}


def sign(value: float) -> int:
    return 0 if value == 0 else (1 if value > 0 else -1)


def independently_aggregate(
    rows: list[dict[str, str]],
    minimum_support: int = 2,
    maximum_opposition: int = 0,
) -> tuple[
    dict[tuple[str, str, str], dict[str, object]],
    dict[tuple[str, str], dict[str, object]],
    dict[str, dict[str, object]],
]:
    """Recompute primary aggregation without importing or calling ProvFold."""

    voting = [
        row for row in rows
        if row["eligibility"] == "eligible"
        and row["mapping_state"] in {"direct_genus_ncbi", "folded_lower_rank_ncbi"}
    ]
    unit_rows: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in voting:
        unit_rows[(row["family_id"], row["comparison_id"], row["target_taxon"])].append(row)
    units: dict[tuple[str, str, str], dict[str, object]] = {}
    for key, members in unit_rows.items():
        directions = {row["direction"] for row in members}
        units[key] = {
            "direction_state": next(iter(directions)) if len(directions) == 1 else "conflict",
            "record_count": len(members),
        }

    family_units: dict[tuple[str, str], list[tuple[tuple[str, str, str], dict[str, object]]]] = defaultdict(list)
    for key, unit in units.items():
        family_units[(key[0], key[2])].append((key, unit))
    family_states: dict[tuple[str, str], dict[str, object]] = {}
    for key, members in family_units.items():
        states = {str(unit["direction_state"]) for _, unit in members}
        state = next(iter(states)) if len(states) == 1 and "conflict" not in states else "conflict"
        family_states[key] = {
            "family_direction_state": state,
            "comparison_count": len({unit_key[1] for unit_key, _ in members}),
            "evidence_unit_count": len(members),
        }

    taxon_states: dict[str, list[str]] = defaultdict(list)
    for (_, taxon), state in family_states.items():
        taxon_states[taxon].append(str(state["family_direction_state"]))
    summaries: dict[str, dict[str, object]] = {}
    for taxon, states in taxon_states.items():
        increased = states.count("increased")
        decreased = states.count("decreased")
        conflicts = states.count("conflict")
        selected_direction = ""
        for direction, support, opposition in (
            ("increased", increased, decreased),
            ("decreased", decreased, increased),
        ):
            if support >= minimum_support and opposition <= maximum_opposition and support > opposition:
                if selected_direction:
                    selected_direction = ""
                    break
                selected_direction = direction
        summaries[taxon] = {
            "increased": increased,
            "decreased": decreased,
            "conflicts": conflicts,
            "selected_direction": selected_direction,
        }
    return units, family_states, summaries


def main() -> None:
    errors: list[str] = []
    manifest = read_csv(OUTPUT / "artifact_manifest.csv")
    for row in manifest:
        path = OUTPUT / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
            errors.append(f"artifact identity mismatch: {row['relative_path']}")

    result = json.loads((OUTPUT / "results.json").read_text(encoding="utf-8"))
    if result["arms_pooled"] or result["publication_families_called_cohort_families"]:
        errors.append("an asymmetric-arm semantic prohibition was violated")
    if result["taxonomy"] != {
        "alias_resolved_count": 4,
        "genus_lineage_available_count": 354,
        "query_count": 503,
        "resolved_count": 503,
    }:
        errors.append("taxonomy identity counts differ from the frozen official mapping")
    unresolved = result["unresolved_literature"]
    if (unresolved["comparison_count"], unresolved["publication_count"], unresolved["genus_direction_row_count"], unresolved["quantitative_pool_count"]) != (16, 14, 57, 0):
        errors.append("unresolved-identifier audit counts differ from the frozen corpus contract")
    disbiome = result["disbiome"]
    if (disbiome["strict_association_rows"], disbiome["strict_publications"], disbiome["experiment_level_explicit_family_rows"], disbiome["quantitative_recurrence_rows"]) != (5011, 547, 0, 0):
        errors.append("Disbiome determinability counts differ from the frozen corpus contract")

    records = read_csv(GENERATED / "adapted_records.csv")
    evidence_units = read_csv(GENERATED / "evidence_units.csv")
    family_states = read_csv(GENERATED / "family_taxon_states.csv")
    taxon_summaries = read_csv(GENERATED / "taxon_summaries.csv")
    methods = read_csv(GENERATED / "method_comparison.csv")
    surfaces = read_csv(GENERATED / "threshold_surface.csv")
    omissions = read_csv(GENERATED / "leave_one_family_out.csv")
    disease_summary = read_csv(OUTPUT / "arm_disease_summary.csv")
    method_summary = read_csv(OUTPUT / "primary_method_summary.csv")

    quantitative_records = [row for row in records if row["pool_status"] == "quantitative"]
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in quantitative_records:
        groups[(row["arm"], row["disease_id"])].append(row)
    if len(groups) != 17 or len(disease_summary) != 17:
        errors.append("expected 17 non-pooled arm-disease groups")
    summary_by_key = {(row["arm"], row["disease_id"]): row for row in disease_summary}
    for key, rows in groups.items():
        summary = summary_by_key.get(key)
        if summary is None:
            errors.append(f"missing disease summary: {key}")
            continue
        voting = [row for row in rows if row["eligibility"] == "eligible" and row["mapping_state"] in {"direct_genus_ncbi", "folded_lower_rank_ncbi"}]
        families = {row["family_id"] for row in voting}
        checks = {
            "directional_source_row_count": len(rows),
            "voting_row_count": len(voting),
            "abstention_row_count": len(rows) - len(voting),
            "evaluable_family_count": len(families),
        }
        for field, value in checks.items():
            if int(summary[field]) != value:
                errors.append(f"{key} {field} mismatch")
        if summary["evaluation_status"] == "evaluable":
            independent_units, independent_states, independent_summaries = independently_aggregate(rows)
            stored_units = {
                (row["family_id"], row["comparison_id"], row["target_taxon"]): {
                    "direction_state": row["direction_state"],
                    "record_count": int(row["record_count"]),
                }
                for row in evidence_units if (row["arm"], row["disease_id"]) == key
            }
            if independent_units != stored_units:
                errors.append(f"{key} independent evidence-unit recomputation mismatch")
            stored_states = {
                (row["family_id"], row["target_taxon"]): {
                    "family_direction_state": row["family_direction_state"],
                    "comparison_count": int(row["comparison_count"]),
                    "evidence_unit_count": int(row["evidence_unit_count"]),
                }
                for row in family_states if (row["arm"], row["disease_id"]) == key
            }
            if independent_states != stored_states:
                errors.append(f"{key} independent family-state recomputation mismatch")
            target_rows = [row for row in taxon_summaries if (row["arm"], row["disease_id"]) == key]
            selected_pairs = selected(target_rows, "target_taxon")
            if int(summary["primary_selected_pair_count"]) != len(selected_pairs):
                errors.append(f"{key} selected-pair count mismatch")
            stored_summaries = {
                row["target_taxon"]: {
                    "increased": int(json.loads(row["direction_counts_json"])["increased"]),
                    "decreased": int(json.loads(row["direction_counts_json"])["decreased"]),
                    "conflicts": int(row["conflict_family_count"]),
                    "selected_direction": row["selected_direction"],
                }
                for row in target_rows
            }
            if independent_summaries != stored_summaries:
                errors.append(f"{key} independent primary-summary recomputation mismatch")
            surface_rows = [row for row in surfaces if (row["arm"], row["disease_id"]) == key]
            settings = {(int(row["minimum_support"]), int(row["maximum_opposition"])) for row in surface_rows}
            if len(settings) != len(families) ** 2 or int(summary["threshold_setting_count"]) != len(settings):
                errors.append(f"{key} threshold surface is incomplete")
            for minimum_support, maximum_opposition in settings:
                _, _, independent_setting = independently_aggregate(rows, minimum_support, maximum_opposition)
                stored_setting = {
                    row["target_taxon"]: row["selected_direction"]
                    for row in surface_rows
                    if int(row["minimum_support"]) == minimum_support
                    and int(row["maximum_opposition"]) == maximum_opposition
                }
                expected_setting = {taxon: value["selected_direction"] for taxon, value in independent_setting.items()}
                if stored_setting != expected_setting:
                    errors.append(f"{key} independent threshold recomputation mismatch at {(minimum_support, maximum_opposition)}")
            method_rows = [row for row in method_summary if (row["arm"], row["disease_id"]) == key]
            if len(method_rows) != 5:
                errors.append(f"{key} does not have all five method summaries")
            detailed_method_rows = [row for row in methods if (row["arm"], row["disease_id"]) == key]
            for method_row in method_rows:
                subset = [row for row in detailed_method_rows if row["method"] == method_row["method"]]
                if int(method_row["selected_pair_count"]) != len(selected(subset, "taxon_key")):
                    errors.append(f"{key} {method_row['method']} selected count mismatch")
            omitted = {row["omitted_family_id"] for row in omissions if (row["arm"], row["disease_id"]) == key}
            if omitted != families:
                errors.append(f"{key} family omission coverage mismatch")
            for omitted_family in families:
                _, _, independent_omission = independently_aggregate(
                    [row for row in rows if row["family_id"] != omitted_family]
                )
                expected_pairs = {
                    f"{taxon}|{value['selected_direction']}"
                    for taxon, value in independent_omission.items() if value["selected_direction"]
                }
                stored_omission = selected([
                    row for row in omissions
                    if (row["arm"], row["disease_id"]) == key
                    and row["omitted_family_id"] == omitted_family
                ], "target_taxon")
                if expected_pairs != stored_omission:
                    errors.append(f"{key} independent omission recomputation mismatch for {omitted_family}")

    balance_rows = read_csv(OUTPUT / "cross_arm_direction_balance.csv")
    balance_by_key = {(row["disease_id"], row["genus"]): row for row in balance_rows}
    states_by_key: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in family_states:
        states_by_key[(row["arm"], row["disease_id"], row["target_taxon"])].append(row)
    overlap_ids = {"D003093", "D003424", "D010300", "D015179"}
    for disease_id in overlap_ids:
        raw_genera = {key[2] for key in states_by_key if key[:2] == ("rawdata", disease_id)}
        literature_genera = {key[2] for key in states_by_key if key[:2] == ("literature", disease_id)}
        for genus in raw_genera & literature_genera:
            stored = balance_by_key.get((disease_id, genus))
            if stored is None:
                errors.append(f"missing cross-arm balance row: {disease_id} {genus}")
                continue
            values = []
            for arm in ("rawdata", "literature"):
                states = states_by_key[(arm, disease_id, genus)]
                inc = sum(row["family_direction_state"] == "increased" for row in states)
                dec = sum(row["family_direction_state"] == "decreased" for row in states)
                values.append((inc - dec) / (inc + dec) if inc + dec else None)
            if values[0] is not None and abs(float(stored["rawdata_direction_balance"]) - values[0]) > 1e-8:
                errors.append(f"rawdata balance mismatch: {disease_id} {genus}")
            if values[1] is not None and abs(float(stored["literature_direction_balance"]) - values[1]) > 1e-8:
                errors.append(f"literature balance mismatch: {disease_id} {genus}")
            expected = "excluded_zero_balance" if values[0] is None or values[1] is None or sign(values[0]) == 0 or sign(values[1]) == 0 else ("agree" if sign(values[0]) == sign(values[1]) else "disagree")
            if stored["sign_agreement_status"] != expected:
                errors.append(f"sign agreement mismatch: {disease_id} {genus}")

    jaccard_rows = read_csv(OUTPUT / "cross_arm_jaccard_surface.csv")
    surface_by_key: dict[tuple[str, str, int, int], list[dict[str, str]]] = defaultdict(list)
    for row in surfaces:
        surface_by_key[(row["arm"], row["disease_id"], int(row["minimum_support"]), int(row["maximum_opposition"]))].append(row)
    for row in jaccard_rows:
        setting = (int(row["minimum_support"]), int(row["maximum_opposition"]))
        first = selected(surface_by_key[("rawdata", row["disease_id"], *setting)], "target_taxon")
        second = selected(surface_by_key[("literature", row["disease_id"], *setting)], "target_taxon")
        union = first | second
        expected = "" if not union else round(len(first & second) / len(union), 8)
        observed = row["Jaccard"]
        if (expected == "" and observed != "") or (expected != "" and abs(float(observed) - expected) > 1e-8):
            errors.append(f"cross-arm Jaccard mismatch: {row['disease_id']} {setting}")

    report = {
        "schema_version": "1.0",
        "date": "2026-08-12",
        "status": "PASS" if not errors else "FAIL",
        "error_count": len(errors),
        "errors": errors[:100],
        "artifact_identity_count": len(manifest),
        "adapted_record_count": len(records),
        "independent_primary_aggregation_group_count": sum(
            row["evaluation_status"] == "evaluable" for row in disease_summary
        ),
        "arm_disease_group_count": len(groups),
        "family_taxon_state_count": len(family_states),
        "threshold_row_count": len(surfaces),
        "omission_row_count": len(omissions),
        "cross_arm_balance_row_count": len(balance_rows),
        "cross_arm_jaccard_row_count": len(jaccard_rows),
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    REPORT_MD.write_text(
        "# public-resource independent validation\n\n"
        f"Status: `{report['status']}`\n\n"
        f"Validated {report['adapted_record_count']} derived records across {report['arm_disease_group_count']} non-pooled arm-disease groups, "
        f"including independent reconstruction of evidence units, family states, primary summaries, all threshold settings and every family omission "
        f"for {report['independent_primary_aggregation_group_count']} evaluable groups. The audit covered {report['threshold_row_count']} threshold rows, "
        f"{report['omission_row_count']} omission rows and {report['artifact_identity_count']} artefact identities. Errors: {report['error_count']}.\n",
        encoding="utf-8",
    )
    completion = {
        "schema_version": "1.0",
        "date": "2026-08-12",
        "status": "COMPLETE_INDEPENDENTLY_VALIDATED" if not errors else "VALIDATION_FAILED",
        "contract_sha256": result["contract_sha256"],
        "results_sha256": sha256(OUTPUT / "results.json"),
        "artifact_manifest_sha256": sha256(OUTPUT / "artifact_manifest.csv"),
        "independent_validation_sha256": sha256(REPORT_JSON),
        "independent_validation_report_sha256": sha256(REPORT_MD),
        "error_count": len(errors),
    }
    COMPLETION_JSON.write_text(json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if errors:
        raise SystemExit("; ".join(errors[:20]))


if __name__ == "__main__":
    main()
