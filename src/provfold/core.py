"""Core validation and provenance-aware evidence aggregation."""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from typing import Iterable, Mapping

from .config import AggregationConfig
from .models import AggregationResult, ValidationReport


REQUIRED_INPUT_FIELDS = (
    "record_id",
    "report_id",
    "family_id",
    "family_state",
    "comparison_id",
    "evidence_tier",
    "source_taxon",
    "source_rank",
    "harmonized_taxon",
    "harmonized_rank",
    "target_taxon",
    "mapping_state",
    "direction",
    "eligibility",
    "source_location",
)


def as_text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _issue(level: str, code: str, record_id: str, field: str, message: str) -> dict[str, str]:
    return {
        "level": level,
        "code": code,
        "record_id": record_id,
        "field": field,
        "message": message,
    }


def validate_records(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> ValidationReport:
    """Validate identifiers, required states and conservative rank-folding inputs."""

    rows = [dict(row) for row in records]
    issues: list[dict[str, str]] = []
    identifiers: list[str] = []
    family_states_by_id: dict[str, set[str]] = defaultdict(set)
    family_ids_by_report: dict[str, set[str]] = defaultdict(set)
    rank_index = {rank: index for index, rank in enumerate(config.rank_order)}

    for position, row in enumerate(rows, start=1):
        record_id = as_text(row.get("record_id")) or f"row:{position}"
        missing = [field for field in REQUIRED_INPUT_FIELDS if field not in row]
        for field in missing:
            issues.append(_issue("error", "missing_field", record_id, field, "required field is absent"))
        if missing:
            continue

        if not as_text(row["record_id"]):
            issues.append(_issue("error", "empty_identifier", record_id, "record_id", "record_id is empty"))
        else:
            identifiers.append(as_text(row["record_id"]))

        for field in ("report_id", "comparison_id", "evidence_tier", "source_taxon", "source_rank", "source_location"):
            if not as_text(row[field]):
                issues.append(_issue("error", "empty_required_value", record_id, field, f"{field} is empty"))

        family_state = as_text(row["family_state"])
        family_id = as_text(row["family_id"])
        if family_state in config.eligible_family_states and not family_id:
            issues.append(_issue("error", "eligible_family_without_id", record_id, "family_id", "an eligible family state requires family_id"))
        if not family_state:
            issues.append(_issue("error", "empty_required_value", record_id, "family_state", "family_state is empty"))
        if family_id and family_state:
            family_states_by_id[family_id].add(family_state)
        report_id = as_text(row["report_id"])
        if report_id and family_id:
            family_ids_by_report[report_id].add(family_id)

        eligibility = as_text(row["eligibility"])
        if not eligibility:
            issues.append(_issue("error", "empty_required_value", record_id, "eligibility", "eligibility is empty"))

        direction = as_text(row["direction"])
        if eligibility in config.eligible_record_states and direction and direction not in config.direction_states:
            issues.append(_issue("error", "unsupported_direction", record_id, "direction", f"direction is not one of {config.direction_states}"))

        harmonized_rank = as_text(row["harmonized_rank"])
        if harmonized_rank and harmonized_rank not in rank_index:
            issues.append(_issue("error", "unsupported_rank", record_id, "harmonized_rank", "harmonized_rank is absent from rank_order"))

        if eligibility in config.eligible_record_states and as_text(row["mapping_state"]) in config.eligible_mapping_states:
            if not as_text(row["harmonized_taxon"]) or not harmonized_rank:
                issues.append(_issue("error", "eligible_record_without_harmonized_taxon", record_id, "harmonized_taxon", "eligible mapped record lacks a harmonized taxon and rank"))
            elif harmonized_rank in rank_index:
                source_position = rank_index[harmonized_rank]
                target_position = rank_index[config.target_rank]
                if (
                    source_position == target_position
                    and as_text(row["target_taxon"])
                    and as_text(row["target_taxon"]) != as_text(row["harmonized_taxon"])
                ):
                    issues.append(_issue("error", "target_taxon_mismatch", record_id, "target_taxon", "target-rank taxon differs from the harmonized taxon"))
                if source_position < target_position:
                    issues.append(_issue("warning", "higher_rank_cannot_map_down", record_id, "target_taxon", "record will abstain because higher ranks cannot map down to the target rank"))
                elif source_position > target_position and not as_text(row["target_taxon"]):
                    issues.append(_issue("warning", "lower_rank_without_target_lineage", record_id, "target_taxon", "record will abstain because no frozen target-rank lineage is supplied"))

    duplicate_ids = sorted(item for item, count in Counter(identifiers).items() if count > 1)
    for record_id in duplicate_ids:
        issues.append(_issue("error", "duplicate_identifier", record_id, "record_id", "record_id occurs more than once"))
    for family_id, states in sorted(family_states_by_id.items()):
        if len(states) > 1:
            issues.append(_issue("error", "inconsistent_family_state", family_id, "family_state", f"family_id has multiple states: {', '.join(sorted(states))}"))
    for report_id, family_ids in sorted(family_ids_by_report.items()):
        if len(family_ids) > 1:
            issues.append(_issue(
                "error",
                "report_maps_to_multiple_families",
                report_id,
                "family_id",
                f"report_id maps to multiple family_id values: {', '.join(sorted(family_ids))}",
            ))

    issues.sort(key=lambda item: (item["level"], item["record_id"], item["field"], item["code"]))
    return ValidationReport(
        valid=not any(item["level"] == "error" for item in issues),
        record_count=len(rows),
        issues=issues,
    )


def resolve_target_taxon(row: Mapping[str, object], config: AggregationConfig) -> tuple[str, str]:
    """Return a frozen target-rank taxon or an explicit abstention reason."""

    if as_text(row["eligibility"]) not in config.eligible_record_states:
        return "", f"eligibility={as_text(row['eligibility']) or 'missing'}"
    if as_text(row["family_state"]) not in config.eligible_family_states:
        return "", f"family_state={as_text(row['family_state']) or 'missing'}"
    if not as_text(row["family_id"]):
        return "", "missing_family_id"
    if as_text(row["mapping_state"]) not in config.eligible_mapping_states:
        return "", f"mapping_state={as_text(row['mapping_state']) or 'missing'}"
    if as_text(row["direction"]) not in config.direction_states:
        return "", f"direction={as_text(row['direction']) or 'missing'}"

    rank_index = {rank: index for index, rank in enumerate(config.rank_order)}
    rank = as_text(row["harmonized_rank"])
    if rank not in rank_index:
        return "", f"unsupported_harmonized_rank={rank or 'missing'}"
    source_position = rank_index[rank]
    target_position = rank_index[config.target_rank]
    if source_position < target_position:
        return "", "higher_rank_cannot_map_down"
    if source_position == target_position:
        taxon = as_text(row["target_taxon"]) or as_text(row["harmonized_taxon"])
        return (taxon, "") if taxon else ("", "missing_target_rank_taxon")
    taxon = as_text(row["target_taxon"])
    return (taxon, "") if taxon else ("", "lower_rank_without_target_lineage")


def prepare_records(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> tuple[list[dict[str, object]], list[dict[str, object]], ValidationReport]:
    """Validate and separate voting records from explicit abstentions."""

    rows = [dict(row) for row in records]
    validation = validate_records(rows, config)
    if not validation.valid:
        messages = [f"{item['record_id']}:{item['field']}:{item['code']}" for item in validation.issues if item["level"] == "error"]
        raise ValueError("input validation failed: " + "; ".join(messages))

    eligible: list[dict[str, object]] = []
    abstentions: list[dict[str, object]] = []
    for row in rows:
        target_taxon, reason = resolve_target_taxon(row, config)
        if reason:
            abstentions.append({
                "record_id": as_text(row["record_id"]),
                "report_id": as_text(row["report_id"]),
                "family_id": as_text(row["family_id"]),
                "family_state": as_text(row["family_state"]),
                "comparison_id": as_text(row["comparison_id"]),
                "source_taxon": as_text(row["source_taxon"]),
                "harmonized_taxon": as_text(row["harmonized_taxon"]),
                "mapping_state": as_text(row["mapping_state"]),
                "direction": as_text(row["direction"]),
                "eligibility": as_text(row["eligibility"]),
                "abstention_reason": reason,
                "source_location": as_text(row["source_location"]),
                "policy_hash": config.policy_hash,
            })
            continue
        prepared = dict(row)
        prepared["resolved_target_taxon"] = target_taxon
        eligible.append(prepared)
    return eligible, sorted(abstentions, key=lambda item: str(item["record_id"])), validation


def evaluate_direction_counts(
    direction_counts: Mapping[str, int],
    conflict_count: int,
    config: AggregationConfig,
) -> tuple[str, str, dict[str, dict[str, object]]]:
    """Evaluate every direction under the same support/opposition rule."""

    evaluations: dict[str, dict[str, object]] = {}
    passed: list[str] = []
    for direction in config.direction_states:
        support = int(direction_counts.get(direction, 0))
        opposition = sum(int(direction_counts.get(other, 0)) for other in config.direction_states if other != direction)
        if config.conflict_policy == "opposition":
            opposition += conflict_count
        retained = (
            support >= config.minimum_support
            and opposition <= config.maximum_opposition
            and support > opposition
        )
        evaluations[direction] = {
            "support": support,
            "opposition": opposition,
            "conflicts": conflict_count,
            "passes": retained,
        }
        if retained:
            passed.append(direction)
    if len(passed) == 1:
        return passed[0], "retained", evaluations
    if len(passed) > 1:
        raise RuntimeError("more than one direction passed under a binary strict-majority configuration")
    if sum(direction_counts.values()) == 0 and conflict_count == 0:
        return "", "not_evaluable", evaluations
    if conflict_count and not sum(direction_counts.values()):
        return "", "conflict_only", evaluations
    return "", "not_retained", evaluations


def aggregate(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> AggregationResult:
    """Collapse records and assess recurrence without row or panel inflation."""

    eligible, abstentions, _ = prepare_records(records, config)
    primary_eligible = [
        row for row in eligible
        if as_text(row["comparison_id"]) in config.primary_comparisons
        and as_text(row["evidence_tier"]) in config.primary_evidence_tiers
    ]
    excluded_from_primary = [
        row for row in eligible
        if as_text(row["comparison_id"]) not in config.primary_comparisons
        or as_text(row["evidence_tier"]) not in config.primary_evidence_tiers
    ]
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in primary_eligible:
        key = (
            as_text(row["family_id"]),
            as_text(row["comparison_id"]),
            as_text(row["resolved_target_taxon"]),
        )
        grouped[key].append(row)

    record_map: list[dict[str, object]] = [
        {
            "record_id": as_text(row["record_id"]),
            "evidence_unit_id": "",
            "family_id": as_text(row["family_id"]),
            "comparison_id": as_text(row["comparison_id"]),
            "target_taxon": as_text(row["resolved_target_taxon"]),
            "direction": as_text(row["direction"]),
            "mapping_state": as_text(row["mapping_state"]),
            "aggregation_status": "excluded_from_primary_scope",
            "non_voting_reason": ";".join(
                reason
                for reason, excluded in (
                    ("comparison_not_primary", as_text(row["comparison_id"]) not in config.primary_comparisons),
                    ("evidence_tier_not_primary", as_text(row["evidence_tier"]) not in config.primary_evidence_tiers),
                )
                if excluded
            ),
            "policy_hash": config.policy_hash,
        }
        for row in excluded_from_primary
    ]
    evidence_units: list[dict[str, object]] = []
    for sequence, key in enumerate(sorted(grouped), start=1):
        rows = grouped[key]
        unit_id = f"PFU-{sequence:06d}"
        directions = sorted({as_text(row["direction"]) for row in rows})
        direction_state = directions[0] if len(directions) == 1 else "conflict"
        tiers = sorted({as_text(row["evidence_tier"]) for row in rows})
        family_states = sorted({as_text(row["family_state"]) for row in rows})
        unit = {
            "evidence_unit_id": unit_id,
            "family_id": key[0],
            "family_state": ";".join(family_states),
            "comparison_id": key[1],
            "target_taxon": key[2],
            "target_rank": config.target_rank,
            "direction_state": direction_state,
            "record_count": len(rows),
            "report_count": len({as_text(row["report_id"]) for row in rows}),
            "evidence_tiers": ";".join(tiers),
            "record_ids": ";".join(sorted(as_text(row["record_id"]) for row in rows)),
            "report_ids": ";".join(sorted({as_text(row["report_id"]) for row in rows})),
            "source_taxa": ";".join(sorted({as_text(row["source_taxon"]) for row in rows})),
            "lineage_components": ";".join(sorted({as_text(row.get("lineage_component")) for row in rows if as_text(row.get("lineage_component"))})),
            "policy_hash": config.policy_hash,
        }
        evidence_units.append(unit)
        for row in rows:
            record_map.append({
                "record_id": as_text(row["record_id"]),
                "evidence_unit_id": unit_id,
                "family_id": key[0],
                "comparison_id": key[1],
                "target_taxon": key[2],
                "direction": as_text(row["direction"]),
                "mapping_state": as_text(row["mapping_state"]),
                "aggregation_status": "collapsed_to_evidence_unit",
                "non_voting_reason": "",
                "policy_hash": config.policy_hash,
            })

    family_groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for unit in evidence_units:
        family_groups[(str(unit["family_id"]), str(unit["target_taxon"]))].append(unit)

    family_taxon_states: list[dict[str, object]] = []
    for family_id, taxon in sorted(family_groups):
        units = family_groups[(family_id, taxon)]
        directions = {str(unit["direction_state"]) for unit in units}
        state = next(iter(directions)) if len(directions) == 1 and "conflict" not in directions else "conflict"
        family_taxon_states.append({
            "family_id": family_id,
            "target_taxon": taxon,
            "target_rank": config.target_rank,
            "family_direction_state": state,
            "comparison_count": len({str(unit["comparison_id"]) for unit in units}),
            "evidence_unit_count": len(units),
            "evidence_unit_ids": ";".join(sorted(str(unit["evidence_unit_id"]) for unit in units)),
            "policy_hash": config.policy_hash,
        })

    states_by_taxon: dict[str, list[dict[str, object]]] = defaultdict(list)
    for state in family_taxon_states:
        states_by_taxon[str(state["target_taxon"])].append(state)

    taxon_summaries: list[dict[str, object]] = []
    for taxon in sorted(states_by_taxon):
        states = states_by_taxon[taxon]
        counts = {direction: sum(str(item["family_direction_state"]) == direction for item in states) for direction in config.direction_states}
        conflicts = sum(str(item["family_direction_state"]) == "conflict" for item in states)
        selected, classification, evaluations = evaluate_direction_counts(counts, conflicts, config)
        taxon_summaries.append({
            "target_taxon": taxon,
            "target_rank": config.target_rank,
            "family_count": len(states),
            "evaluable_direction_family_count": sum(counts.values()),
            "conflict_family_count": conflicts,
            "direction_counts_json": json.dumps(counts, sort_keys=True, separators=(",", ":")),
            "direction_evaluations_json": json.dumps(evaluations, sort_keys=True, separators=(",", ":")),
            "selected_direction": selected,
            "signal_classification": classification,
            "family_ids": ";".join(sorted(str(item["family_id"]) for item in states)),
            "minimum_support": config.minimum_support,
            "maximum_opposition": config.maximum_opposition,
            "conflict_policy": config.conflict_policy,
            "policy_hash": config.policy_hash,
        })

    return AggregationResult(
        record_map=sorted(record_map, key=lambda item: str(item["record_id"])),
        evidence_units=evidence_units,
        family_taxon_states=family_taxon_states,
        taxon_summaries=taxon_summaries,
        abstentions=abstentions,
        policy_hash=config.policy_hash,
    )
