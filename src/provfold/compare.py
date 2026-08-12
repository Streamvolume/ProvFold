"""Naive and provenance-aware comparison methods on common inputs."""

from __future__ import annotations

from collections import defaultdict
import json
from typing import Iterable, Mapping

from .config import AggregationConfig
from .core import aggregate, as_text, evaluate_direction_counts, prepare_records


METHOD_DEFINITIONS = {
    "source_row": "each eligible source row is a vote after target-rank mapping",
    "unique_publication": "rows collapse by report or source-document key and target taxon",
    "family_without_comparison_taxonomy": "rows collapse by declared family and immutable source taxon",
    "source_era_taxon": "rows collapse by family, comparison and immutable source taxon",
    "provenance_rank_folded": "family, comparison and frozen target-rank lineage controls are applied",
}

# Stable reader-facing labels. The machine keys above remain unchanged for API
# compatibility and are intentionally more compact than these display names.
METHOD_DISPLAY_NAMES = {
    "source_row": "source-row counting",
    "unique_publication": "report/source-document counting",
    "family_without_comparison_taxonomy": "family-only source-label counting",
    "source_era_taxon": "family–comparison source-label counting",
    "provenance_rank_folded": "provenance-aware rank-folded aggregation",
}

FACTORIAL_VOTING_UNITS = ("row", "report", "family")
FACTORIAL_TAXON_KEYS = ("source_label", "resolved_target")
FACTORIAL_COMPARISON_POLICIES = ("ignore", "separate_votes", "family_state")


def factorial_method_key(voting_unit: str, taxon_key: str, comparison_policy: str) -> str:
    """Return the stable machine key for one scope-matched factorial ablation."""

    return f"vu_{voting_unit}__taxon_{taxon_key}__comparison_{comparison_policy}"


FACTORIAL_METHODS = tuple(
    factorial_method_key(voting_unit, taxon_key, comparison_policy)
    for voting_unit in FACTORIAL_VOTING_UNITS
    for taxon_key in FACTORIAL_TAXON_KEYS
    for comparison_policy in FACTORIAL_COMPARISON_POLICIES
)

PROVFOLD_FACTORIAL_METHOD = factorial_method_key("family", "resolved_target", "family_state")


def factorial_method_display_name(voting_unit: str, taxon_key: str, comparison_policy: str) -> str:
    """Return a compact controlled label assembled from the three declared factors."""

    unit = {"row": "row", "report": "report", "family": "family"}[voting_unit]
    taxon = {"source_label": "source label", "resolved_target": "resolved target taxon"}[taxon_key]
    comparison = {
        "ignore": "comparisons ignored",
        "separate_votes": "comparison-specific votes",
        "family_state": "family-state collapse",
    }[comparison_policy]
    return f"{unit}; {taxon}; {comparison}"


def _in_scope(row: Mapping[str, object], config: AggregationConfig) -> bool:
    if config.comparator_scope == "all_eligible":
        return True
    return (
        as_text(row["comparison_id"]) in config.primary_comparisons
        and as_text(row["evidence_tier"]) in config.primary_evidence_tiers
    )


def _collapse_votes(
    rows: list[dict[str, object]],
    key_fields: tuple[str, ...],
    taxon_field: str,
) -> list[dict[str, str]]:
    groups: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = tuple(as_text(row[field]) for field in key_fields) + (as_text(row[taxon_field]),)
        groups[key].append(row)
    votes: list[dict[str, str]] = []
    for key in sorted(groups):
        directions = sorted({as_text(row["direction"]) for row in groups[key]})
        votes.append({
            "vote_id": "|".join(key),
            "taxon": key[-1],
            "direction": directions[0] if len(directions) == 1 else "conflict",
        })
    return votes


def _summarize(method: str, votes: list[dict[str, str]], config: AggregationConfig) -> list[dict[str, object]]:
    by_taxon: dict[str, list[dict[str, str]]] = defaultdict(list)
    for vote in votes:
        by_taxon[vote["taxon"]].append(vote)
    rows: list[dict[str, object]] = []
    policy_hash = config.policy_hash
    for taxon in sorted(by_taxon):
        taxon_votes = by_taxon[taxon]
        counts = {direction: sum(vote["direction"] == direction for vote in taxon_votes) for direction in config.direction_states}
        conflicts = sum(vote["direction"] == "conflict" for vote in taxon_votes)
        selected, classification, evaluations = evaluate_direction_counts(counts, conflicts, config)
        rows.append({
            "method": method,
            "method_definition": METHOD_DEFINITIONS[method],
            "taxon_key": taxon,
            "vote_count": len(taxon_votes),
            "conflict_vote_count": conflicts,
            "direction_counts_json": json.dumps(counts, sort_keys=True, separators=(",", ":")),
            "direction_evaluations_json": json.dumps(evaluations, sort_keys=True, separators=(",", ":")),
            "selected_direction": selected,
            "signal_classification": classification,
            "minimum_support": config.minimum_support,
            "maximum_opposition": config.maximum_opposition,
            "policy_hash": policy_hash,
        })
    return rows


def _factorial_votes(
    rows: list[dict[str, object]],
    voting_unit: str,
    taxon_key: str,
    comparison_policy: str,
) -> list[dict[str, str]]:
    """Create votes for one voting-unit × taxon-key × comparison-policy cell."""

    if voting_unit not in FACTORIAL_VOTING_UNITS:
        raise ValueError(f"unsupported voting unit: {voting_unit}")
    if taxon_key not in FACTORIAL_TAXON_KEYS:
        raise ValueError(f"unsupported taxon key: {taxon_key}")
    if comparison_policy not in FACTORIAL_COMPARISON_POLICIES:
        raise ValueError(f"unsupported comparison policy: {comparison_policy}")

    unit_field = {"row": "record_id", "report": "report_id", "family": "family_id"}[voting_unit]
    taxon_field = {"source_label": "source_taxon", "resolved_target": "resolved_target_taxon"}[taxon_key]

    if comparison_policy == "ignore":
        return _collapse_votes(rows, (unit_field,), taxon_field)

    comparison_groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (
            as_text(row[unit_field]),
            as_text(row["comparison_id"]),
            as_text(row[taxon_field]),
        )
        comparison_groups[key].append(row)

    comparison_votes: list[dict[str, str]] = []
    for key in sorted(comparison_groups):
        group = comparison_groups[key]
        directions = sorted({as_text(row["direction"]) for row in group})
        family_ids = sorted({as_text(row["family_id"]) for row in group})
        comparison_votes.append({
            "vote_id": "|".join(key),
            "family_id": ";".join(family_ids),
            "taxon": key[-1],
            "direction": directions[0] if len(directions) == 1 else "conflict",
        })

    if comparison_policy == "separate_votes":
        return comparison_votes

    family_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for vote in comparison_votes:
        family_groups[(vote["family_id"], vote["taxon"])].append(vote)
    family_votes: list[dict[str, str]] = []
    for key in sorted(family_groups):
        directions = {vote["direction"] for vote in family_groups[key]}
        family_votes.append({
            "vote_id": "|".join(key),
            "taxon": key[-1],
            "direction": next(iter(directions)) if len(directions) == 1 and "conflict" not in directions else "conflict",
        })
    return family_votes


def compare_factorial_methods(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Evaluate a scope-matched 3 × 2 × 3 ablation on the same primary records."""

    eligible, _, _ = prepare_records(records, config)
    scoped = [
        row for row in eligible
        if as_text(row["comparison_id"]) in config.primary_comparisons
        and as_text(row["evidence_tier"]) in config.primary_evidence_tiers
    ]
    output: list[dict[str, object]] = []
    for voting_unit in FACTORIAL_VOTING_UNITS:
        for taxon_key in FACTORIAL_TAXON_KEYS:
            for comparison_policy in FACTORIAL_COMPARISON_POLICIES:
                method = factorial_method_key(voting_unit, taxon_key, comparison_policy)
                votes = _factorial_votes(scoped, voting_unit, taxon_key, comparison_policy)
                rows = _summarize_factorial(method, votes, config)
                for row in rows:
                    row.update({
                        "method_display_name": factorial_method_display_name(voting_unit, taxon_key, comparison_policy),
                        "voting_unit": voting_unit,
                        "taxon_key_policy": taxon_key,
                        "comparison_policy": comparison_policy,
                        "scope_policy": "primary_comparisons_and_primary_evidence_tiers",
                    })
                output.extend(rows)
    return sorted(output, key=lambda item: (str(item["method"]), str(item["taxon_key"])))


def _summarize_factorial(
    method: str,
    votes: list[dict[str, str]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Summarise a factorial cell without depending on the legacy method catalogue."""

    by_taxon: dict[str, list[dict[str, str]]] = defaultdict(list)
    for vote in votes:
        by_taxon[vote["taxon"]].append(vote)
    rows: list[dict[str, object]] = []
    policy_hash = config.policy_hash
    for taxon in sorted(by_taxon):
        taxon_votes = by_taxon[taxon]
        counts = {
            direction: sum(vote["direction"] == direction for vote in taxon_votes)
            for direction in config.direction_states
        }
        conflicts = sum(vote["direction"] == "conflict" for vote in taxon_votes)
        selected, classification, evaluations = evaluate_direction_counts(counts, conflicts, config)
        rows.append({
            "method": method,
            "taxon_key": taxon,
            "vote_count": len(taxon_votes),
            "conflict_vote_count": conflicts,
            "direction_counts_json": json.dumps(counts, sort_keys=True, separators=(",", ":")),
            "direction_evaluations_json": json.dumps(evaluations, sort_keys=True, separators=(",", ":")),
            "selected_direction": selected,
            "signal_classification": classification,
            "minimum_support": config.minimum_support,
            "maximum_opposition": config.maximum_opposition,
            "policy_hash": policy_hash,
        })
    return rows


def compare_methods(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Return all declared comparators without treating their counts as attrition."""

    eligible, _, _ = prepare_records(records, config)
    scoped = [row for row in eligible if _in_scope(row, config)]
    source_rows = [
        {"vote_id": as_text(row["record_id"]), "taxon": as_text(row["resolved_target_taxon"]), "direction": as_text(row["direction"])}
        for row in scoped
    ]
    publication = _collapse_votes(scoped, ("report_id",), "resolved_target_taxon")
    family_no_controls = _collapse_votes(scoped, ("family_id",), "source_taxon")
    source_era = _collapse_votes(scoped, ("family_id", "comparison_id"), "source_taxon")

    full_result = aggregate(records, config)
    full_rows = []
    for row in full_result.taxon_summaries:
        full_rows.append({
            "method": "provenance_rank_folded",
            "method_definition": METHOD_DEFINITIONS["provenance_rank_folded"],
            "taxon_key": row["target_taxon"],
            "vote_count": row["family_count"],
            "conflict_vote_count": row["conflict_family_count"],
            "direction_counts_json": row["direction_counts_json"],
            "direction_evaluations_json": row["direction_evaluations_json"],
            "selected_direction": row["selected_direction"],
            "signal_classification": row["signal_classification"],
            "minimum_support": row["minimum_support"],
            "maximum_opposition": row["maximum_opposition"],
            "policy_hash": row["policy_hash"],
        })

    output: list[dict[str, object]] = []
    for method, votes in (
        ("source_row", source_rows),
        ("unique_publication", publication),
        ("family_without_comparison_taxonomy", family_no_controls),
        ("source_era_taxon", source_era),
    ):
        output.extend(_summarize(method, votes, config))
    output.extend(full_rows)
    return sorted(output, key=lambda item: (str(item["method"]), str(item["taxon_key"])))
