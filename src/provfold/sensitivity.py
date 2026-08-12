"""Threshold-surface and leave-one-family-out analyses."""

from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Mapping

from .config import AggregationConfig
from .core import aggregate, as_text, prepare_records


def _primary_families(records: Iterable[Mapping[str, object]], config: AggregationConfig) -> list[str]:
    eligible, _, _ = prepare_records(records, config)
    return sorted({
        as_text(row["family_id"])
        for row in eligible
        if as_text(row["comparison_id"]) in config.primary_comparisons
        and as_text(row["evidence_tier"]) in config.primary_evidence_tiers
    })


def threshold_surface(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
    minimum_support_values: Iterable[int] | None = None,
    maximum_opposition_values: Iterable[int] | None = None,
) -> list[dict[str, object]]:
    """Evaluate the complete declared support/opposition grid."""

    materialised = [dict(row) for row in records]
    family_count = len(_primary_families(materialised, config))
    support_values = list(minimum_support_values or range(1, family_count + 1))
    opposition_values = list(maximum_opposition_values or range(0, max(1, family_count)))
    if any(value < 1 or value > family_count for value in support_values):
        raise ValueError("minimum support values must lie between 1 and the number of primary families")
    if any(value < 0 or value >= max(1, family_count) for value in opposition_values):
        raise ValueError("maximum opposition values must lie between 0 and family_count - 1")

    rows: list[dict[str, object]] = []
    for support in sorted(set(support_values)):
        for opposition in sorted(set(opposition_values)):
            setting = replace(config, minimum_support=support, maximum_opposition=opposition)
            result = aggregate(materialised, setting)
            for summary in result.taxon_summaries:
                rows.append({
                    "minimum_support": support,
                    "maximum_opposition": opposition,
                    "target_taxon": summary["target_taxon"],
                    "selected_direction": summary["selected_direction"],
                    "signal_classification": summary["signal_classification"],
                    "family_count": summary["family_count"],
                    "conflict_family_count": summary["conflict_family_count"],
                    "direction_counts_json": summary["direction_counts_json"],
                    "policy_hash": setting.policy_hash,
                })
    return rows


def leave_one_family_out(
    records: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Omit each eligible primary family exactly once and recompute summaries."""

    materialised = [dict(row) for row in records]
    families = _primary_families(materialised, config)
    rows: list[dict[str, object]] = []
    for family_id in families:
        subset = [row for row in materialised if as_text(row.get("family_id")) != family_id]
        result = aggregate(subset, config)
        if not result.taxon_summaries:
            rows.append({
                "omitted_family_id": family_id,
                "target_taxon": "",
                "selected_direction": "",
                "signal_classification": "no_evaluable_taxon",
                "family_count": 0,
                "conflict_family_count": 0,
                "direction_counts_json": "{}",
                "policy_hash": config.policy_hash,
            })
        for summary in result.taxon_summaries:
            rows.append({
                "omitted_family_id": family_id,
                "target_taxon": summary["target_taxon"],
                "selected_direction": summary["selected_direction"],
                "signal_classification": summary["signal_classification"],
                "family_count": summary["family_count"],
                "conflict_family_count": summary["conflict_family_count"],
                "direction_counts_json": summary["direction_counts_json"],
                "policy_hash": config.policy_hash,
            })
    return rows
