"""Schema-validated configuration models for ProvFold."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Any, Mapping


DEFAULT_RANK_ORDER = (
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
)


def _tuple_text(value: Any, name: str) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError(f"{name} must be a non-empty list")
    result = tuple(str(item).strip() for item in value)
    if any(not item for item in result) or len(result) != len(set(result)):
        raise ValueError(f"{name} must contain unique non-empty strings")
    return result


@dataclass(frozen=True)
class AggregationConfig:
    """Configuration for validation and evidence aggregation."""

    analysis_id: str
    primary_comparisons: tuple[str, ...]
    primary_evidence_tiers: tuple[str, ...]
    direction_states: tuple[str, ...]
    eligible_family_states: tuple[str, ...]
    eligible_mapping_states: tuple[str, ...]
    eligible_record_states: tuple[str, ...]
    target_rank: str
    minimum_support: int = 2
    maximum_opposition: int = 0
    conflict_policy: str = "abstain"
    comparator_scope: str = "primary"
    rank_order: tuple[str, ...] = field(default=DEFAULT_RANK_ORDER)
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if not self.analysis_id.strip():
            raise ValueError("analysis_id must not be empty")
        if len(self.direction_states) != 2:
            raise ValueError("direction_states must contain exactly two states")
        if self.minimum_support < 1:
            raise ValueError("minimum_support must be at least 1")
        if self.maximum_opposition < 0:
            raise ValueError("maximum_opposition must be non-negative")
        if self.conflict_policy not in {"abstain", "opposition"}:
            raise ValueError("conflict_policy must be abstain or opposition")
        if self.comparator_scope not in {"primary", "all_eligible"}:
            raise ValueError("comparator_scope must be primary or all_eligible")
        if self.target_rank not in self.rank_order:
            raise ValueError("target_rank is absent from rank_order")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "AggregationConfig":
        allowed = {
            "schema_version",
            "analysis_id",
            "primary_comparisons",
            "primary_evidence_tiers",
            "direction_states",
            "eligible_family_states",
            "eligible_mapping_states",
            "eligible_record_states",
            "target_rank",
            "minimum_support",
            "maximum_opposition",
            "conflict_policy",
            "comparator_scope",
            "rank_order",
        }
        unknown = sorted(set(data) - allowed)
        if unknown:
            raise ValueError(f"unsupported configuration fields: {', '.join(unknown)}")
        required = allowed - {
            "schema_version",
            "minimum_support",
            "maximum_opposition",
            "conflict_policy",
            "comparator_scope",
            "rank_order",
        }
        missing = sorted(required - set(data))
        if missing:
            raise ValueError(f"missing configuration fields: {', '.join(missing)}")
        return cls(
            schema_version=str(data.get("schema_version", "1.0")),
            analysis_id=str(data["analysis_id"]),
            primary_comparisons=_tuple_text(data["primary_comparisons"], "primary_comparisons"),
            primary_evidence_tiers=_tuple_text(data["primary_evidence_tiers"], "primary_evidence_tiers"),
            direction_states=_tuple_text(data["direction_states"], "direction_states"),
            eligible_family_states=_tuple_text(data["eligible_family_states"], "eligible_family_states"),
            eligible_mapping_states=_tuple_text(data["eligible_mapping_states"], "eligible_mapping_states"),
            eligible_record_states=_tuple_text(data["eligible_record_states"], "eligible_record_states"),
            target_rank=str(data["target_rank"]).strip(),
            minimum_support=int(data.get("minimum_support", 2)),
            maximum_opposition=int(data.get("maximum_opposition", 0)),
            conflict_policy=str(data.get("conflict_policy", "abstain")),
            comparator_scope=str(data.get("comparator_scope", "primary")),
            rank_order=_tuple_text(data.get("rank_order", DEFAULT_RANK_ORDER), "rank_order"),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "analysis_id": self.analysis_id,
            "primary_comparisons": list(self.primary_comparisons),
            "primary_evidence_tiers": list(self.primary_evidence_tiers),
            "direction_states": list(self.direction_states),
            "eligible_family_states": list(self.eligible_family_states),
            "eligible_mapping_states": list(self.eligible_mapping_states),
            "eligible_record_states": list(self.eligible_record_states),
            "target_rank": self.target_rank,
            "minimum_support": self.minimum_support,
            "maximum_opposition": self.maximum_opposition,
            "conflict_policy": self.conflict_policy,
            "comparator_scope": self.comparator_scope,
            "rank_order": list(self.rank_order),
        }

    @property
    def policy_hash(self) -> str:
        payload = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SimulationConfig:
    """Independent hierarchical generator settings."""

    simulation_id: str
    seed: int
    n_taxa: int = 40
    n_families: int = 8
    participants_case: int = 30
    participants_control: int = 30
    reports_per_family: int = 2
    primary_comparisons_per_family: int = 1
    primary_comparison_fraction: float = 1.0
    persistent_fraction: float = 0.25
    heterogeneous_fraction: float = 0.15
    effect_size: float = 0.8
    between_family_sd: float = 0.35
    measurement_sd: float = 1.0
    zero_inflation_probability: float = 0.0
    detection_threshold: float = 0.25
    selective_reporting_probability: float = 0.8
    sign_error_probability: float = 0.03
    direction_missing_probability: float = 0.0
    rank_inflation_probability: float = 0.3
    rank_replacement_probability: float = 0.0
    rank_lineage_missing_probability: float = 0.0
    ambiguous_taxonomy_probability: float = 0.0
    comparator_leakage_probability: float = 0.15
    nested_report_probability: float = 0.5
    nested_report_fraction: float = 0.7
    missing_family_probability: float = 0.0
    family_merge_probability: float = 0.0
    family_split_probability: float = 0.0
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if not self.simulation_id.strip():
            raise ValueError("simulation_id must not be empty")
        for name in (
            "n_taxa",
            "n_families",
            "participants_case",
            "participants_control",
            "reports_per_family",
            "primary_comparisons_per_family",
        ):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be at least 1")
        for name in (
            "persistent_fraction",
            "heterogeneous_fraction",
            "selective_reporting_probability",
            "sign_error_probability",
            "direction_missing_probability",
            "rank_inflation_probability",
            "rank_replacement_probability",
            "rank_lineage_missing_probability",
            "ambiguous_taxonomy_probability",
            "comparator_leakage_probability",
            "nested_report_probability",
            "nested_report_fraction",
            "missing_family_probability",
            "family_merge_probability",
            "family_split_probability",
            "zero_inflation_probability",
            "primary_comparison_fraction",
        ):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1]")
        if self.persistent_fraction + self.heterogeneous_fraction > 1.0:
            raise ValueError("persistent_fraction + heterogeneous_fraction must not exceed 1")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SimulationConfig":
        fields = {item.name for item in __import__("dataclasses").fields(cls)}
        unknown = sorted(set(data) - fields)
        if unknown:
            raise ValueError(f"unsupported simulation fields: {', '.join(unknown)}")
        if "simulation_id" not in data or "seed" not in data:
            raise ValueError("simulation_id and seed are required")
        return cls(**dict(data))

    def to_dict(self) -> dict[str, Any]:
        return {item.name: getattr(self, item.name) for item in __import__("dataclasses").fields(self)}
