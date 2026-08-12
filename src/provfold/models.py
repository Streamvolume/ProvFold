"""Result containers used by the public ProvFold API."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ValidationIssue:
    level: str
    code: str
    record_id: str
    field: str
    message: str


@dataclass(frozen=True)
class ValidationReport:
    valid: bool
    record_count: int
    issues: list[dict[str, str]]


@dataclass(frozen=True)
class AggregationResult:
    record_map: list[dict[str, object]]
    evidence_units: list[dict[str, object]]
    family_taxon_states: list[dict[str, object]]
    taxon_summaries: list[dict[str, object]]
    abstentions: list[dict[str, object]]
    policy_hash: str


@dataclass(frozen=True)
class SimulationResult:
    truth: list[dict[str, object]]
    cohort_effects: list[dict[str, object]]
    records: list[dict[str, object]]
    simulation_hash: str
