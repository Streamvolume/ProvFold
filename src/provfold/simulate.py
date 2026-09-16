"""Independent cohort-to-participant-to-report evidence generator."""

from __future__ import annotations

import hashlib
import json
import random

from .config import SimulationConfig
from .models import SimulationResult


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _draw_measurement(rng: random.Random, mean: float, sd: float, zero_probability: float) -> float:
    return 0.0 if rng.random() < zero_probability else rng.gauss(mean, sd)


def simulate_records(config: SimulationConfig) -> SimulationResult:
    """Generate latent truth before participant measurements and report rows."""

    rng = random.Random(config.seed)
    truth: list[dict[str, object]] = []
    states: dict[str, str] = {}
    taxa = [f"Taxon_{index:03d}" for index in range(1, config.n_taxa + 1)]
    persistent_count = int(round(config.n_taxa * config.persistent_fraction))
    heterogeneous_count = int(round(config.n_taxa * config.heterogeneous_fraction))
    shuffled = list(taxa)
    rng.shuffle(shuffled)
    persistent = set(shuffled[:persistent_count])
    heterogeneous = set(shuffled[persistent_count:persistent_count + heterogeneous_count])
    for taxon in taxa:
        if taxon in persistent:
            state = rng.choice(("increased", "decreased"))
        elif taxon in heterogeneous:
            state = "heterogeneous"
        else:
            state = "null"
        states[taxon] = state
        truth.append({
            "taxon": taxon,
            "latent_state": state,
            "truth_direction": state if state in {"increased", "decreased"} else "",
            "truth_is_directional": state in {"increased", "decreased"},
        })

    cohort_effects: list[dict[str, object]] = []
    participant_values: dict[tuple[str, str], tuple[list[float], list[float]]] = {}
    for family_index in range(1, config.n_families + 1):
        family_id = f"FAM-{family_index:03d}"
        for taxon in taxa:
            state = states[taxon]
            if state == "increased":
                centre = config.effect_size
            elif state == "decreased":
                centre = -config.effect_size
            elif state == "heterogeneous":
                centre = config.effect_size * rng.choice((-1.0, 1.0))
            else:
                centre = 0.0
            effect = rng.gauss(centre, config.between_family_sd)
            controls = [
                _draw_measurement(rng, 0.0, config.measurement_sd, config.zero_inflation_probability)
                for _ in range(config.participants_control)
            ]
            cases = [
                _draw_measurement(rng, effect, config.measurement_sd, config.zero_inflation_probability)
                for _ in range(config.participants_case)
            ]
            observed_difference = _mean(cases) - _mean(controls)
            participant_values[(family_id, taxon)] = (cases, controls)
            cohort_effects.append({
                "family_id": family_id,
                "taxon": taxon,
                "latent_state": state,
                "cohort_effect": round(effect, 8),
                "case_mean": round(_mean(cases), 8),
                "control_mean": round(_mean(controls), 8),
                "observed_difference": round(observed_difference, 8),
                "participants_case": config.participants_case,
                "participants_control": config.participants_control,
            })

    observed_family_ids = {f"FAM-{index:03d}": f"FAM-{index:03d}" for index in range(1, config.n_families + 1)}
    for family_index in range(2, config.n_families + 1):
        family_id = f"FAM-{family_index:03d}"
        if rng.random() < config.family_merge_probability:
            observed_family_ids[family_id] = observed_family_ids[f"FAM-{family_index - 1:03d}"]

    records: list[dict[str, object]] = []
    record_number = 0
    for family_index in range(1, config.n_families + 1):
        true_family_id = f"FAM-{family_index:03d}"
        for report_index in range(1, config.reports_per_family + 1):
            report_id = f"RPT-{family_index:03d}-{report_index:02d}"
            indeterminate = rng.random() < config.missing_family_probability
            family_id = "" if indeterminate else observed_family_ids[true_family_id]
            if family_id and report_index > 1 and rng.random() < config.family_split_probability:
                family_id = f"{family_id}-S{report_index:02d}"
            family_state = "indeterminate" if indeterminate else "cohort_family"
            for taxon in taxa:
                cases, controls = participant_values[(true_family_id, taxon)]
                if report_index > 1 and rng.random() < config.nested_report_probability:
                    case_n = max(2, int(round(len(cases) * config.nested_report_fraction)))
                    control_n = max(2, int(round(len(controls) * config.nested_report_fraction)))
                    report_cases = rng.sample(cases, min(case_n, len(cases)))
                    report_controls = rng.sample(controls, min(control_n, len(controls)))
                    report_scope = "nested_subset"
                else:
                    report_cases = cases
                    report_controls = controls
                    report_scope = "full_cohort"
                for comparison_index in range(1, config.primary_comparisons_per_family + 1):
                    if config.primary_comparisons_per_family == 1:
                        comparison_id = "case_vs_control"
                        comparison_cases = report_cases
                        comparison_controls = report_controls
                    else:
                        comparison_id = f"case_vs_control_panel_{comparison_index:02d}"
                        case_n = max(2, int(round(len(report_cases) * config.primary_comparison_fraction)))
                        control_n = max(2, int(round(len(report_controls) * config.primary_comparison_fraction)))
                        comparison_cases = rng.sample(report_cases, min(case_n, len(report_cases)))
                        comparison_controls = rng.sample(report_controls, min(control_n, len(report_controls)))
                    difference = _mean(comparison_cases) - _mean(comparison_controls)
                    detected = abs(difference) >= config.detection_threshold
                    if not detected or rng.random() > config.selective_reporting_probability:
                        continue
                    direction = "increased" if difference > 0 else "decreased"
                    if rng.random() < config.sign_error_probability:
                        direction = "decreased" if direction == "increased" else "increased"
                    if rng.random() < config.direction_missing_probability:
                        direction = ""
                    ambiguous = rng.random() < config.ambiguous_taxonomy_probability
                    # Do not advance the base random stream when the new
                    # replacement perturbation is disabled.
                    replacement = (
                        config.rank_replacement_probability > 0
                        and rng.random() < config.rank_replacement_probability
                    )
                    if replacement:
                        lineage_known = (
                            config.rank_lineage_missing_probability <= 0
                            or rng.random() >= config.rank_lineage_missing_probability
                        )
                        source_rows = [(
                            f"{taxon}_variant_F{family_index:03d}_R{report_index:02d}",
                            "species",
                            taxon if lineage_known else "",
                            "rank_folded" if lineage_known else "ambiguous",
                        )]
                    else:
                        source_rows = [(taxon, "genus", taxon, "exact")]
                    if not replacement and rng.random() < config.rank_inflation_probability:
                        source_rows.append((
                            f"{taxon}_variant_F{family_index:03d}_R{report_index:02d}",
                            "species",
                            taxon,
                            "rank_folded",
                        ))
                    for source_taxon, source_rank, target_taxon, mapping_state in source_rows:
                        if ambiguous:
                            target_taxon = ""
                            mapping_state = "ambiguous"
                        record_number += 1
                        records.append({
                            "record_id": f"SIM-{record_number:07d}",
                            "report_id": report_id,
                            "family_id": family_id,
                            "family_state": family_state,
                            "comparison_id": comparison_id,
                            "evidence_tier": "primary",
                            "source_taxon": source_taxon,
                            "source_rank": source_rank,
                            "harmonized_taxon": source_taxon,
                            "harmonized_rank": source_rank,
                            "target_taxon": target_taxon,
                            "mapping_state": mapping_state,
                            "direction": direction,
                            "lineage_component": taxon,
                            "eligibility": "eligible",
                            "source_location": (
                                f"simulation:{config.simulation_id}:{report_id}:{taxon}:"
                                f"{report_scope}:{comparison_id}"
                            ),
                        })
                        if rng.random() < config.comparator_leakage_probability:
                            record_number += 1
                            records.append({
                                **records[-1],
                                "record_id": f"SIM-{record_number:07d}",
                                "comparison_id": f"secondary_comparator_{comparison_index:02d}",
                                "evidence_tier": "secondary",
                                "source_location": (
                                    f"simulation:{config.simulation_id}:{report_id}:{taxon}:"
                                    f"{report_scope}:secondary:{comparison_index:02d}"
                                ),
                            })

    truth.sort(key=lambda item: str(item["taxon"]))
    cohort_effects.sort(key=lambda item: (str(item["family_id"]), str(item["taxon"])))
    records.sort(key=lambda item: str(item["record_id"]))
    payload = {
        "config": config.to_dict(),
        "truth": truth,
        "cohort_effects": cohort_effects,
        "records": records,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    return SimulationResult(truth=truth, cohort_effects=cohort_effects, records=records, simulation_hash=digest)
