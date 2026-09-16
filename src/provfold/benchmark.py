"""Benchmark metrics against independently generated latent truth."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping

from .compare import (
    FACTORIAL_METHODS,
    METHOD_DEFINITIONS,
    compare_factorial_methods,
    compare_methods,
)
from .config import AggregationConfig


def _safe_divide(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def _defined_divide(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _score_predictions(
    comparisons: Iterable[Mapping[str, object]],
    methods: Iterable[str],
    truth: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Score predictions on a shared truth space and separate out-of-space selections."""

    truth_map = {str(row["taxon"]): str(row["latent_state"]) for row in truth}
    if not truth_map:
        raise ValueError("truth is empty")
    predictions: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)
    for row in comparisons:
        predictions[str(row["method"])][str(row["taxon_key"])] = (
            str(row["selected_direction"]),
            str(row["signal_classification"]),
        )

    output: list[dict[str, object]] = []
    policy_hash = config.policy_hash
    for method in sorted(methods):
        predicted = predictions.get(method, {})
        tp = within_truth_fp = tn = fn = 0
        heterogeneity_preserved = heterogeneity_total = abstained = 0
        direction_states = set(config.direction_states)
        for taxon, latent_state in truth_map.items():
            selected, classification = predicted.get(taxon, ("", "not_evaluable"))
            true_states = {latent_state} if latent_state in direction_states else set()
            predicted_states = {selected} if selected in direction_states else set()
            tp_i = len(true_states & predicted_states)
            fp_i = len(predicted_states - true_states)
            fn_i = len(true_states - predicted_states)
            tp += tp_i
            within_truth_fp += fp_i
            fn += fn_i
            tn += len(direction_states) - tp_i - fp_i - fn_i
            predicted_is_directional = bool(predicted_states)
            if latent_state == "heterogeneous":
                heterogeneity_total += 1
                heterogeneity_preserved += not predicted_is_directional
            abstained += classification in {"not_evaluable", "conflict_only"}
        out_of_truth_selected = sum(
            taxon not in truth_map and selected in direction_states
            for taxon, (selected, _) in predicted.items()
        )
        false_positive = within_truth_fp + out_of_truth_selected
        selected_count = tp + false_positive
        precision = _defined_divide(tp, selected_count)
        false_discovery_proportion = _defined_divide(false_positive, selected_count)
        recall = _safe_divide(tp, tp + fn)
        output.append({
            "method": method,
            "taxon_count": len(truth_map),
            "true_positive": tp,
            "false_positive": false_positive,
            "within_truth_false_positive": within_truth_fp,
            "out_of_truth_selected_pair_count": out_of_truth_selected,
            "true_negative": tn,
            "false_negative": fn,
            "selected_pair_count": selected_count,
            "precision_defined": precision is not None,
            "precision": None if precision is None else round(precision, 8),
            "false_discovery_proportion": None if false_discovery_proportion is None else round(false_discovery_proportion, 8),
            "recall": round(recall, 8),
            "specificity": round(_safe_divide(tn, tn + within_truth_fp), 8),
            "f1": round(_safe_divide(2 * tp, 2 * tp + false_positive + fn), 8),
            "false_positive_rate": round(_safe_divide(within_truth_fp, within_truth_fp + tn), 8),
            "false_recurrent_signal_rate": round(_safe_divide(within_truth_fp, within_truth_fp + tn), 8),
            "heterogeneity_preservation_rate": round(_safe_divide(heterogeneity_preserved, heterogeneity_total), 8),
            "abstention_rate": round(_safe_divide(abstained, len(truth_map)), 8),
            "policy_hash": policy_hash,
        })
    return output


def benchmark_predictions(
    records: Iterable[Mapping[str, object]],
    truth: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Score the five-rule comparison catalogue on a shared truth space."""

    return _score_predictions(compare_methods(records, config), METHOD_DEFINITIONS, truth, config)


def benchmark_factorial_predictions(
    records: Iterable[Mapping[str, object]],
    truth: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Score the scope-matched factorial ablation on a shared truth space."""

    comparisons = compare_factorial_methods(records, config)
    return benchmark_factorial_comparisons(comparisons, truth, config)


def benchmark_factorial_comparisons(
    comparisons: Iterable[Mapping[str, object]],
    truth: Iterable[Mapping[str, object]],
    config: AggregationConfig,
) -> list[dict[str, object]]:
    """Score already-computed factorial rows without repeating the ablation."""

    comparisons = list(comparisons)
    factor_by_method = {
        str(row["method"]): {
            "method_display_name": row["method_display_name"],
            "voting_unit": row["voting_unit"],
            "taxon_key_policy": row["taxon_key_policy"],
            "comparison_policy": row["comparison_policy"],
            "scope_policy": row["scope_policy"],
        }
        for row in comparisons
    }
    output = _score_predictions(comparisons, FACTORIAL_METHODS, truth, config)
    for row in output:
        row.update(factor_by_method[str(row["method"])])
    return output
