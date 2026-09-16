#!/usr/bin/env python3
"""Recalculate the cohort-deduplicated AH alpha-diversity branch.

The script retains the DerSimonian--Laird point estimator for
compatibility with the recovered RevMan analysis and adds Hartung--Knapp
small-sample inference.  It also exports study-level random-effects weights,
leave-one-cohort-out results and report-replacement sensitivity results.

Only the registered full-precision SMD and standard-error values are used.
No source graph is digitised by this program.
"""

from __future__ import annotations

import csv
import os
import json
import math
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "input.csv"
OUTPUT = Path(os.environ.get("PROVFOLD_ALPHA_OUTPUT", HERE / "generated_run")).resolve()

OUTCOMES = ("Shannon", "Observed OTUs/ACE", "Chao1", "Simpson")
DISPLAY_OUTCOME = {
    "Observed OTUs/ACE": "ACE richness",
}
PRIMARY_ROLE = "primary_cohort_deduplicated"
REPLACEMENT_ROLE = "replacement_sensitivity"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    material = list(rows)
    if not material:
        raise ValueError(f"refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(material[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(material)


def _beta_fraction(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function."""
    maximum_iterations = 400
    epsilon = 3e-14
    tiny = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for iteration in range(1, maximum_iterations + 1):
        m2 = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < epsilon:
            return h
    raise RuntimeError("incomplete-beta continued fraction did not converge")


def regularised_incomplete_beta(x: float, a: float, b: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_term = (
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    term = math.exp(log_term)
    if x < (a + 1.0) / (a + b + 2.0):
        return term * _beta_fraction(a, b, x) / a
    return 1.0 - term * _beta_fraction(b, a, 1.0 - x) / b


def student_t_cdf(value: float, df: int) -> float:
    if df < 1:
        raise ValueError("degrees of freedom must be positive")
    if value == 0:
        return 0.5
    x = df / (df + value * value)
    tail_component = 0.5 * regularised_incomplete_beta(x, df / 2.0, 0.5)
    return 1.0 - tail_component if value > 0 else tail_component


def student_t_quantile(probability: float, df: int) -> float:
    if not 0.5 < probability < 1.0:
        raise ValueError("this implementation expects an upper-half probability")
    lower, upper = 0.0, 1.0
    while student_t_cdf(upper, df) < probability:
        upper *= 2.0
    for _ in range(160):
        midpoint = (lower + upper) / 2.0
        if student_t_cdf(midpoint, df) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def dl_hk(rows: list[dict[str, str]]) -> tuple[dict[str, float | int], list[float]]:
    if len(rows) < 2:
        raise ValueError("a pooled model requires at least two study estimates")
    effects = [float(row["smd"]) for row in rows]
    variances = [float(row["standard_error"]) ** 2 for row in rows]
    fixed_weights = [1.0 / variance for variance in variances]
    fixed_mean = sum(weight * effect for weight, effect in zip(fixed_weights, effects)) / sum(fixed_weights)
    q = sum(weight * (effect - fixed_mean) ** 2 for weight, effect in zip(fixed_weights, effects))
    df = len(rows) - 1
    c_value = sum(fixed_weights) - sum(weight * weight for weight in fixed_weights) / sum(fixed_weights)
    tau2 = max(0.0, (q - df) / c_value)
    random_weights = [1.0 / (variance + tau2) for variance in variances]
    weight_sum = sum(random_weights)
    pooled = sum(weight * effect for weight, effect in zip(random_weights, effects)) / weight_sum
    dl_se = math.sqrt(1.0 / weight_sum)
    dl_z = pooled / dl_se
    dl_p = math.erfc(abs(dl_z) / math.sqrt(2.0))
    hk_scale = sum(weight * (effect - pooled) ** 2 for weight, effect in zip(random_weights, effects)) / df
    hk_unadjusted_se = math.sqrt(hk_scale / weight_sum)
    hk_unadjusted_t = pooled / hk_unadjusted_se if hk_unadjusted_se else math.copysign(math.inf, pooled)
    hk_unadjusted_p = 2.0 * (1.0 - student_t_cdf(abs(hk_unadjusted_t), df))
    # The modified Knapp--Hartung safeguard prevents Hartung--Knapp variance
    # from shrinking below the conventional random-effects variance when
    # very small k and an unusually small Q produce q < 1.
    hk_se = max(hk_unadjusted_se, dl_se)
    hk_t = pooled / hk_se if hk_se else math.copysign(math.inf, pooled)
    hk_p = 2.0 * (1.0 - student_t_cdf(abs(hk_t), df))
    hk_critical = student_t_quantile(0.975, df)
    i2 = max(0.0, (q - df) / q * 100.0) if q > 0 else 0.0
    result: dict[str, float | int] = {
        "k": len(rows),
        "pooled_smd": pooled,
        "dl_standard_error": dl_se,
        "dl_ci_low": pooled - 1.959963984540054 * dl_se,
        "dl_ci_high": pooled + 1.959963984540054 * dl_se,
        "dl_z": dl_z,
        "dl_p_value": dl_p,
        "hk_scale": hk_scale,
        "hk_unadjusted_standard_error": hk_unadjusted_se,
        "hk_unadjusted_t": hk_unadjusted_t,
        "hk_unadjusted_ci_low": pooled - hk_critical * hk_unadjusted_se,
        "hk_unadjusted_ci_high": pooled + hk_critical * hk_unadjusted_se,
        "hk_unadjusted_p_value": hk_unadjusted_p,
        "hk_variance_safeguard_applied": hk_unadjusted_se < dl_se,
        "hk_standard_error": hk_se,
        "hk_df": df,
        "hk_t": hk_t,
        "hk_ci_low": pooled - hk_critical * hk_se,
        "hk_ci_high": pooled + hk_critical * hk_se,
        "hk_p_value": hk_p,
        "q": q,
        "q_df": df,
        "tau2_dl": tau2,
        "i2_pct": i2,
    }
    return result, random_weights


def model_row(model_id: str, outcome: str, role: str, rows: list[dict[str, str]]) -> dict[str, object]:
    result, _ = dl_hk(rows)
    return {
        "model_id": model_id,
        "outcome": DISPLAY_OUTCOME.get(outcome, outcome),
        "source_outcome": outcome,
        "analysis_role": role,
        "member_reports": ";".join(row["report_id"] for row in rows),
        "member_cohort_ids": ";".join(sorted({row["cohort_id"] for row in rows})),
        **result,
        "estimator": "inverse-variance DerSimonian-Laird random effects",
        "primary_inference": (
            "Hartung-Knapp t interval and two-sided t test with the modified Knapp-Hartung "
            "variance safeguard (SE=max[SE_HK, SE_DL])"
        ),
        "compatibility_inference": "normal interval and two-sided normal test",
    }


def main() -> None:
    rows = read_csv(SOURCE)
    primary_models: dict[str, list[dict[str, str]]] = {
        outcome: [row for row in rows if row["outcome"] == outcome and row["analysis_role"] == PRIMARY_ROLE]
        for outcome in OUTCOMES
    }
    replacement_rows = {
        row["outcome"]: row for row in rows if row["analysis_role"] == REPLACEMENT_ROLE
    }

    model_rows: list[dict[str, object]] = []
    study_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    replacement_results: list[dict[str, object]] = []

    for outcome_index, outcome in enumerate(OUTCOMES, start=1):
        members = primary_models[outcome]
        result, weights = dl_hk(members)
        model_id = f"AH-ALPHA-{outcome_index:02d}-PRIMARY"
        model_rows.append(model_row(model_id, outcome, PRIMARY_ROLE, members))
        display_outcome = DISPLAY_OUTCOME.get(outcome, outcome)
        total_weight = sum(weights)
        for member, weight in zip(members, weights):
            se = float(member["standard_error"])
            smd = float(member["smd"])
            study_rows.append({
                "model_id": model_id,
                "outcome": display_outcome,
                "source_outcome": outcome,
                "report_id": member["report_id"],
                "cohort_id": member["cohort_id"],
                "study_label": member["study_label"],
                "smd": smd,
                "standard_error": se,
                "study_ci_low": smd - 1.959963984540054 * se,
                "study_ci_high": smd + 1.959963984540054 * se,
                "random_effect_weight_pct": 100.0 * weight / total_weight,
                "reported_case_n": member["reported_case_n"],
                "reported_comparator_n": member["reported_comparator_n"],
                "source_location": member["source_location"],
                "numerical_provenance": (
                    "Group summaries were digitised from the cited source figure with WebPlotDigitizer 4.6 "
                    "for the registered review; this analysis uses the full-precision "
                    "SMD and standard error without a new digitisation pass."
                ),
            })
        if len(members) >= 3:
            for omitted in members:
                remaining = [row for row in members if row is not omitted]
                loo = model_row(f"{model_id}-OMIT-{omitted['report_id']}", outcome, "leave_one_cohort_out", remaining)
                loo_rows.append({
                    "parent_model_id": model_id,
                    "outcome": display_outcome,
                    "source_outcome": outcome,
                    "omitted_report_id": omitted["report_id"],
                    "omitted_study_label": omitted["study_label"],
                    **{key: value for key, value in loo.items() if key not in {"model_id", "outcome", "analysis_role"}},
                })

        replacement = replacement_rows.get(outcome)
        if replacement is not None:
            eom = next(row for row in members if row["report_id"] == "RPT0015")
            alternate_members = [replacement if row is eom else row for row in members]
            alternate_id = f"AH-ALPHA-{outcome_index:02d}-GANESAN-REPLACEMENT"
            alternate = model_row(alternate_id, outcome, "report_replacement", alternate_members)
            model_rows.append(alternate)
            replacement_results.append({
                "outcome": display_outcome,
                "source_outcome": outcome,
                "replacement_metric": replacement["outcome"].replace("Observed OTUs/ACE", "Observed OTUs"),
                "primary_model_id": model_id,
                "replacement_model_id": alternate_id,
                "replaced_report_id": "RPT0015",
                "replacement_report_id": replacement["report_id"],
                "primary_pooled_smd": result["pooled_smd"],
                "replacement_pooled_smd": alternate["pooled_smd"],
                "absolute_smd_difference": abs(float(alternate["pooled_smd"]) - float(result["pooled_smd"])),
                "primary_hk_ci_low": result["hk_ci_low"],
                "primary_hk_ci_high": result["hk_ci_high"],
                "replacement_hk_ci_low": alternate["hk_ci_low"],
                "replacement_hk_ci_high": alternate["hk_ci_high"],
                "primary_hk_p_value": result["hk_p_value"],
                "replacement_hk_p_value": alternate["hk_p_value"],
            })

    write_csv(OUTPUT / "alpha_diversity_study_effects_and_weights.csv", study_rows)
    write_csv(OUTPUT / "alpha_diversity_model_summary_dl_hk.csv", model_rows)
    write_csv(OUTPUT / "alpha_diversity_leave_one_cohort_out_dl_hk.csv", loo_rows)
    write_csv(OUTPUT / "alpha_diversity_report_replacement_dl_hk.csv", replacement_results)
    summary = {
        "source": SOURCE.name,
        "primary_models": len(OUTCOMES),
        "primary_study_rows": len(study_rows),
        "leave_one_cohort_out_models": len(loo_rows),
        "report_replacement_models": len(replacement_results),
        "method": (
            "DerSimonian-Laird tau-squared and point estimate; Hartung-Knapp t inference with "
            "the modified Knapp-Hartung variance safeguard SE=max(SE_HK, SE_DL)"
        ),
        "simpson_orientation": "higher values denote greater dominance and therefore lower diversity in the registered source metrics",
    }
    (OUTPUT / "alpha_diversity_recalculation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
