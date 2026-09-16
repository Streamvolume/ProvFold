#!/usr/bin/env python3
"""Independently validate the alpha-diversity DL and Hartung--Knapp results.

This validator deliberately uses numerical integration of the Student t
density rather than the incomplete-beta implementation used by the analysis
script.
"""

from __future__ import annotations

import csv
import os
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "input.csv"
OUTPUT = Path(os.environ.get("PROVFOLD_ALPHA_OUTPUT", HERE / "generated_run")).resolve()
RESULTS = OUTPUT / "alpha_diversity_model_summary_dl_hk.csv"
REPORT = OUTPUT / "alpha_diversity_independent_validation.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def t_density(x: float, df: int) -> float:
    coefficient = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    return coefficient * (1 + x * x / df) ** (-(df + 1) / 2)


def simpson_integral(function, lower: float, upper: float, intervals: int = 200000) -> float:
    if intervals % 2:
        intervals += 1
    width = (upper - lower) / intervals
    total = function(lower) + function(upper)
    total += 4 * sum(function(lower + index * width) for index in range(1, intervals, 2))
    total += 2 * sum(function(lower + index * width) for index in range(2, intervals, 2))
    return total * width / 3


def t_cdf(value: float, df: int) -> float:
    if value == 0:
        return 0.5
    area = simpson_integral(lambda x: t_density(x, df), 0.0, abs(value), 120000)
    return 0.5 + area if value > 0 else 0.5 - area


def t_quantile(probability: float, df: int) -> float:
    lower, upper = 0.0, 32.0
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if t_cdf(midpoint, df) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


def independent_model(rows: list[dict[str, str]]) -> dict[str, float | int]:
    effects = [float(row["smd"]) for row in rows]
    variances = [float(row["standard_error"]) ** 2 for row in rows]
    fixed = [1 / variance for variance in variances]
    fixed_mean = sum(w * y for w, y in zip(fixed, effects)) / sum(fixed)
    q = sum(w * (y - fixed_mean) ** 2 for w, y in zip(fixed, effects))
    df = len(rows) - 1
    c = sum(fixed) - sum(w * w for w in fixed) / sum(fixed)
    tau2 = max(0.0, (q - df) / c)
    weights = [1 / (variance + tau2) for variance in variances]
    pooled = sum(w * y for w, y in zip(weights, effects)) / sum(weights)
    dl_se = math.sqrt(1 / sum(weights))
    hk_scale = sum(w * (y - pooled) ** 2 for w, y in zip(weights, effects)) / df
    hk_unadjusted_se = math.sqrt(hk_scale / sum(weights))
    hk_se = max(hk_unadjusted_se, dl_se)
    critical = t_quantile(0.975, df)
    t_value = pooled / hk_se
    return {
        "k": len(rows), "pooled_smd": pooled, "tau2_dl": tau2,
        "dl_ci_low": pooled - 1.959963984540054 * dl_se,
        "dl_ci_high": pooled + 1.959963984540054 * dl_se,
        "hk_unadjusted_standard_error": hk_unadjusted_se,
        "hk_unadjusted_ci_low": pooled - critical * hk_unadjusted_se,
        "hk_unadjusted_ci_high": pooled + critical * hk_unadjusted_se,
        "hk_standard_error": hk_se,
        "hk_ci_low": pooled - critical * hk_se,
        "hk_ci_high": pooled + critical * hk_se,
        "hk_p_value": 2 * (1 - t_cdf(abs(t_value), df)),
    }


def main() -> None:
    if not RESULTS.is_file():
        raise SystemExit(
            "Alpha-diversity results are missing. Run recalculate_alpha_diversity.py first, "
            "then rerun this validator."
        )
    inputs = read_csv(SOURCE)
    observed = {
        row.get("source_outcome", row["outcome"]): row for row in read_csv(RESULTS)
        if row["analysis_role"] == "primary_cohort_deduplicated"
    }
    checks = []
    fields = (
        "pooled_smd", "tau2_dl", "dl_ci_low", "dl_ci_high",
        "hk_unadjusted_standard_error", "hk_unadjusted_ci_low", "hk_unadjusted_ci_high",
        "hk_standard_error", "hk_ci_low", "hk_ci_high", "hk_p_value",
    )
    for outcome, saved in observed.items():
        rows = [row for row in inputs if row["outcome"] == outcome and row["analysis_role"] == "primary_cohort_deduplicated"]
        recomputed = independent_model(rows)
        differences = {field: abs(float(saved[field]) - float(recomputed[field])) for field in fields}
        checks.append({
            "outcome": outcome,
            "k": len(rows),
            "maximum_absolute_difference": max(differences.values()),
            "status": "PASS" if max(differences.values()) < 5e-7 else "FAIL",
            "differences": differences,
        })
    payload = {
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "implementation_independence": "Student t probabilities and critical values recomputed by numerical integration",
        "checks": checks,
    }
    REPORT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
