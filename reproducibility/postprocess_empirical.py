#!/usr/bin/env python3
"""Rebuild HCC support provenance and public-arm agreement summaries."""
from __future__ import annotations
import argparse, csv, json, math, os, sys
from pathlib import Path
from dataclasses import replace
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from provfold.config import AggregationConfig
from provfold.compare import compare_factorial_methods, factorial_method_key
HCC = Path(__file__).resolve().parent/'hcc_repeated_report_case'
PUBLIC = Path(os.environ.get('PROVFOLD_PUBLIC_OUTPUT', Path(__file__).resolve().parent/'public_evaluation/generated_run')).resolve()
def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))

def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        if fields:
            writer.writeheader()
            writer.writerows(rows)

def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def hcc_reanalysis() -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    records = read_csv(HCC / "input.csv")
    config = AggregationConfig.from_dict(json.loads((HCC / "aggregation_config.json").read_text(encoding="utf-8")))
    factorial = compare_factorial_methods(records, config)
    selected_separate = [
        row for row in factorial
        if row["method"] == factorial_method_key("family", "resolved_target", "separate_votes")
        and row["selected_direction"]
    ]
    provenance_rows: list[dict[str, object]] = []
    for result in selected_separate:
        taxon = str(result["taxon_key"])
        direction = str(result["selected_direction"])
        contributing = [
            row for row in records
            if row["target_taxon"] == taxon
            and row["direction"] == direction
            and row["comparison_id"] in config.primary_comparisons
            and row["evidence_tier"] in config.primary_evidence_tiers
        ]
        groups = sorted({row["family_id"] for row in contributing})
        units = sorted({f"{row['family_id']}|{row['comparison_id']}" for row in contributing})
        provenance_rows.append({
            "taxon": taxon,
            "direction": direction,
            "supporting_provenance_group_count": len(groups),
            "supporting_provenance_groups": ";".join(groups),
            "supporting_group_comparison_unit_count": len(units),
            "supporting_group_comparison_units": ";".join(units),
            "single_group_only": len(groups) == 1,
        })

    grid_rows: list[dict[str, object]] = []
    for unit in ("report", "family"):
        for policy in ("ignore", "separate_votes", "family_state"):
            method = factorial_method_key(unit, "resolved_target", policy)
            selected = sorted(
                f"{row['taxon_key']}|{row['selected_direction']}"
                for row in factorial
                if row["method"] == method and row["selected_direction"]
            )
            grid_rows.append({
                "voting_unit": "report" if unit == "report" else "provenance_group",
                "comparison_policy": policy,
                "selected_pair_count": len(selected),
                "selected_pairs": ";".join(selected),
            })

    comparison_unit_count = len({(row["family_id"], row["comparison_id"]) for row in records})
    normalised_support = math.ceil((2 / 3) * comparison_unit_count)
    normalised = compare_factorial_methods(records, replace(config, minimum_support=normalised_support))
    normalised_selected = sorted(
        f"{row['taxon_key']}|{row['selected_direction']}"
        for row in normalised
        if row["method"] == factorial_method_key("family", "resolved_target", "separate_votes")
        and row["selected_direction"]
    )
    summary = {
        "record_count": len(records),
        "report_count": len({row["report_id"] for row in records}),
        "provenance_group_count": len({row["family_id"] for row in records}),
        "group_comparison_unit_count": comparison_unit_count,
        "fixed_support_2_selected_count": len(selected_separate),
        "fixed_support_2_single_group_selected_count": sum(bool(row["single_group_only"]) for row in provenance_rows),
        "normalised_support_formula": "ceil((2/3) * total group-comparison units)",
        "normalised_minimum_support": normalised_support,
        "normalised_selected_count": len(normalised_selected),
        "normalised_selected_pairs": normalised_selected,
    }
    return provenance_rows, grid_rows, summary

def binomial_cdf(k: int, n: int, p: float) -> float:
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(0, k + 1))

def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n == 0:
        return math.nan, math.nan
    if k == 0:
        lower = 0.0
    else:
        lo, hi = 0.0, k / n
        for _ in range(80):
            mid = (lo + hi) / 2
            survival = 1 - binomial_cdf(k - 1, n, mid)
            if survival > alpha / 2:
                hi = mid
            else:
                lo = mid
        lower = (lo + hi) / 2
    if k == n:
        upper = 1.0
    else:
        lo, hi = k / n, 1.0
        for _ in range(80):
            mid = (lo + hi) / 2
            cdf = binomial_cdf(k, n, mid)
            if cdf > alpha / 2:
                lo = mid
            else:
                hi = mid
        upper = (lo + hi) / 2
    return lower, upper

def cross_arm_reanalysis() -> list[dict[str, object]]:
    balance = read_csv(PUBLIC / "cross_arm_direction_balance.csv")
    jaccard = read_csv(PUBLIC / "cross_arm_jaccard_surface.csv")
    labels = {row["disease_id"]: row["disease_label"] for row in balance}
    output = []
    for disease_id in sorted(labels):
        rows = [
            row for row in balance
            if row["disease_id"] == disease_id and row["sign_agreement_status"] in {"agree", "disagree"}
        ]
        if not rows:
            continue
        raw_positive = sum(float(row["rawdata_direction_balance"]) > 0 for row in rows)
        literature_positive = sum(float(row["literature_direction_balance"]) > 0 for row in rows)
        agreement = sum(row["sign_agreement_status"] == "agree" for row in rows)
        n = len(rows)
        p_raw = raw_positive / n
        p_literature = literature_positive / n
        expected = p_raw * p_literature + (1 - p_raw) * (1 - p_literature)
        observed = agreement / n
        lower, upper = clopper_pearson(agreement, n)
        kappa = None if expected == 1 else (observed - expected) / (1 - expected)
        selected = next(
            row for row in jaccard
            if row["disease_id"] == disease_id
            and row["minimum_support"] == "2"
            and row["maximum_opposition"] == "0"
        )
        output.append({
            "disease_id": disease_id,
            "disease_label": labels[disease_id],
            "comparable_direction_count": n,
            "observed_sign_agreement_count": agreement,
            "observed_sign_agreement_rate": round(observed, 8),
            "exact_95_ci_lower": round(lower, 8),
            "exact_95_ci_upper": round(upper, 8),
            "rawdata_positive_count": raw_positive,
            "literature_positive_count": literature_positive,
            "marginal_null_expected_agreement": round(expected, 8),
            "cohens_kappa": None if kappa is None else round(kappa, 8),
            "support_2_rawdata_selected_count": int(selected["rawdata_selected_pair_count"]),
            "support_2_literature_selected_count": int(selected["literature_selected_pair_count"]),
            "support_2_intersection_count": int(selected["intersection_count"]),
            "support_2_union_count": int(selected["union_count"]),
            "support_2_jaccard": float(selected["Jaccard"]) if selected["Jaccard"] else None,
        })
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path('results/empirical_postprocessing'))
    parser.add_argument('--include-public', action='store_true')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    provenance, grid, summary = hcc_reanalysis()
    write_csv(args.output_dir/'support_provenance.csv', provenance)
    write_csv(args.output_dir/'voting_unit_by_comparison_policy.csv', grid)
    write_json(args.output_dir/'hcc_summary.json', summary)
    if args.include_public:
        write_csv(PUBLIC/'cross_arm_agreement_with_null.csv', cross_arm_reanalysis())
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
