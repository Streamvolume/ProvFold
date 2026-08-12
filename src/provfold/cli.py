"""Command-line interface calling the public ProvFold API."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .benchmark import benchmark_predictions
from .compare import compare_methods
from .core import aggregate, validate_records
from .io import (
    load_aggregation_config,
    load_simulation_config,
    read_csv,
    read_json,
    write_aggregation_result,
    write_csv,
    write_json,
    write_manifest,
)
from .sensitivity import leave_one_family_out, threshold_surface
from .simulate import simulate_records


def _common_input_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", required=True, help="Input CSV following schema version 1.0")
    parser.add_argument("--config", required=True, help="Aggregation configuration JSON")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="provfold", description="Provenance-aware qualitative evidence aggregation")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate", help="Validate input and configuration")
    _common_input_parser(validate)
    validate.add_argument("--output", help="Optional JSON validation report")

    aggregate_parser = subparsers.add_parser("aggregate", help="Build evidence units and taxon summaries")
    _common_input_parser(aggregate_parser)
    aggregate_parser.add_argument("--output-dir", required=True)

    compare = subparsers.add_parser("compare", help="Run declared comparator methods")
    _common_input_parser(compare)
    compare.add_argument("--output-dir", required=True)

    sensitivity = subparsers.add_parser("sensitivity", help="Run the complete support/opposition surface")
    _common_input_parser(sensitivity)
    sensitivity.add_argument("--output-dir", required=True)

    omission = subparsers.add_parser("omit-family", help="Omit every eligible family once")
    _common_input_parser(omission)
    omission.add_argument("--output-dir", required=True)

    simulate = subparsers.add_parser("simulate", help="Generate independent hierarchical truth and flat records")
    simulate.add_argument("--simulation-config", required=True)
    simulate.add_argument("--output-dir", required=True)

    benchmark = subparsers.add_parser("benchmark", help="Score methods against stored independent truth")
    _common_input_parser(benchmark)
    benchmark.add_argument("--truth", required=True)
    benchmark.add_argument("--output-dir", required=True)

    reproduce = subparsers.add_parser("reproduce", help="Run a frozen recipe through the public API")
    reproduce.add_argument("--recipe", required=True)
    reproduce.add_argument("--output-dir", required=True)
    return parser


def _load_common(args: argparse.Namespace):
    records = read_csv(args.input)
    config = load_aggregation_config(args.config)
    return records, config


def _run_recipe(recipe_path: Path, output_dir: Path) -> None:
    payload = read_json(recipe_path)
    if not isinstance(payload, dict):
        raise ValueError("recipe must be a JSON object")
    allowed = {"schema_version", "input", "config", "operations"}
    unknown = sorted(set(payload) - allowed)
    if unknown:
        raise ValueError(f"unsupported recipe fields: {', '.join(unknown)}")
    operations = payload.get("operations")
    if not isinstance(operations, list) or not operations:
        raise ValueError("recipe operations must be a non-empty list")
    permitted = {"validate", "aggregate", "compare", "sensitivity", "omit-family"}
    unsupported = sorted(set(str(item) for item in operations) - permitted)
    if unsupported:
        raise ValueError(f"unsupported recipe operations: {', '.join(unsupported)}")
    base = recipe_path.parent
    input_path = (base / str(payload["input"])).resolve()
    config_path = (base / str(payload["config"])).resolve()
    records = read_csv(input_path)
    config = load_aggregation_config(config_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    if "validate" in operations:
        report = validate_records(records, config)
        write_json(output_dir / "validation.json", report.__dict__)
        if not report.valid:
            raise ValueError("recipe validation failed")
    if "aggregate" in operations:
        write_aggregation_result(output_dir / "aggregate", aggregate(records, config))
    if "compare" in operations:
        write_csv(output_dir / "compare" / "method_comparison.csv", compare_methods(records, config))
    if "sensitivity" in operations:
        write_csv(output_dir / "sensitivity" / "threshold_surface.csv", threshold_surface(records, config))
    if "omit-family" in operations:
        write_csv(output_dir / "omit_family" / "leave_one_family_out.csv", leave_one_family_out(records, config))
    write_manifest(output_dir, [input_path, config_path, recipe_path], config.policy_hash, {"operations": operations})


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "validate":
            records, config = _load_common(args)
            report = validate_records(records, config)
            payload = report.__dict__
            if args.output:
                write_json(args.output, payload)
            else:
                print(json.dumps(payload, indent=2, sort_keys=True))
            return 0 if report.valid else 2

        if args.command == "simulate":
            simulation_config = load_simulation_config(args.simulation_config)
            result = simulate_records(simulation_config)
            output = Path(args.output_dir)
            write_csv(output / "truth.csv", result.truth)
            write_csv(output / "cohort_effects.csv", result.cohort_effects)
            write_csv(output / "records.csv", result.records)
            write_json(output / "simulation_summary.json", {"simulation_hash": result.simulation_hash, "config": simulation_config.to_dict()})
            write_manifest(output, [args.simulation_config], result.simulation_hash, {"operation": "simulate"})
            return 0

        if args.command == "reproduce":
            _run_recipe(Path(args.recipe).resolve(), Path(args.output_dir))
            return 0

        records, config = _load_common(args)
        output = Path(args.output_dir) if hasattr(args, "output_dir") else None
        if args.command == "aggregate":
            write_aggregation_result(output, aggregate(records, config))
        elif args.command == "compare":
            write_csv(output / "method_comparison.csv", compare_methods(records, config))
        elif args.command == "sensitivity":
            write_csv(output / "threshold_surface.csv", threshold_surface(records, config))
        elif args.command == "omit-family":
            write_csv(output / "leave_one_family_out.csv", leave_one_family_out(records, config))
        elif args.command == "benchmark":
            truth = read_csv(args.truth)
            write_csv(output / "benchmark_metrics.csv", benchmark_predictions(records, truth, config))
        else:
            raise ValueError(f"unsupported command: {args.command}")
        input_files = [args.input, args.config]
        if args.command == "benchmark":
            input_files.append(args.truth)
        write_manifest(output, input_files, config.policy_hash, {"operation": args.command})
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(f"provfold: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
