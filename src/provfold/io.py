"""Deterministic file input, output and manifest helpers."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping

from . import __version__
from .config import AggregationConfig, SimulationConfig


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: str | Path) -> list[dict[str, str]]:
    source = Path(path)
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def read_json(path: str | Path) -> object:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_aggregation_config(path: str | Path) -> AggregationConfig:
    payload = read_json(path)
    if not isinstance(payload, dict):
        raise ValueError("aggregation configuration must be a JSON object")
    return AggregationConfig.from_dict(payload)


def load_simulation_config(path: str | Path) -> SimulationConfig:
    payload = read_json(path)
    if not isinstance(payload, dict):
        raise ValueError("simulation configuration must be a JSON object")
    return SimulationConfig.from_dict(payload)


def write_json(path: str | Path, value: object) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def write_csv(path: str | Path, rows: Iterable[Mapping[str, object]], fieldnames: Iterable[str] | None = None) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    materialised = [dict(row) for row in rows]
    columns = list(fieldnames or (materialised[0].keys() if materialised else ()))
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="raise")
        if columns:
            writer.writeheader()
            for row in materialised:
                writer.writerow({column: row.get(column, "") for column in columns})


def write_manifest(
    output_dir: str | Path,
    input_files: Iterable[str | Path],
    config_hash: str,
    extra: Mapping[str, object] | None = None,
) -> dict[str, object]:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    input_entries = []
    for path in sorted((Path(item) for item in input_files), key=lambda item: str(item)):
        input_entries.append({"path": path.name, "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    output_entries = []
    for path in sorted(directory.rglob("*"), key=lambda item: str(item.relative_to(directory))):
        if path.is_file() and path.name != "run_manifest.json":
            output_entries.append({"path": str(path.relative_to(directory)), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "software": "provfold",
        "software_version": __version__,
        "policy_hash": config_hash,
        "inputs": input_entries,
        "outputs": output_entries,
    }
    if extra:
        manifest["run_metadata"] = dict(extra)
    write_json(directory / "run_manifest.json", manifest)
    return manifest


def write_aggregation_result(output_dir: str | Path, result: object) -> None:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    write_csv(directory / "record_map.csv", result.record_map)
    write_csv(directory / "evidence_units.csv", result.evidence_units)
    write_csv(directory / "family_taxon_states.csv", result.family_taxon_states)
    write_csv(directory / "taxon_summaries.csv", result.taxon_summaries)
    write_csv(directory / "abstentions.csv", result.abstentions)
    retained = [row for row in result.taxon_summaries if row["signal_classification"] == "retained"]
    lines = [
        "# ProvFold aggregation report",
        "",
        f"Policy hash: `{result.policy_hash}`",
        f"Evidence units: {len(result.evidence_units)}",
        f"Family–taxon states: {len(result.family_taxon_states)}",
        f"Taxa evaluated: {len(result.taxon_summaries)}",
        f"Records abstaining: {len(result.abstentions)}",
        f"Retained directional signals: {len(retained)}",
        "",
        "Counts use different evidence units in comparator outputs and must not be read as an attrition series.",
    ]
    (directory / "audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
