#!/usr/bin/env python3
"""Fetch and parse the frozen public-resource NCBI Taxonomy query without aggregation."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


OUTPUT = Path(os.environ.get("PROVFOLD_PUBLIC_INPUT_DIR", Path(__file__).resolve().parent / "source_inputs")).resolve()
QUERY = Path(__file__).resolve().parent / "taxonomy_query_taxids.txt"
CONTRACT = Path(__file__).resolve().parent / "execution_contract.json"
ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 100


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def text(node: ET.Element, path: str) -> str:
    value = node.findtext(path)
    return "" if value is None else value.strip()


def main() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    taxids = [line.strip() for line in QUERY.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(taxids) != 503 or sha256(QUERY) != contract["source_identities"]["taxonomy_query_sha256"]:
        raise SystemExit("frozen taxonomy query identity mismatch")
    if len(taxids) != len(set(taxids)) or any(not item.isdigit() for item in taxids):
        raise SystemExit("taxonomy query must contain unique numeric identifiers")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    response_files: list[dict[str, object]] = []
    nodes_by_query: dict[str, dict[str, str]] = {}
    for start in range(0, len(taxids), BATCH_SIZE):
        batch = taxids[start : start + BATCH_SIZE]
        request = urllib.request.Request(
            ENDPOINT,
            data=urllib.parse.urlencode({
                "db": "taxonomy",
                "id": ",".join(batch),
                "retmode": "xml",
                "tool": "provfold",
            }).encode("ascii"),
            headers={"User-Agent": "ProvFold/0.1 taxonomy-metadata-audit"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
        batch_number = start // BATCH_SIZE + 1
        path = OUTPUT / f"taxonomy_efetch_batch_{batch_number:02d}.xml"
        path.write_bytes(payload)
        root = ET.fromstring(payload)
        for node in root.findall("./Taxon"):
            current_taxid = text(node, "TaxId")
            rank = text(node, "Rank").lower()
            scientific_name = text(node, "ScientificName")
            genus_taxid = current_taxid if rank == "genus" else ""
            genus_name = scientific_name if rank == "genus" else ""
            for ancestor in node.findall("./LineageEx/Taxon"):
                if text(ancestor, "Rank").lower() == "genus":
                    genus_taxid = text(ancestor, "TaxId")
                    genus_name = text(ancestor, "ScientificName")
            aliases = {current_taxid}
            aliases.update(value.text.strip() for value in node.findall("./AkaTaxIds/TaxId") if value.text)
            row = {
                "current_taxid": current_taxid,
                "scientific_name": scientific_name,
                "returned_rank": rank,
                "genus_taxid": genus_taxid,
                "genus_name": genus_name,
            }
            for alias in aliases:
                nodes_by_query[alias] = row
        response_files.append({
            "batch": batch_number,
            "query_count": len(batch),
            "filename": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
        if start + BATCH_SIZE < len(taxids):
            time.sleep(0.4)

    rows = []
    for query_taxid in taxids:
        resolved = nodes_by_query.get(query_taxid)
        if resolved is None:
            rows.append({
                "query_taxid": query_taxid,
                "resolution_status": "not_returned_abstain",
                "current_taxid": "",
                "scientific_name": "",
                "returned_rank": "",
                "genus_taxid": "",
                "genus_name": "",
            })
        else:
            rows.append({
                "query_taxid": query_taxid,
                "resolution_status": "resolved_current" if query_taxid == resolved["current_taxid"] else "resolved_aka_taxid",
                **resolved,
            })
    mapping = OUTPUT / "ncbi_taxonomy_mapping.csv"
    with mapping.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "schema_version": "1.0",
        "retrieval_date": datetime.now(timezone.utc).date().isoformat(),
        "input_acquisition_date": "2026-08-13",
        "input_acquisition_date_basis": "author-confirmed",
        "official_source": "NCBI Taxonomy EFetch",
        "endpoint": ENDPOINT,
        "query_count": len(taxids),
        "query_sha256": sha256(QUERY),
        "batch_size": BATCH_SIZE,
        "batches": response_files,
        "mapping_file": mapping.name,
        "mapping_size_bytes": mapping.stat().st_size,
        "mapping_sha256": sha256(mapping),
        "resolved_current_count": sum(row["resolution_status"] == "resolved_current" for row in rows),
        "resolved_alias_count": sum(row["resolution_status"] == "resolved_aka_taxid" for row in rows),
        "unreturned_count": sum(row["resolution_status"] == "not_returned_abstain" for row in rows),
        "genus_lineage_available_count": sum(bool(row["genus_taxid"]) for row in rows),
        "recurrent_taxa_computed": False,
    }
    manifest_path = OUTPUT / "ncbi_taxonomy_fetch_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
