# Troubleshooting

Run `provfold validate --input input.csv --config config.json` before aggregation. The report gives an `_issue` code, record or family identifier, field and message. Preserve the original source value in an extra provenance column when an adapter normalises terminology.

| `_issue` code | Typical example | Corrective action |
|---|---|---|
| `missing_field` | The CSV has no `comparison_id` column. | Add every required schema column; do not silently substitute a similarly named source column. |
| `empty_identifier` | `record_id` is blank. | Assign a stable unique identifier derived from the extraction ledger. |
| `empty_required_value` | `source_taxon` or another required value is empty. | Restore the source value, or declare a non-voting state if it cannot be recovered. |
| `eligible_family_without_id` | `family_state=cohort` but `family_id` is blank. | Supply a defensible provenance-group identifier or change the state to `indeterminate`. |
| `unsupported_direction` | The source says `enriched`, while the configuration accepts `increased`/`decreased`. | Map the term in an adapter, preserve the original term and use a configured direction. |
| `unsupported_rank` | `harmonized_rank=OTU` is absent from `rank_order`. | Add a justified rank or retain an unresolved mapping; do not guess a rank. |
| `eligible_record_without_harmonized_taxon` | A mapped eligible row lacks `harmonized_taxon`. | Supply the frozen mapped taxon and rank, or change `mapping_state` so the row abstains. |
| `target_taxon_mismatch` | Genus-ranked `Bacteroides` has target taxon `Prevotella`. | Correct the frozen mapping; target-rank records must map to themselves. |
| `higher_rank_cannot_map_down` | A family-ranked record is supplied for a genus target. | Choose a compatible higher target rank or accept an abstention. |
| `lower_rank_without_target_lineage` | A species row lacks its frozen genus ancestor. | Supply `target_taxon` from a registered lineage source or accept an abstention. |
| `duplicate_identifier` | Two rows use `record_id=R017`. | Assign distinct stable identifiers; duplicate IDs are not automatically duplicate evidence. |
| `inconsistent_family_state` | One `family_id` is labelled both `cohort` and `publication`. | Resolve the registry inconsistency or split identifiers when the units are genuinely different. |
| `report_maps_to_multiple_families` | One `report_id` points to more than one `family_id`. | Correct the provenance registry. A report cannot create a synthetic compound group. |

If no signal is retained, inspect `family_taxon_states.csv`, the threshold surface and provenance-group omission output before changing a rule. A `not_retained` cell may reflect inadequate support, excess opposition or the strict `support > opposition` condition.

Different factorial or five-rule comparator counts are expected when rows, reports, families, comparison panels and taxon keys are not one-to-one. They are parallel sensitivity outputs, not sequential filtering stages or estimates of one biological quantity.
