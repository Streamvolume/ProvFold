# Field dictionary

| Field | Requirement | Meaning |
|---|---|---|
| `record_id` | required, unique | Immutable identifier for one extracted source row. |
| `report_id` | required | Publication or report identity; it is not assumed to equal a cohort. |
| `family_id` | conditional | Declared provenance-group identity. The field name is retained for file-format compatibility and is required when `family_state` is eligible. |
| `family_state` | required | Provenance state such as `cohort_family`, `project_family`, `publication_family` or `indeterminate`. Eligibility is configured per analysis. |
| `comparison_id` | required | Comparator stratum. Primary strata are configured rather than hard-coded. |
| `evidence_tier` | required | Primary or labelled corroborative evidence role. |
| `source_taxon` | required | Immutable source-reported taxon label. |
| `source_rank` | required | Immutable source-reported rank or rank-like label. Rank eligibility is evaluated from `harmonized_rank`; `source_rank` is retained for audit. |
| `harmonized_taxon` | conditional | Taxon after an externally frozen nomenclature mapping. The source label is never overwritten. |
| `harmonized_rank` | conditional | Rank of `harmonized_taxon`. |
| `target_taxon` | conditional | Explicit target-rank lineage used for conservative upward folding. It is never inferred by name similarity. |
| `mapping_state` | required | Mapping provenance; eligible values are configured. |
| `direction` | conditional | One of the two configured directional states. |
| `lineage_component` | optional | Audit identifier for correlated ancestor/descendant source rows. It is retained in outputs and does not define a vote. |
| `eligibility` | required | Extraction eligibility state; voting states are configured. |
| `source_location` | required | File, page, table, figure or database-row provenance. |

Additional source-specific fields are retained by the loader but are not used unless an adapter explicitly maps them into this contract.

## Output and analysis fields

The companion reproducibility archive's `DATA_DICTIONARY.md` describes output fields, units and missing-state conventions by analysis class; `DATA_CATALOG.csv` records actual headers, and `machine_readable/TABLE_MAP.csv` links Tables S1–S47 to their sources and generation entry points. In software summaries, a blank selected direction is not zero support: inspect the signal classification and directional counts. `family_id` retains its input identity throughout the record map, evidence-unit and group-state outputs.
