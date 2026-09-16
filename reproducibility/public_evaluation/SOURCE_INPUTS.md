# Public-resource inputs for ProvFold 0.1.2

The analysis uses five source files with exact identities listed in `source_input_identities.csv` and enforced by `execution_contract.json`. The authors confirm an input acquisition date of **13 August 2026** for all five files. This is an input acquisition date, not the date of a software execution or an AH literature search.

The original database exports are kept outside the public distribution because affirmative permission to redistribute complete exports has not been established. Access from an official service does not itself grant redistribution rights. Users must obtain the identified files from the providers or another authorised source. The adapter fails rather than silently accepting different input bytes.

## Required files

- `gutMDisorder_v3_Rawdata-based_Disorder_Health.xlsx`
- `gutMDisorder_v3_Literature-based_Disorder_Health.xlsx`
- `ncbi_taxonomy_mapping.csv`
- `experiment.json`
- `publication_accession_numbers.json`

## Providers

gutMDisorder resource page: https://bio-computing.hrbmu.edu.cn/gutMDisorder/resource.dhtml . The two workbooks are described by their service-provided filenames containing `v3`; that label is not used to infer a separate database-paper release. `retrieve_gutmdisorder_inputs.py` requests the official workbooks and checks their identities.

NCBI Taxonomy E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25499/ . The distribution includes the 503-identifier query and `fetch_ncbi_taxonomy.py`. Only numeric identifiers and explicit genus lineage are used; no approximate name matching is performed. Verify the mapping checksum before use.

Disbiome export: https://disbiome.ugent.be/export . The analysis uses the experiment export and a complete publication-level accession lookup. Experiment-family identity is not inferred from publication identity: unresolved experiment families abstain in the strict analysis.

## Execute

```bash
export PROVFOLD_PUBLIC_INPUT_DIR=/absolute/path/to/source_inputs
export PROVFOLD_PUBLIC_OUTPUT=/absolute/path/to/results/public_evaluation
python reproducibility/public_evaluation/run_public_evaluation.py
python reproducibility/public_evaluation/validate_public_evaluation.py
```

The raw-data arm uses BioProject-resolved project groups; the literature arm uses publication groups; the Disbiome strict arm measures determinability and abstention. These arms are not pooled.
