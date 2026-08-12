# Official source inputs for the public-resource reconstruction

The public-resource analyses used independently downloaded official sources.
These source exports are not redistributed in the ProvFold archive because an
affirmative database-export licence was not located. Retrieve the five files
listed in `source_input_identities.csv`, preserve the listed file names and
place them in one directory. Exact reconstruction requires the recorded
SHA-256 identities; a current file with a different identity is a new source
state and should be analysed as such.

Set that directory and a separate output directory before running the scripts:

```bash
export PROVFOLD_PUBLIC_INPUT_DIR=/path/to/source_inputs
export PROVFOLD_PUBLIC_OUTPUT=/path/to/public_evaluation_results
python reproducibility/public_evaluation/run_public_evaluation.py
python reproducibility/public_evaluation/validate_public_evaluation.py
```

## gutMDisorder exports labelled v3 by the service

Official resource page:
<https://bio-computing.hrbmu.edu.cn/gutMDisorder/resource.dhtml>

The official Resource page listed the disorder–health raw-data and literature
workbooks with filenames containing `v3` on 11 August 2026. The site's Release
& Version page still listed v2.0 as its latest named release, and the most
recent database paper describes v2.0; no separate v3 release note or paper was
located. The article therefore describes these as files labelled v3 by the
service, not as a documented v3 publication release.

Retrieve and verify the two registered files with:

```bash
python reproducibility/public_evaluation/retrieve_gutmdisorder_inputs.py
```

The script uses the exact official download endpoint exposed by the Resource
page and fails loudly on a SHA-256 mismatch. A current file with a different
identity is a new source state and must not be described as an exact replay.

## NCBI Taxonomy

Official E-utilities documentation:
<https://www.ncbi.nlm.nih.gov/books/NBK25499/>

The registered analysis queried numeric Taxonomy identifiers only and used an
explicit genus lineage. No approximate name matching was performed. The
redistributable 503-identifier query is included. To reconstruct the mapping
from NCBI EFetch, run:

```bash
python reproducibility/public_evaluation/fetch_ncbi_taxonomy.py
```

The script writes `ncbi_taxonomy_mapping.csv` into
`PROVFOLD_PUBLIC_INPUT_DIR`. A changed mapping hash records taxonomy drift and
must not be silently substituted for the recorded mapping.

## Disbiome

Official export page:
<https://disbiome.ugent.be/export>

The strict analysis required human faecal disease–healthy-control comparisons
with a qualitative direction. It used the official experiment export and a
complete publication-level accession lookup saved as
`publication_accession_numbers.json`. Experiment-family identity was not
recoverable from either source; those records therefore abstained from
recurrence analysis. If the current official service returns a different
snapshot, retain the new date and identity rather than describing the run as
an exact reconstruction.

## Interpretation

The gutMDisorder raw-data arm uses a BioProject-resolved project family. The
literature arm uses a publication family, not a cohort family. The Disbiome arm
is a determinability and abstention analysis. The three arms must not be
pooled.
