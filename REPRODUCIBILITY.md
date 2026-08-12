# Reproducibility guide

ProvFold separates the generic software from the inputs used to evaluate it.
The compact archive contains all package source, tests, examples, synthetic
benchmark parameters and redistribution-safe derived results.

## Core and minimal example

```bash
python -m unittest discover -s tests -v
python -m provfold.cli reproduce --recipe examples/minimal/recipe.json --output-dir results/minimal
```

## Independent benchmark

The benchmark grid preserves the 18 scenarios, seed schedule, 1,720
replicates and all 18 configured cells used in the article. The
equivalence audit maps these configurations to 10 distinct estimators.
The runner writes
four deterministic synthetic streams and the complete metric table. Runtime
and peak-memory columns are descriptive machine measurements and are not
expected to be byte-identical across systems.

```bash
python reproducibility/independent_benchmark/run_independent_benchmark.py
python reproducibility/independent_benchmark/validate_metrics.py
python reproducibility/independent_benchmark/validate_replay.py
```

## AH provenance-registry example

The qualitative example begins from 137 source-linked records and provides a
complete threshold and provenance-group omission analysis. The parallel
`reproducibility/ah_alpha_diversity/` directory contains the source-linked
cohort-deduplicated alpha-diversity inputs and derived model, leave-one-out and
report-replacement outputs from the same cohort registry. It also includes the
recorded search strategies, PRISMA counts, study-quality/risk-of-bias tables,
executable DerSimonian–Laird/Hartung–Knapp analysis and an independent numerical
validator.

```bash
python -m provfold.cli reproduce \
  --recipe reproducibility/ah_case/recipe.json \
  --output-dir results/ah_case
python reproducibility/ah_alpha_diversity/recalculate_alpha_diversity.py
python reproducibility/ah_alpha_diversity/validate_alpha_diversity_independently.py
```

## HCC repeated-report stress test

The purposively selected corpus contains 69 eligible genus-direction records
from 4 reports and 3 BioProject-resolved provenance groups. Two reports reuse
PRJNA428932 in the same HCC–healthy comparison strata. The example tests
report, separate-panel and collapsed-group counting under a fixed policy. Its
support ledger identifies selections generated within one provenance group.

```bash
python -m provfold.cli reproduce \
  --recipe reproducibility/hcc_repeated_report_case/recipe.json \
  --output-dir results/hcc_case
```

## Public-resource evaluation

The three source arms were analysed separately. Official gutMDisorder, NCBI
Taxonomy and Disbiome snapshots are not redistributed because an affirmative
export licence was not located. The archive includes the executable adapter,
independent validator, registered execution contract, numeric Taxonomy query, derived
summaries, exact source identities and retrieval instructions. After placing
identity-matched official inputs in a local directory, follow
`reproducibility/public_evaluation/SOURCE_INPUTS.md` to reconstruct the full
evaluation.

## Figure reproduction

The four article figures can be rebuilt from the bundled registered input
tables. Install Matplotlib, NumPy and Pillow, then run:

```bash
python reproducibility/figures/build_figures.py
```

## Integrity

`SHA256SUMS.csv` covers every release file except itself. The release checks
verify this manifest before running analytical tests. They also reproduce the
AH and HCC examples and compare each example's four expected CSV
outputs row by row and field by field, including method descriptions and
policy hashes; summary signal-count checks are retained as secondary guards.

The code is licensed under MIT. The redistributed reproduction inputs and
derived data tables are licensed under CC BY 4.0; official database exports
listed in `public_evaluation/SOURCE_INPUTS.md` are not included.
