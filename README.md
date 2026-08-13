# ProvFold

ProvFold is a dependency-light Python package for auditing qualitative
taxon–direction synthesis when source rows, reports, provenance groups,
comparison strata and taxonomic labels represent different analytical units.

The software validates a provenance schema, retains voting and abstention
paths, compares 10 distinct estimators represented by 18 compatibility
configurations, and generates complete support–opposition, provenance-group
omission and configuration-equivalence analyses.

## Installation

ProvFold requires Python 3.11 or later and has no third-party runtime
dependency.

```bash
python -m pip install .
```

For a network-isolated installation, build or obtain the wheel first and then
install it without consulting a package index:

```bash
python -m pip wheel --no-build-isolation --no-deps . -w wheelhouse
python -m pip install --no-index --find-links wheelhouse provfold==0.1.1
```

## Minimal run

```bash
provfold validate --input examples/minimal/input.csv --config examples/minimal/config.json
provfold aggregate --input examples/minimal/input.csv --config examples/minimal/config.json --output-dir results/aggregate
provfold compare --input examples/minimal/input.csv --config examples/minimal/config.json --output-dir results/compare
provfold sensitivity --input examples/minimal/input.csv --config examples/minimal/config.json --output-dir results/sensitivity
provfold omit-family --input examples/minimal/input.csv --config examples/minimal/config.json --output-dir results/omit_family
provfold reproduce --recipe examples/minimal/recipe.json --output-dir results/reproduction
```

The same operations are available from Python:

```python
from provfold import AggregationConfig, aggregate
from provfold.io import read_csv, read_json

config = AggregationConfig.from_dict(read_json("config.json"))
result = aggregate(read_csv("input.csv"), config)
print(result.taxon_summaries)
```

Every analytical output includes a policy hash: the SHA-256 identity of the
complete validated aggregation configuration for this software release. It
identifies the exact within-release configuration; schema changes mean that
cross-version policy equivalence must be checked field by field rather than by
hash alone. A reproduction manifest separately records the software version
and SHA-256 identities of the inputs, configuration and outputs.

## Evidence contract

The machine field `family_id` stores an analyst-declared provenance group. A
configurable voting unit may combine this group with comparison and target
taxon. Records with an indeterminate group, unsupported direction, unresolved
taxonomy or ineligible source state abstain. Opposite directions inside a unit
remain a conflict. Under the current binary unanimity rule, two-stage
comparison collapse and direct provenance-group collapse are algebraically
equivalent final estimators. The two-stage path retains comparison-level audit
detail. Separate-panel voting is available as an explicit sensitivity and can
reveal when repeated panels inflate nominal support.

Taxonomic folding is conservative. A source record at the target rank can map
directly. A lower-rank record can fold upwards only when the registered input
supplies an explicit target taxon. Higher-rank records cannot be mapped
downwards, and unresolved or ambiguous mappings abstain. ProvFold does not
resolve names by string similarity.

See `docs/method_guide.md`, `docs/field_dictionary.md` and the JSON schemas for
details.

## Reproducibility materials

The `reproducibility/` directory contains:

- the independent hierarchical benchmark and expected outputs;
- the alcohol-associated hepatitis provenance registry, qualitative
  sensitivity recipe, and source-linked cohort-deduplicated alpha-diversity
  inputs and summaries;
- the purposively selected hepatocellular carcinoma provenance stress test,
  its report–group registry and complete analytical recipe;
- a clean public-resource adapter, independent validator and registered controls;
- a checksum-enforcing helper for the two gutMDisorder workbooks labelled v3
  by the official resource service;
- redistribution-safe summaries from the non-pooled public-resource
  evaluation; and
- retrieval locations and expected identities for official source files that
  are not redistributed.

Run the package-local checks with:

```bash
python reproducibility/run_release_checks.py
```

The unit-test suite can also be run directly:

```bash
python -m unittest discover -s tests -v
```

From an uninstalled source tree, run the same standard-library suite with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
```

The complete 1,720-replicate benchmark is intentionally a separate command
because it creates approximately 45 MB of deterministic synthetic streams:

```bash
python reproducibility/independent_benchmark/run_independent_benchmark.py
python reproducibility/independent_benchmark/validate_metrics.py
python reproducibility/independent_benchmark/validate_replay.py
```

The public-resource reconstruction requires separately retrieved official
inputs and `openpyxl`; see `reproducibility/public_evaluation/SOURCE_INPUTS.md`.

The AH quantitative branch is reproduced independently of the ProvFold
aggregation core:

```bash
python reproducibility/ah_alpha_diversity/recalculate_alpha_diversity.py
python reproducibility/ah_alpha_diversity/validate_alpha_diversity_independently.py
```

## Intended use

Use ProvFold after extraction to declare comparison scope, provenance groups,
taxon mappings, recurrence thresholds and omission criteria. The primary
deliverables are the record-to-result map, estimator crosswalk, threshold
surface and group-dependence audit. Literature screening, scientific
provenance assignment, quantitative effect-size modelling and taxonomic
curation remain upstream analytical tasks.

## Licence and citation

The software is licensed under the MIT licence. The accompanying
reproducibility datasets are released under CC BY 4.0. Citation metadata are
provided in `CITATION.cff`. Version 0.1.1 is archived at
https://doi.org/10.5281/zenodo.21915120, and the continuing source repository is
https://github.com/Streamvolume/ProvFold.
