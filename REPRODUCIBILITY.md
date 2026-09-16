# Reproducing ProvFold 0.1.2

Run commands from the source root with Python 3.11 or later. The core has no third-party runtime dependencies. Install the package with `python -m pip install .`; the public-resource adapter additionally requires `python -m pip install '.[reproduction]'`.

## Complete run

```bash
python reproducibility/reproduce_all.py --output-dir results/complete --public-input-dir /absolute/path/to/source_inputs
```

The command runs unit tests, release checks, the minimal recipe, AH and HCC recipes, alpha-diversity recalculation and its independent numerical validator, the complete synthetic benchmark, metric validation, deterministic replay samples, all-scenario threshold calibration and estimator-equivalence checks, and the three-arm public-resource analysis and its independent validator. It exits non-zero on failure. With no `--public-input-dir`, public-resource analysis is explicitly recorded as NOT_RUN rather than passed. Figure rendering is a separate optional step requiring Matplotlib, NumPy and Pillow.

`run_status.json` records commands, real execution timestamps and exit codes. Generated streams and analytical outputs are written to the chosen output directory. Benchmark generation comprises 18 scenarios, 1,720 replicates and 30,960 configured-cell evaluations. All 18 scenarios receive threshold post-processing; the 18 configurations represent 10 distinct estimators. Metric validation independently calculates arithmetic from stored predictions; replay uses the package generator and is not an independent generator implementation.

## Input boundary

The distribution contains all synthetic parameters, the seed schedule, source-linked AH/HCC extraction inputs, analytical recipes, validators, figure inputs and expected analytical outputs. The five public-resource source files are specified by filename, size and SHA-256 in `reproducibility/public_evaluation/source_input_identities.csv`. The adapter rejects any identity mismatch. Keep source exports in a separate local input directory, not in the public repository. See `SOURCE_INPUTS.md` for acquisition and licence boundaries. All public-resource arms remain separate.

## Expected-output contract

`python reproducibility/validate_expected_outputs.py` compares the generated benchmark and public-resource tables against the bundled expectations. Set `PROVFOLD_BENCHMARK_OUTPUT` and `PROVFOLD_PUBLIC_OUTPUT` to their output directories. Rows are compared by content, independent of row ordering. Deterministic policy/simulation identities and decompressed synthetic stream identities are checked. Runtime, peak memory, execution timestamps and machine paths are observational metadata, not deterministic scientific outputs. Gzip container bytes may depend on the Python/zlib implementation; canonical decompressed SHA-256 values are supplied.

## Individual entry points

- `python reproducibility/run_release_checks.py`
- `provfold reproduce --recipe examples/minimal/recipe.json --output-dir results/minimal`
- `provfold reproduce --recipe reproducibility/ah_case/recipe.json --output-dir results/ah_case`
- `provfold reproduce --recipe reproducibility/hcc_repeated_report_case/recipe.json --output-dir results/hcc_case`
- `python reproducibility/ah_alpha_diversity/recalculate_alpha_diversity.py`
- `python reproducibility/ah_alpha_diversity/validate_alpha_diversity_independently.py`
- `python reproducibility/independent_benchmark/run_independent_benchmark.py`
- `python reproducibility/independent_benchmark/validate_metrics.py`
- `python reproducibility/independent_benchmark/validate_replay.py`
- `python reproducibility/independent_benchmark/postprocess_benchmark.py`
- `python reproducibility/public_evaluation/run_public_evaluation.py`
- `python reproducibility/public_evaluation/validate_public_evaluation.py`
- `python reproducibility/figures/build_figures.py`
- `python reproducibility/figures/build_prisma_figure.py`

## Licences and integrity

Software is MIT licensed. Author-produced reproduction inputs and derived tables are CC BY 4.0. Third-party source exports retain their own rights; no blanket CC BY licence is applied to them. `SHA256SUMS.csv` covers distribution files except itself. The release checker validates that manifest before analytical checks.
