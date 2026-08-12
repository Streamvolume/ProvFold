# Alcohol-associated hepatitis registry example

This directory contains 137 literature-derived records already mapped to the
generic 16-field ProvFold schema. The source-location column links each record
to a report identifier and a page, table or figure location; `report_key.csv`
provides the bibliographic identity of every report.

The primary display point requires support from 2 declared provenance groups
and no opposing group. It retains 2 genus directions. Reducing minimum support
to 1 retains 18 directions, and omitting either contributing group removes
both reference directions. These outputs quantify the selected set's threshold
and provenance-group dependence.

Run all operations with:

```bash
python -m provfold.cli reproduce --recipe reproducibility/ah_case/recipe.json --output-dir results/ah_case
```
