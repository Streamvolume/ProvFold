# HCC repeated-report stress test

This purposively selected corpus contains 69 eligible genus-direction records
from 4 reports and 3 BioProject-resolved provenance groups. Two reports reuse
PRJNA428932 and contribute concordant decreased-Megamonas records in the same
non-HBV-HCC-versus-healthy comparison stratum.

At minimum support 2 and maximum opposition 0, report voting with comparisons
ignored retains 3 directions, separate provenance-group-by-comparison votes
retain 16, and collapsed provenance-group voting retains 2. Fourteen of the 16
separate-panel selections receive support from only one provenance group. The
ledger therefore identifies within-group panel pseudoreplication directly.

Run all operations with:

```bash
python -m provfold.cli reproduce --recipe reproducibility/hcc_repeated_report_case/recipe.json --output-dir results/hcc_case
```
