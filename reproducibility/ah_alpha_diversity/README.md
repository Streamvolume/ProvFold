# AH alpha-diversity synthesis

This directory contains source-linked, full-precision standardised mean
differences and standard errors for the cohort-deduplicated alpha-diversity
models, together with executable analysis and independent-validation scripts.
Eom and Ganesan represent the same NCT04339725 cohort: Eom contributes to the
primary models and Ganesan replaces it in the registered sensitivity analysis.
The calculations use inverse-variance DerSimonian-Laird random effects and
Hartung-Knapp t inference with the modified Knapp-Hartung variance safeguard,
SE=max(SE_HK, SE_DL). Unadjusted Hartung-Knapp and normal intervals are retained
as explicitly labelled sensitivity summaries. Every row retains its report, cohort and
source-location identity. Run:

```bash
python reproducibility/ah_alpha_diversity/recalculate_alpha_diversity.py
python reproducibility/ah_alpha_diversity/validate_alpha_diversity_independently.py
```
