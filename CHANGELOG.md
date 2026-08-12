# Changelog

## 0.1.0 — 2026-08-13

- Added the generic schema, configuration model and deterministic aggregation API.
- Added explicit family, comparison, taxonomy, conflict and abstention states.
- Added naive comparators, threshold surfaces and leave-one-family-out analysis.
- Added hierarchical simulation, benchmark and recipe-based reproduction interfaces.
- Added a minimal example, an AH provenance-registry reproduction fixture and automated tests.
- Filtered non-primary evidence tiers before evidence-unit construction so that a corroborative record cannot veto a primary record in the same unit.
- Added the complete voting-unit × taxon-key × comparison-policy factorial comparison and a multi-primary-comparison simulation scenario.
- Put benchmark methods on the same primary scope, added precision and false-discovery proportion, and fixed the false recurrent-signal denominator to the shared truth space while reporting out-of-truth selections separately.
- Added explicit lower-rank replacement scenarios with and without registered lineage and tests for their voting and abstention behaviour.
