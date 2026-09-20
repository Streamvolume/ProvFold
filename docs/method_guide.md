# Method guide

ProvFold separates five operations that flat vote counting often conflates.

1. Validate source and mapping states without inferring missing provenance.
2. Resolve only explicit target-rank lineages; higher ranks never map downwards.
3. Collapse eligible rows within provenance group × comparison × target taxon.
4. Collapse multiple primary comparisons to one provenance-group taxon state.
5. Evaluate both directions against a declared minimum-support and maximum-opposition rule.

Here `family_id` names an analyst-declared provenance group, not the taxonomic family rank. A record abstains when its provenance group is indeterminate, its direction is unsupported, its mapping is ineligible, a lower-rank record lacks a frozen target lineage, or a higher-rank record would require downward inference. Otherwise eligible records outside the configured comparison or evidence tier remain out of scope. Opposite directions in the same unit or provenance group remain explicit conflicts. With the default `abstain` conflict policy, conflicts do not support either direction; with `opposition`, each conflict counts against both directional hypotheses.

A direction is retained only if support is at least `minimum_support`, opposition is at most `maximum_opposition`, and support is strictly greater than opposition. With exactly two configured directions, these conditions allow at most one direction to pass.

The code exposes comparison-level units and then one state per provenance group and target taxon. Under the current binary direction and unanimity rules, this result is identical to direct provenance-group collapse with comparisons ignored: both retain a direction exactly when every eligible row for that provenance group and target taxon has that direction. The comparison-level output remains useful for identifying the panels that created agreement or conflict, but it is not a separate estimator under this configuration.

Comparator counts use different evidence units and must not be presented as an attrition series. The factorial comparison uses controlled display names that state three factors in order: voting unit (`row`, `report` or `family`), taxon key (`source label` or `resolved target taxon`) and comparison handling (`ignored`, `separate votes` or `family-state collapse`). Stable machine keys are retained in the output crosswalk. The five-rule interface remains available for compatibility, but methodological attribution should use the complete factorial comparison on a common primary scope.

The `policy_hash` is the SHA-256 identity of the complete validated aggregation configuration in the current release. It is intended for exact within-release reproduction. When software versions expose different configuration fields, policy equivalence must be assessed by comparing the serialised fields rather than by assuming that equal analytical intent produces an equal hash.

`comparator_scope=primary` restricts five-rule comparators to configured primary comparisons and evidence tiers. Scope-matched method comparisons should use this setting. `comparator_scope=all_eligible` is a diagnostic option that exposes otherwise eligible strata to five-rule comparators; it must not be compared causally with a method restricted to primary scope.

The threshold surface reports every support value from 1 to the number of evaluable primary families and every opposition allowance from 0 to one less than that number. Leave-one-family-out analysis removes all records belonging to the omitted family, including secondary panels.

The simulation interface generates latent population states before participant measurements and reports. It is therefore suitable for method benchmarking, whereas a scenario whose expected output is defined by the aggregation rule is only a structural regression test.

Example conflict-policy configuration:

```json
{
  "minimum_support": 2,
  "maximum_opposition": 0,
  "conflict_policy": "opposition"
}
```

Only the policy field changes; the input ledger and family declarations remain fixed so that the output can be compared with the default `abstain` run.
