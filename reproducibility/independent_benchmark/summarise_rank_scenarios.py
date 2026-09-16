#!/usr/bin/env python3
"""Extract the two lower-rank scenario summaries from the complete benchmark run."""
from __future__ import annotations
import csv,gzip,json,os
from collections import defaultdict
from pathlib import Path
B=Path(os.environ.get('PROVFOLD_BENCHMARK_OUTPUT',Path(__file__).resolve().parent/'generated_run')).resolve()
def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def write(p,rs):
    with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),lineterminator='\n');w.writeheader();w.writerows(rs)
def main():
    names={'rank_replacement_with_lineage','rank_replacement_without_lineage'}
    records=defaultdict(list);pred=defaultdict(lambda:defaultdict(set))
    with gzip.open(B/'generated/records.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line)
            if r['scenario_id'] in names:records[(r['scenario_id'],r['replicate_index'])].append(r)
    with gzip.open(B/'generated/predictions.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line)
            if r['scenario_id'] in names and r['selected_direction']:pred[(r['scenario_id'],r['replicate_index'])][r['method']].add((r['taxon_key'],r['selected_direction']))
    classes=rows(B/'distinct_estimator_equivalence_classes.csv');output=[]
    for r in rows(B/'replicate_audit.csv'):
        if r['scenario_id'] not in names:continue
        key=(r['scenario_id'],int(r['replicate_index']));rs=records[key];diff=0
        for c in classes:
            sets=[pred[key][m] for m in c['configured_cells'].split(';')]
            diff+=any(x!=sets[0] for x in sets[1:])
        output.append({**{k:r[k] for k in ['scenario_id','scenario_role','replicate_index','seed','truth_rows','record_rows']},'rank_replacement_rows':sum(x['source_rank']=='species' for x in rs),'lineage_missing_rows':sum(not x['target_taxon'] for x in rs),'metric_rows':r['metric_rows'],'within_class_selected_set_differences':diff,'excluded':r['excluded']})
    write(B/'rank_replacement_replicate_audit.csv',output)
    write(B/'rank_replacement_scenario_performance.csv',[r for r in rows(B/'scenario_performance_summary.csv') if r['scenario_id'] in names])
    if len(output)!=200 or any(r['within_class_selected_set_differences'] for r in output):raise SystemExit('Rank-scenario validation failed')
    print('PASS: 200 rank-scenario replicates')
if __name__=='__main__':main()
