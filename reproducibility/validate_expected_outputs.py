#!/usr/bin/env python3
"""Compare current generated scientific outputs with the distribution's expected fixtures."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BENCHMARK=Path(os.environ.get('PROVFOLD_BENCHMARK_OUTPUT',ROOT/'reproducibility/independent_benchmark/generated_run')).resolve()
PUBLIC=Path(os.environ.get('PROVFOLD_PUBLIC_OUTPUT',ROOT/'reproducibility/public_evaluation/generated_run')).resolve()

def records(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f);return r.fieldnames,list(r)

def compare(observed,expected):
    if not observed.is_file():return {'status':'FAIL','reason':'Generated file is missing'}
    af,a=records(observed);bf,b=records(expected)
    ca=Counter(tuple(sorted(r.items())) for r in a);cb=Counter(tuple(sorted(r.items())) for r in b)
    return {'status':'PASS' if af==bf and ca==cb else 'FAIL','rows':len(a),'expected_rows':len(b),'different_generated_rows':sum((ca-cb).values()),'missing_expected_rows':sum((cb-ca).values()),'headers_equal':af==bf}

def decompressed_sha(path):
    h=hashlib.sha256()
    with gzip.open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--skip-benchmark',action='store_true');p.add_argument('--skip-public',action='store_true');p.add_argument('--report',type=Path,default=Path('results/expected_validation.json'));p.add_argument('--empirical-output',type=Path)
    args=p.parse_args();checks={};not_run=[]
    if not args.skip_benchmark:
        e=ROOT/'reproducibility/independent_benchmark/expected'
        for n in ['replicate_metrics.csv','replicate_audit.csv','scenario_performance_summary.csv','failure_region_summary.csv','support_threshold_calibration.csv','distinct_estimator_equivalence_classes.csv','factorial_attribution_summary.csv','rank_replacement_scenario_performance.csv','rank_replacement_replicate_audit.csv']:
            checks['benchmark/'+n]=compare(BENCHMARK/n,e/n)
        for r in records(e/'generated_stream_identities.csv')[1]:
            f=BENCHMARK/'generated'/r['file']
            checks['benchmark/streams/'+r['file']]={'status':'PASS' if f.is_file() and decompressed_sha(f)==r['decompressed_sha256'] else 'FAIL'}
    else:not_run.append('benchmark')
    if not args.skip_public:
        e=ROOT/'reproducibility/public_evaluation/expected'
        for f in sorted(e.glob('*.csv')):checks['public/'+f.name]=compare(PUBLIC/f.name,f)
        if not e.exists():checks['public/expected_directory']={'status':'FAIL','reason':'Missing expected fixtures'}
    else:not_run.append('public_resource_analysis')
    if args.empirical_output:
        e=ROOT/'reproducibility/hcc_repeated_report_case'
        for n in ['support_provenance.csv','voting_unit_by_comparison_policy.csv']:checks['hcc/'+n]=compare(args.empirical_output/n,e/n)
    else:not_run.append('HCC support-provenance post-processing')
    failed=[k for k,v in checks.items() if v['status']!='PASS']
    result={'executed_at':datetime.now(timezone.utc).isoformat(),'status':'FAIL' if failed else ('PARTIAL' if not_run else 'PASS'),'failed_checks':failed,'not_run':not_run,'checks':checks,'scope':'Scientific CSV outputs including policy/simulation hashes and decompressed stream bytes. Timestamps, runtime, peak memory and machine paths are not deterministic acceptance criteria.'}
    args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if failed:raise SystemExit(1)
if __name__=='__main__':main()
