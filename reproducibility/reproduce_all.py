#!/usr/bin/env python3
"""Run the packaged reproducibility workflow without network access or data substitution."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json, os, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--public-input-dir',type=Path)
    p.add_argument('--skip-benchmark',action='store_true',help='Explicitly skip full synthetic generation and its validators')
    p.add_argument('--figures',action='store_true',help='Requires Matplotlib, NumPy and Pillow')
    a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    env={**os.environ,'PYTHONPATH':str(ROOT/'src'),'PROVFOLD_BENCHMARK_OUTPUT':str(out/'benchmark'),'PROVFOLD_PUBLIC_OUTPUT':str(out/'public_evaluation'),'PROVFOLD_ALPHA_OUTPUT':str(out/'alpha_diversity'),'PROVFOLD_FIGURE_OUTPUT':str(out/'figures')}
    if a.public_input_dir:env['PROVFOLD_PUBLIC_INPUT_DIR']=str(a.public_input_dir.resolve())
    checks=[]
    def save():
        (out/'run_status.json').write_text(json.dumps({'software_version':'0.1.2','status':'FAIL' if any(c.get('exit_code',0) for c in checks) else ('PARTIAL' if any(c['status']=='NOT_RUN' for c in checks) else 'PASS'),'checks':checks},indent=2)+'\n')
    def run(label,args):
        cmd=[sys.executable,*args];started=datetime.now(timezone.utc).isoformat()
        with (out/(label+'.log')).open('w') as log:r=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        checks.append({'check':label,'command':cmd,'started_at':started,'completed_at':datetime.now(timezone.utc).isoformat(),'exit_code':r.returncode,'status':'PASS' if not r.returncode else 'FAIL'});save()
        print(label,checks[-1]['status'],flush=True)
        if r.returncode:raise SystemExit(r.returncode)
    run('unit_tests',['-m','unittest','discover','-s','tests','-v'])
    run('release_checks',['reproducibility/run_release_checks.py'])
    for name,recipe in [('minimal','examples/minimal/recipe.json'),('ah','reproducibility/ah_case/recipe.json'),('hcc','reproducibility/hcc_repeated_report_case/recipe.json')]:
        run(name,['-m','provfold.cli','reproduce','--recipe',recipe,'--output-dir',str(out/name)])
    # Release checks execute alpha-diversity recalculation and independent validation.
    if not a.skip_benchmark:
        for name,script in [('benchmark','run_independent_benchmark.py'),('metrics','validate_metrics.py'),('replay','validate_replay.py'),('thresholds_and_equivalence','postprocess_benchmark.py'),('attribution','validate_attribution.py'),('rank_scenarios','summarise_rank_scenarios.py')]:
            run(name,['reproducibility/independent_benchmark/'+script])
    else:checks.append({'check':'benchmark_and_validators','status':'NOT_RUN','reason':'--skip-benchmark requested'})
    empirical=['reproducibility/postprocess_empirical.py','--output-dir',str(out/'empirical_postprocessing')]
    if a.public_input_dir:
        for name,script in [('public','run_public_evaluation.py'),('public_validation','validate_public_evaluation.py')]:run(name,['reproducibility/public_evaluation/'+script])
        empirical.append('--include-public')
    else:checks.append({'check':'public_resource_analysis','status':'NOT_RUN','reason':'No --public-input-dir supplied'})
    run('empirical_postprocessing',empirical)
    validation=['reproducibility/validate_expected_outputs.py','--report',str(out/'expected_validation.json'),'--empirical-output',str(out/'empirical_postprocessing')]
    if a.skip_benchmark:validation.append('--skip-benchmark')
    if not a.public_input_dir:validation.append('--skip-public')
    run('expected_outputs',validation)
    if a.figures:
        run('figures',['reproducibility/figures/build_figures.py']);run('prisma',['reproducibility/figures/build_prisma_figure.py'])
    else:checks.append({'check':'figure_rendering','status':'NOT_RUN','reason':'Optional --figures not requested'})
    save()
if __name__=='__main__':main()
