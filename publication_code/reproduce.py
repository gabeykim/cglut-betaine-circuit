"""Re-execute the manuscript stages in a NEW directory using cached inputs.
Usage: python publication_code/reproduce.py /absolute/path/to/new-output-directory
Install reproduction/requirements_reproduction.txt first; hmmscan must be on PATH.
Does not refetch inputs or repeat DeepTMHMM neural inference.
"""
from pathlib import Path
import csv, hashlib, json, os, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parents[1]
dest=Path(sys.argv[1]).resolve()
if dest.exists():raise SystemExit('Choose a new output directory; existing paths are not overwritten.')
if dest==ROOT or ROOT in dest.parents:raise SystemExit('Use an output directory outside the release package.')
if shutil.which('hmmscan') is None:raise SystemExit('hmmscan is not on PATH. Activate the scientific environment first.')
dest.mkdir(parents=True)
repo=dest/'project';shutil.copytree(ROOT/'project',repo)
logdir=dest/'logs';logdir.mkdir()

def run(name,args):
    with (logdir/(name+'.log')).open('w') as f:
        p=subprocess.run([sys.executable]+[str(a) for a in args],cwd=repo,stdout=f,stderr=subprocess.STDOUT)
    print(name, 'PASS' if p.returncode==0 else 'FAILED',flush=True)
    if p.returncode:raise SystemExit('See '+str(logdir/(name+'.log')))

run('protein_comparison',['run_pipeline.py','--skip-fetch'])
run('domain_annotation',['run_domain_annotation_v2.py'])
for rel in ['designs/junctions_v4.csv','designs/junctions_v4.fasta','results/tss_distance_distribution.csv','results/psyn_tss_assessment.md']:
    p=repo/rel
    if p.exists():p.rename(p.with_name(p.name+'.archived'))
run('sensor_generation',['run_junction_design_v4.py'])
run('transcript_start',['scripts/psyn_tss_report.py'])
run('baseline_consolidation',['scripts/consolidate_v4.py',dest/'consolidated_v4'])
run('operator_scan',['run_scan_pipeline.py','--skip-fetch'])
run('background_followup',['run_followup.py'])
run('p1_generation',[ROOT/'publication_code/build_revision_proposal.py',repo,dest/'p1'])
run('independent_audit',[ROOT/'publication_code/audit_candidates.py',dest/'p1',dest/'independent_audit'])
for name,script in [('sequence_tests','test_seqchecks.py'),('register_tests','verify_v4_register.py'),('report_tests','verify_crosstalk_report.py')]:
    run(name,['scripts/'+script])
comparisons=[]
with (ROOT/'reproduction/reproduction_comparisons.csv').open() as f:
    for row in csv.DictReader(f):
        stage=row['stage'];rel=row['file']
        p=(dest/'consolidated_v4'/rel if stage=='Baseline consolidation' else dest/'p1'/rel if stage=='P1 proposal generation' else repo/rel)
        digest=hashlib.sha256(p.read_bytes()).hexdigest()
        comparisons.append({**row,'current_sha256':digest,'current_status':'BYTE_IDENTICAL' if digest==row['original_sha256'] else 'DIFFERENT'})
with (dest/'comparison_results.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
different=[r for r in comparisons if r['current_status']!='BYTE_IDENTICAL']
print(f'{len(comparisons)-len(different)}/{len(comparisons)} selected files byte-identical')
print('Numerical-equivalence exception: inspect data/cache/intergenic_null_summary.json separately; see release notes.')
if different:raise SystemExit('Differences found; inspect comparison_results.csv. Do not equate differing output with a failed biological hypothesis.')
