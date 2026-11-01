"""One-command complete frozen release replay; creates fresh outputs only."""
from pathlib import Path
import argparse,re,subprocess,sys,json,time,os
ROOT=Path(__file__).resolve().parents[1]
def main():
    if not __debug__:raise RuntimeError('Python -O is forbidden')
    ap=argparse.ArgumentParser();ap.add_argument('--tag',default='replay_complete');a=ap.parse_args()
    if not re.fullmatch('[A-Za-z0-9_-]+',a.tag):raise ValueError('tag must contain only letters, numbers, underscore or hyphen')
    logs=ROOT/'replay_logs'/a.tag
    outputs=[logs,ROOT/'external'/a.tag,ROOT/'certificates'/a.tag,ROOT/'monthus'/(a.tag+'_primary'),ROOT/'monthus'/(a.tag+'_independent'),ROOT/'baseline'/a.tag]
    if any(p.exists() for p in outputs):raise FileExistsError('Choose a fresh --tag; no existing output is overwritten')
    logs.mkdir(parents=True)
    jobs=[('manifest',['code/check_manifest.py']),('external',['external/code/replay_external.py','--output',a.tag]),
      ('all_original',['certificates/code/replay_submission.py','--scope','all','--output',a.tag,'--dps','120']),
      ('monthus_primary',['monthus/code/replay_monthus.py','--output',a.tag+'_primary']),
      ('monthus_independent',['monthus/code/replay_independent.py','--output',a.tag+'_independent']),
      ('historical_cost_accuracy',['certificates/code/check_cost_accuracy.py']),
      ('closed_form',['certificates/checks/closed_form_checks.py']),('reversible',['certificates/checks/reversible_checks.py']),
      ('baseline',['baseline/code/reproduce_v14.py','--output','baseline/'+a.tag])]
    rows=[];env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONIOENCODING='utf-8')
    for name,args in jobs:
        st=time.perf_counter()
        with (logs/(name+'.log')).open('w',encoding='utf8') as f:
            p=subprocess.run([sys.executable,'-B',*args],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        rows.append(dict(check=name,returncode=p.returncode,seconds=time.perf_counter()-st,log=name+'.log'))
        report=dict(status='FAILED' if p.returncode else 'RUNNING',checks=rows)
        (logs/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(name,'PASS' if not p.returncode else 'FAIL',flush=True)
        if p.returncode:raise RuntimeError('Replay failed; inspect '+str(logs/(name+'.log')))
    (logs/'REPORT.json').write_text(json.dumps(dict(status='PASS',checks=rows),indent=2)+'\n');print('All frozen-release checks passed.',flush=True)
if __name__=='__main__':main()
