"""Sequential isolated workers; no overlapping timed solves."""
from pathlib import Path
import sys,os,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
jobs=[]
for case in ['strong','moderate','weak']:
    for r in [2,3,4]:jobs.append((f'js_{case}_R{r}',['--method','js','--case',case,'--m',str(r)]))
    for m in [2,4,6,8]:jobs.append((f'anderson_{case}_m{m}',['--method','anderson','--case',case,'--m',str(m)]))
    for m in [2,3,4,5,6]:jobs.append((f'lti_{case}_m{m}',['--method','lti','--case',case,'--m',str(m)]))
    for m,t in [(2,1),(3,1),(3,2),(4,1),(4,2)]:jobs.append((f'eisert_{case}_m{m}_t{t}',['--method','eisert','--case',case,'--m',str(m),'--overlap',str(t)]))
    for last in [3,5,7]:jobs.append((f'kull_rg_{case}_k{last}_gauged',['--method','kull_rg','--case',case,'--last-k',str(last)]))
    for m in [3,4]:jobs.append((f'ffs_t0_{case}_m{m}',['--method','ffs_t0','--case',case,'--m',str(m)]))
(ROOT/'PRIMARY_JOBS.json').write_text(json.dumps(jobs,indent=2))
for name,args in jobs:
    out=ROOT/'results'/name
    if out.exists():print('EXISTS',name,flush=True);continue
    start=time.perf_counter()
    try:
        p=subprocess.run([sys.executable,'-B',str(ROOT/'code/benchmark.py'),*args,'--output',str(out)],capture_output=True,text=True,encoding='utf8',timeout=180)
        (out/'process.log').write_text(p.stdout+'\n'+p.stderr,encoding='utf8')
        print(name,p.returncode,round(time.perf_counter()-start,3),p.stdout.strip()[-450:],flush=True)
    except subprocess.TimeoutExpired as e:
        out.mkdir(exist_ok=True);(out/'failure.json').write_text(json.dumps(dict(status='TIMEOUT',seconds=180,args=args)))
        print('TIMEOUT',name,flush=True)
