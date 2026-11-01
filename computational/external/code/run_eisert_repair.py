"""Apply the same second-solver repair to every added overlap-three program."""
from pathlib import Path
import sys,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];jobs=[]
for case in ['strong','moderate','weak']:
    for m,t in [(3,3),(4,3),(6,3)]:
        jobs.append((f'eisert_{case}_m{m}_t{t}_scs',['--method','eisert','--case',case,'--m',str(m),'--overlap',str(t),'--solver','SCS']))
(ROOT/'EISERT_REPAIR_JOBS.json').write_text(json.dumps(jobs,indent=2))
for name,args in jobs:
    out=ROOT/'results'/name
    if out.exists():continue
    st=time.perf_counter()
    try:
        p=subprocess.run([sys.executable,'-B',str(ROOT/'code/benchmark.py'),*args,'--output',str(out)],capture_output=True,text=True,encoding='utf8',timeout=180)
        (out/'process.log').write_text(p.stdout+'\n'+p.stderr,encoding='utf8')
        print(name,p.returncode,round(time.perf_counter()-st,3),p.stdout.strip()[-400:],flush=True)
    except subprocess.TimeoutExpired:
        out.mkdir(exist_ok=True);(out/'failure.json').write_text(json.dumps(dict(status='TIMEOUT',seconds=180,args=args)));print('TIMEOUT',name,flush=True)
