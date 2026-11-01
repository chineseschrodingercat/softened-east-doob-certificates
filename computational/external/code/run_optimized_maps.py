"""Same general D=2 map optimization and k=5 RG level at all three cases."""
from pathlib import Path
import subprocess,sys,json,time
ROOT=Path(__file__).resolve().parents[1];jobs=[]
for case in ['strong','moderate','weak']:
    if not (ROOT/'lower_sources'/f'constant_{case}.json').exists():
        subprocess.run([sys.executable,'-B',str(ROOT/'code/constant_lower.py'),'--case',case],check=True)
    mp=ROOT/'maps'/f'optimized_{case}.json'
    if not mp.exists():subprocess.run([sys.executable,'-B',str(ROOT/'code/optimize_mps.py'),'--case',case,'--output',str(mp)],check=True)
    args=['--method','kull_rg','--case',case,'--last-k','5','--solver','SCS','--mps-map',str(mp)]
    name=f'kull_rg_{case}_k5_optimized_mps';out=ROOT/'results'/name;jobs.append((name,args))
    if out.exists():continue
    st=time.perf_counter()
    try:
        p=subprocess.run([sys.executable,'-B',str(ROOT/'code/benchmark.py'),*args,'--output',str(out)],capture_output=True,text=True,encoding='utf8',timeout=180)
        (out/'process.log').write_text(p.stdout+'\n'+p.stderr,encoding='utf8');print(name,p.returncode,time.perf_counter()-st,p.stdout.strip()[-400:],flush=True)
    except subprocess.TimeoutExpired:
        out.mkdir(exist_ok=True);(out/'failure.json').write_text(json.dumps(dict(status='TIMEOUT',seconds=180,args=args)));print('TIMEOUT',name,flush=True)
(ROOT/'OPTIMIZED_MAP_JOBS.json').write_text(json.dumps(jobs,indent=2))
