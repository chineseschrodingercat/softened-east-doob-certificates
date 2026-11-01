"""Run pinned unmodified Kull MATLAB dual SDP; exact certification in Python.

Generation needs a work root with author_code and matlab_runtime, as documented.
Frozen certificate replay requires neither MATLAB nor these dependencies.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys,json,time,subprocess,argparse,traceback
from fractions import Fraction as F
import numpy as np,psutil
from scipy.io import savemat
import exact_sdp as ex
ROOT=Path(__file__).resolve().parents[1]
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')
def worker(a):
    name=f'kull_native_{a.case}_k{a.last_k}_{a.maps}'
    out=ROOT/'results'/name;out.mkdir(exist_ok=False)
    start=time.perf_counter();peak=0
    try:
        source=ROOT/'results'/f'kull_rg_{a.case}_k5_optimized_mps'/'candidate.json' if a.maps=='optimized' else ROOT/'results'/f'kull_rg_{a.case}_k3_scs'/'candidate.json'
        if not source.exists():
            source=next(p for p in (ROOT/'results').glob(f'*{a.case}*k3*/candidate.json') if json.loads(p.read_text()).get('method')=='kull_rg')
        old=json.loads(source.read_text());spec={k:old[k] for k in ['case','c','kappa','s','mps']}
        spec.update(method='kull_rg',implementation='unmodified MATLAB author dual SDP',k0=3,last_k=a.last_k,n_min=a.last_k+2,m=4,solver='SDPT3 standard model 1',author_commit='2e9015fff5d9bc5b170cdc6cee98fbbb928decda')
        A=[ex.loadmat(x) for x in spec['mps']];V,_,_=ex.mps_maps(A,3)
        U,sv,_=np.linalg.svd(np.asarray(V,float),full_matrices=False);W=np.diag(1/sv)@U.T
        spec['gauge']=ex.dumpmat(np.array([[F(format(float(x),'.14g')) for x in row] for row in W],dtype=object))
        V,L,R=ex.gauged_maps(A,3,spec['gauge']);H,_,meta=ex.local_h(spec['c'],spec['kappa'],spec['s'])
        savemat(out/'input.mat',{k:np.asarray(v,float) for k,v in dict(H=H,V=V,L=L,R=R,last_k=a.last_k).items()})
        preparation=time.perf_counter()-start
        env=os.environ.copy();env.update(EXTERNAL_WORKROOT=str(a.work_root).replace('\\','/'),EXTERNAL_NATIVE_OUT=str(out).replace('\\','/'))
        cmd=[a.matlab,'-wait','-singleCompThread','-batch',"run('"+str(Path(__file__).with_name('native_kull.m')).replace('\\','/')+"')"]
        t=time.perf_counter()
        with (out/'matlab_stdout.log').open('w',encoding='utf8') as log:
            p=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
            while p.poll() is None:
                try:
                    processes=[psutil.Process(),psutil.Process(p.pid)]+psutil.Process(p.pid).children(recursive=True)
                    peak=max(peak,sum(q.memory_info().rss for q in processes if q.is_running()))
                except (psutil.NoSuchProcess,psutil.AccessDenied):pass
                if time.perf_counter()-t>180:
                    for q in psutil.Process(p.pid).children(recursive=True):q.kill()
                    p.kill();p.wait();raise TimeoutError('native worker exceeded 180 seconds')
                time.sleep(.05)
        wall=time.perf_counter()-t
        raw=json.loads((out/'matlab_output.json').read_text())
        if p.returncode:raise RuntimeError('MATLAB exit '+str(p.returncode))
        spec['duals']=[ex.dumpmat(ex.rat(d)) for d in raw['duals']]
        t=time.perf_counter();mats,_,_,_=ex.slacks(spec);spec['psd_certificates']=[ex.make_psd_certificate(x) for x in mats];factor=time.perf_counter()-t
        write(out/'candidate.json',spec)
        t=time.perf_counter();verified=ex.verify_sdp(spec);verification=time.perf_counter()-t
        if a.maps=='optimized':
            mp=json.loads((ROOT/'maps'/f'optimized_{a.case}.json').read_text());mpsgen=mp['total_generation_seconds'];mpspeak=mp['process_peak_memory_bytes']
        else:
            mp=json.loads((ROOT/'results'/f'js_{a.case}_R2'/'result.json').read_text());mpsgen=mp['candidate_seconds'];mpspeak=mp['process_peak_memory_bytes']
        gen=preparation+raw['candidate_seconds']+factor+mpsgen
        result=dict(method='kull_native',case=a.case,m=4,last_k=a.last_k,overlap=1,**verified,
          candidate_seconds=gen,verification_seconds=verification,total_seconds=gen+verification,
          shared_mps_generation_seconds=mpsgen,mps_strategy='general_real_D2_optimized' if a.maps=='optimized' else 'JS_range_two',
          raw_primal_energy=raw['energy'],raw_energy_kind='floating dual objective, diagnostic only',
          solver_status='native_problem_'+str(raw['problem']),yalmip_problem=raw['problem'],solver_info=raw['info'],solver_seconds=raw['solver_seconds'],
          matlab=raw['matlab'],yalmip=raw['yalmip'],coefficient=meta,
          process_peak_memory_bytes=max(peak,mpspeak,psutil.Process().memory_info().peak_wset),
          native_process_wall_seconds=wall,matlab_startup_and_io_seconds=max(0,wall-raw['candidate_seconds']),
          timing_scope='Fresh MATLAB process; mathematical core plus Python preparation, rational factor construction, verification and MPS preparation. MATLAB startup and I/O separately recorded.',
          memory_scope='Maximum sampled aggregate Python worker plus MATLAB launcher and descendants; 50 ms sampling, including startup. Also includes peak Python worker and separate MPS preparation. Could miss shorter transient peaks.')
        write(out/'result.json',result);print(json.dumps({k:result[k] for k in ['case','last_k','mps_strategy','upper_density','total_seconds','solver_status']}),flush=True)
    except Exception as e:
        write(out/'failure.json',dict(method='kull_native',case=a.case,last_k=a.last_k,maps=a.maps,status='FAILED',error=repr(e),traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-start,process_peak_memory_bytes=peak));raise
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--worker',action='store_true');ap.add_argument('--case');ap.add_argument('--last-k',type=int,default=5);ap.add_argument('--maps',default='block');ap.add_argument('--work-root',type=Path,required=True);ap.add_argument('--matlab',default=r'C:\Program Files\MATLAB\R2026b\bin\matlab.exe');a=ap.parse_args()
    if a.worker:worker(a);return
    for case in ['strong','moderate','weak']:
        for k,maps in [(3,'block'),(5,'block'),(7,'block'),(5,'optimized')]:
            cmd=[sys.executable,__file__,'--worker','--case',case,'--last-k',str(k),'--maps',maps,'--work-root',str(a.work_root),'--matlab',a.matlab]
            r=subprocess.run(cmd);print(case,k,maps,'exit',r.returncode,flush=True)
if __name__=='__main__':main()
