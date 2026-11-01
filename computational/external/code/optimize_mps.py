"""Generate general real D=2 MPS maps; no floating energy is a certificate."""
from pathlib import Path
import sys,time,json,argparse,platform
import numpy as np
from scipy.linalg import eig
from scipy.optimize import minimize
from fractions import Fraction as F
import psutil
import exact_sdp as ex
from benchmark import CASES,ROOT

def energy(z,h):
    A=np.asarray(z).reshape(2,2,2);A=A/np.linalg.norm(A)
    T=sum(np.kron(a,a) for a in A)
    vals,l,r=eig(T,left=True,right=True);j=int(np.argmax(vals.real));rho=vals[j].real
    if abs(vals[j].imag)>1e-8 or rho<=0:return 1e3
    lv=l[:,j].real;rv=r[:,j].real;den=lv@rv
    if abs(den)<1e-12:return 1e3
    B=np.zeros((4,4))
    for x in range(4):
        for y in range(4):B+=h[x,y]*np.kron(A[x//2]@A[x%2],A[y//2]@A[y%2])
    return float((lv@B@rv)/(rho*rho*den))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',choices=CASES,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=a.output.resolve()
    if out.exists() or not out.is_relative_to(ROOT):raise ValueError('fresh map path inside external required')
    st=time.perf_counter()
    tr=json.loads((ROOT/'results'/f'js_{a.case}_R2'/'trial.json').read_text());v=np.array(tr['density'],float)
    old=json.loads((ROOT/'results'/f'js_{a.case}_R2'/'result.json').read_text())
    k,s=CASES[a.case];h=np.asarray(ex.local_h('.3',k,s)[0],float);A=np.zeros((2,2,2))
    for x in range(2):
        for y in range(2):A[x,x,y]=np.sqrt(.3 if x else .7)*np.exp(v[x+2*y])
    A/=np.linalg.norm(A);rng=np.random.default_rng(621);starts=[A.ravel(),(A+.02*rng.standard_normal(A.shape)).ravel(),(A+.15*rng.standard_normal(A.shape)).ravel()]
    logs=[];solutions=[(energy(starts[0],h),starts[0])]
    for start in starts:
        sol=minimize(energy,start,args=(h,),method='BFGS',options={'gtol':1e-11,'maxiter':300})
        logs.append(dict(status=int(sol.status),message=str(sol.message),nit=int(sol.nit),energy_diagnostic=float(sol.fun)))
        solutions.append((float(sol.fun),sol.x))
    val,z=min(solutions,key=lambda q:q[0]);A=z.reshape(2,2,2);A/=np.linalg.norm(A)
    elapsed=time.perf_counter()-st
    data=dict(case=a.case,bond_dimension=2,mps=[ex.dumpmat(np.array([[F(format(float(x),'.14g')) for x in row] for row in q],dtype=object)) for q in A],
       energy_diagnostic=val,initial_energy_diagnostic=solutions[0][0],optimizer_logs=logs,
       map_optimization_seconds=elapsed,initial_trial_generation_seconds=old['candidate_seconds'],
       total_generation_seconds=elapsed+old['candidate_seconds'],process_peak_memory_bytes=psutil.Process().memory_info().peak_wset,
       use='Map generation only; neither MPS energy nor optimizer status is a certified bound',python=sys.version)
    out.parent.mkdir(exist_ok=True);out.write_text(json.dumps(data,indent=2));print(json.dumps({k:data[k] for k in ['case','energy_diagnostic','initial_energy_diagnostic','total_generation_seconds']}),flush=True)
if __name__=='__main__':main()
