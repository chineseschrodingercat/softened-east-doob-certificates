"""Exact interval Rayleigh value for h=1, used as a cheap lower-endpoint option."""
from pathlib import Path
from fractions import Fraction as F
import time,json,argparse,sys
from mpmath import iv
import psutil
ROOT=Path(__file__).resolve().parents[1]
CASES={'strong':('1','0.02'),'moderate':('0.1','0.01'),'weak':('0.01','0.001')}
def frac(t):
    sign,m,e,b=t
    if b<0:raise ArithmeticError('nonfinite endpoint')
    return (-1 if sign else 1)*F(m)*F(2)**e
def I(x):
    x=F(x);return iv.mpf(x.numerator)/x.denominator
def dec(x,up=False):
    y=F(x)*10**35;n=-((-y.numerator)//y.denominator) if up else y.numerator//y.denominator
    return ('-' if n<0 else '')+str(abs(n)//10**35)+'.'+str(abs(n)%10**35).zfill(35)
def verify(d):
    c,k,s=map(F,[d['c'],d['kappa'],d['s']]);assert 0<c<1 and k>0 and s>=0;iv.dps=100
    q=2*I(c)*I(1-c)*(I(k)+I(c))*(iv.exp(-I(s))-1)
    return {'lower':dec(frac(q._mpi_[0])),'upper':dec(frac(q._mpi_[1]),True)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',choices=CASES,required=True);a=ap.parse_args()
    out=ROOT/'lower_sources'/f'constant_{a.case}.json'
    if out.exists():raise FileExistsError(out)
    st=time.perf_counter();k,s=CASES[a.case];d=dict(case=a.case,c='0.3',kappa=k,s=s,n_min=3,targets=['limit','n18']);gen=time.perf_counter()-st
    st=time.perf_counter();q=verify(d);ver=time.perf_counter()-st
    d.update(Q=q,candidate_seconds=gen,verification_seconds=ver,total_seconds=gen+ver,
        process_peak_memory_bytes=psutil.Process().memory_info().peak_wset,
        derivation='E_pi[(L_s 1)/n] = 2 c (1-c) (kappa+c) (exp(-s)-1); Rayleigh lower value, plus universal spectral upper zero')
    out.parent.mkdir(exist_ok=True);out.write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
if __name__=='__main__':main()
