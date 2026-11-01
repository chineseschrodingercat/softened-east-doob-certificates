"""Outward interval enclosure of the finite-ring variational center Q_n.

Uses tr(T^n) and a local observable insertion on a de Bruijn transfer matrix.
T has 2^(R-1) states. No Perron vector/eigensolver is needed for finite n.
All interval operations use mpmath.iv; final endpoints and corrections use
exact rational arithmetic. Output scope is a SPECIFIC finite ring length.
"""
import argparse
import json
import hashlib
import time
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv
from candidate_inputs import candidate_path
from phase0 import ROOT, I, bounds, dec, rec, local_tables, direct_ring


def enclose_Q(data,n,dps):
    iv.dps=dps
    R=data['R'];m=2*R-1;D=1<<(R-1)
    assert n>=m
    c=F(data['c']);v=list(map(F,data['density']))
    te=[I(c if word>>(R-1) else 1-c)*iv.exp(I(2*v[word])) for word in range(1<<R)]
    zero=I(0);one=I(1)
    A=[[one if i==j else zero for j in range(D)] for i in range(D)]
    bridge=None
    for t in range(1,n+1):
        A=[[row[((head<<1)&(D-1))]*te[head<<1]+row[(((head<<1)|1)&(D-1))]*te[(head<<1)|1]
            for head in range(D)] for row in A]
        if t==n-R:bridge=A
    den=sum((A[i][i] for i in range(D)),zero)
    assert bounds(den)[0]>0 and bridge is not None
    phi=local_tables(data,True)[0]
    numerator=zero;mass=zero
    for word in range(1<<m):
        path=one
        for j in range(R):path*=te[(word>>j)&((1<<R)-1)]
        weight=path*bridge[word>>R][word&(D-1)]
        p=iv.mpf([I(phi[word][0]).a,I(phi[word][1]).b])
        numerator+=weight*p;mass+=weight
    # The inserted observable 1 must reproduce the partition function.
    mass_error=mass-den
    assert bounds(mass_error)[0]<=0<=bounds(mass_error)[1]
    return bounds(numerator/den),{'states':D,'exponent':n,'partition_function':rec(*bounds(den)),
                                'normalization_residual':rec(*bounds(mass_error))}


def interval_unit_check(dps):
    # Independent enumeration: directly sum pi(x) e^(2U(x)) and V(x).
    data={'R':3,'c':'0.3','kappa':'0.1','s':'-0.03',
          'density':['0','0.01','-0.03','0.02','0.01','-0.02','0.04','-0.01']}
    n=5;c=F(data['c']);k=F(data['kappa']);s=F(data['s']);v=list(map(F,data['density']))
    iv.dps=dps; U=[];pis=[]
    for x in range(1<<n):
        U.append(sum(v[sum(((x>>((i+j)%n))&1)<<j for j in range(3))] for i in range(n)))
        z=x.bit_count();pis.append(c**z*(1-c)**(n-z))
    den=I(0);num=I(0)
    for x in range(1<<n):
        w=I(pis[x])*iv.exp(I(2*U[x]));V=I(0)
        for i in range(n):
            bit=(x>>i)&1;prev=(x>>((i-1)%n))&1
            rate=(k+prev)*(c if bit==0 else 1-c)
            V+=I(rate)*(iv.exp(I(U[x^(1<<i)]-U[x]-s))-1)
        den+=w;num+=w*V/n
    direct=bounds(num/den);tm,_=enclose_Q(data,n,dps)
    assert max(direct[0],tm[0])<=min(direct[1],tm[1])
    return {'status':'PASS_INTERVAL_ENUMERATION_CROSSCHECK','n':n,'direct':rec(*direct),'transfer':rec(*tm)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--constants',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--dps',type=int,default=80)
    ap.add_argument('--case',default='');args=ap.parse_args()
    out=args.output.resolve()
    if not out.is_relative_to(ROOT):raise ValueError('Output must stay inside extension')
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use a new result directory')
    out.mkdir(parents=True,exist_ok=True)
    (out/'interval_crosscheck.json').write_text(json.dumps(interval_unit_check(args.dps),indent=2)+'\n')
    all_results=[]
    for path in sorted(args.constants.glob('*_constants.json')):
        const=json.loads(path.read_text())
        if args.case and args.case not in const['candidate']:continue
        p=candidate_path(const['candidate']);data=json.loads(p.read_text())
        assert const['candidate_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
        assert const['status']=='INTERVAL_CONSTANTS_VERIFIED'
        n=data['reference_N'];t0=time.perf_counter()
        (lo,hi),diagnostic=enclose_Q(data,n,args.dps)
        correction=F(const['correction_upper']);newlo=lo;newhi=hi+correction
        old=json.loads(p.with_name(p.stem+'_verified.json').read_text())
        oldlo=F(old['lower_density']);oldhi=F(old['upper_density'])
        intersection=(max(newlo,oldlo),min(newhi,oldhi))
        assert intersection[0]<=intersection[1]
        qr=rec(lo,hi);sr=rec(newlo,newhi)
        qw=F(qr['upper'])-F(qr['lower']);sw=F(sr['upper'])-F(sr['lower'])
        result={'status':'FINITE_N_VARIATIONAL_ENCLOSURE','candidate':p.name,
                'candidate_sha256':const['candidate_sha256'],'kappa':data['kappa'],'s':data['s'],'R':data['R'],
                'n':n,'interval_dps':args.dps,'Q_n':qr,'Q_width_upper':dec(qw,True),
                'correction_upper':const['correction_upper'],'lambda_n_over_n':sr,
                'spectral_width_upper':dec(sw,True),'original_common_graph_interval':rec(oldlo,oldhi),
                'original_width':dec(oldhi-oldlo,True),'old_width_over_new_width':float((oldhi-oldlo)/sw),
                'intersection_for_this_n_only':rec(*intersection),'transfer_check':diagnostic,
                'seconds':time.perf_counter()-t0,
                'scope':'Specific finite ring n shown here. The numerical correction is uniform, but Q_n and these endpoints depend on n. Theorem D and the finite-center transfer formula apply to this specified ring.'}
        (out/(p.stem+'_finite_Q.json')).write_text(json.dumps(result,indent=2)+'\n')
        all_results.append(result)
        print(json.dumps({k:result[k] for k in ['candidate','n','spectral_width_upper','old_width_over_new_width','seconds']}),flush=True)
    (out/'summary.json').write_text(json.dumps(all_results,indent=2)+'\n')
    (out/'runtime.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dps':args.dps},indent=2)+'\n')


if __name__=='__main__':main()
