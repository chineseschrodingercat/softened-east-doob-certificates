"""Replayable stdlib-only rational certificate, independent of candidate generation.
All arithmetic after JSON parsing uses Fraction, with outward decimal rounding.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json

def decimal(q,places=15,upper=False):
    q=q*10**places
    n=-((-q.numerator)//q.denominator) if upper else q.numerator//q.denominator
    sign='-' if n<0 else '';n=abs(n)
    return f'{sign}{n//10**places}.{n%10**places:0{places}d}'

def verify(path,exhaustive=False):
    raw=path.read_bytes();d=json.loads(raw)
    assert d['model']=='periodic_soft_east_activity'
    assert d['bit_order']=='leftmost_is_least_significant'
    R=d['R'];m=2*R-1;V=1<<(m-1)
    assert d['n_min']==m and R>=2
    p=list(map(F,d['factor']));lo=list(map(F,d['z_lower']));hi=list(map(F,d['z_upper']))
    assert len(p)==1<<R and min(p)>0 and len(lo)==len(hi)==V
    c,k,t=map(F,(d['c'],d['kappa'],d['tilt_factor']))
    assert 0<c<1 and k>0 and t>0
    phi=[];upper_edges=[];lower_edges=[]
    for word in range(2*V):
        ratio=F(1)
        for start in range(R):
            block=(word>>start)&((1<<R)-1)
            ratio*=p[block^(1<<(R-1-start))]/p[block]
        bit=(word>>(R-1))&1;prev=(word>>(R-2))&1
        rate=(k+prev)*(c if bit==0 else 1-c)
        f=rate*(t*ratio-1);phi.append(f)
        tail=word&(V-1);head=word>>1
        lower_edges.append(f+lo[tail]-lo[head]);upper_edges.append(f+hi[tail]-hi[head])
    a=min(lower_edges);b=max(upper_edges)
    out={'status':'PASS_EXACT_RATIONAL_ALL_EDGES','candidate_sha256':hashlib.sha256(raw).hexdigest(),'R':R,'n_min':m,'edges':2*V,'lower_exact':str(a),'upper_exact':str(b),'lower':decimal(a),'upper':decimal(b,upper=True),'scope':'lambda_n(log(2))/n for all integer n>=n_min; same interval for any thermodynamic limit','no_analytic_locality_assumption':True}
    if exhaustive:
        n=m;h=[]
        for x in range(1<<n):
            val=F(1)
            for i in range(n):
                block=sum(((x>>((i+j)%n))&1)<<j for j in range(R))
                val*=p[block]
            h.append(val)
        res=[]
        for x in range(1<<n):
            actual=sum(((k+((x>>((i-1)%n))&1))*(c if not((x>>i)&1) else 1-c)*(t*h[x^(1<<i)]/h[x]-1)) for i in range(n))
            graph=F(0)
            for i in range(n):
                word=sum(((x>>((i+j)%n))&1)<<j for j in range(m))
                graph+=phi[word]
            assert actual==graph
            assert n*a<=actual<=n*b
            res.append(actual/F(n))
        out['direct_ring_identity_check']={'n':n,'states':1<<n,'all_exact':True,'min_residual':decimal(min(res)),'max_residual':decimal(max(res),upper=True)}
    return out

if __name__=='__main__':
    if not __debug__: raise RuntimeError('Do not disable certificate assertions')
    parser=argparse.ArgumentParser();parser.add_argument('candidate',type=Path);parser.add_argument('--exhaustive',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args();out=verify(args.candidate,args.exhaustive)
    if args.output:
        if args.output.exists():raise FileExistsError(args.output)
        args.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('lower_exact','upper_exact')},indent=2))
