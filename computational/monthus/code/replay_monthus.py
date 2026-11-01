"""Reproduce the core Monthus counterexample with Python standard library only.
Run from any directory: python -B monthus/code/replay_monthus.py
Float optimization, diagonalization, and logarithms are not used in validation.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,gcd
from functools import reduce
import json
from verify_uniform_certificate import verify,decimal
ROOT=Path(__file__).resolve().parents[1]

def square_root_interval(q,digits=30):
    s=10**digits
    n=isqrt(q.numerator*s*s//q.denominator)
    lo,hi=F(n,s),F(n+1,s)
    assert lo*lo<=q<hi*hi
    return lo,hi

def run():
    # B = 4 L_{3,log(2)} + 12 I: every entry is an integer >= 0.
    n=3;B=[[0]*8 for _ in range(8)]
    for x in range(8):
        total=0
        for i in range(3):
            prev=(x>>((i-1)%3))&1
            B[x][x^(1<<i)]=1+prev
            total+=2*(1+prev)
        B[x][x]=12-total
    assert all(v>=0 for row in B for v in row)
    assert all(B[x][y]==B[y][x] for x in range(8) for y in range(8))
    v=[1]*8
    for _ in range(100):
        v=[sum(rj*vj for rj,vj in zip(row,v)) for row in B]
        g=reduce(gcd,v);v=[a//g for a in v]
    assert min(v)>0
    q=[F(sum(rj*vj for rj,vj in zip(row,v)),v[i]) for i,row in enumerate(B)]
    lo=(min(q)-12)/12;hi=(max(q)-12)/12
    assert F('-0.350706886146')<lo<=hi<F('-0.350706886143')
    finite={'matrix_convention':'B=4*L_3+12*I','B':B,'positive_integer_witness':[str(t) for t in v],'power_steps':100,'lambda3_over_3_lower':decimal(lo,30),'lambda3_over_3_upper':decimal(hi,30,True),'lower_exact':str(lo),'upper_exact':str(hi)}
    a,b,d=F('0.110'),F('0.187'),F('0.516')
    assert a+2*b+d==1 and min(a,b,d)>0
    l1,u1=square_root_interval(a*b);l2,u2=square_root_interval(b*d)
    fl=l1+l2/2-a-3*b/2-d/2;fh=u1+u2/2-a-3*b/2-d/2
    assert fl>F('-0.3498')>hi
    probabilities=[d-b+a/2]+[b-a/2 if x.bit_count()==1 else a/2 for x in range(1,8)]
    assert sum(probabilities)==1 and min(probabilities)>0
    for i in range(3):
        pair={(left,center):sum(p for x,p in enumerate(probabilities) if ((x>>((i-1)%3))&1)==left and ((x>>i)&1)==center) for left in (0,1) for center in (0,1)}
        assert pair=={(0,0):d,(0,1):b,(1,0):b,(1,1):a}
    assert a*d!=b*b
    trial={'pair':{'rho11':str(a),'rho10':str(b),'rho01':str(b),'rho00':str(d)},'lower':decimal(fl,30),'upper':decimal(fh,30,True),'positive_global_probabilities':[str(v) for v in probabilities],'ad_minus_b2':str(a*d-b*b)}
    short=verify(ROOT/'results/uniform_R3_simple_candidate.json',True)
    sharper=verify(ROOT/'results/uniform_R4_candidate.json',True)
    assert F(short['upper'])<F('-0.3499')<F('-0.3498')<fl
    assert F(sharper['upper'])<F(short['upper'])
    return {'status':'PASS','arithmetic':'Python integers/Fraction; square roots enclosed by integer squares','finite_ring':finite,'simple_trial':trial,'R3_all_n_at_least_5':short,'R4_all_n_at_least_7':sharper,'conclusion':'Local closure maximum exceeds the true finite-N spectral density and also the limsup as N tends to infinity at c=1/2,kappa=1,exp(-s)=1/2. No locality theorem or extrapolation is used.'}

if __name__=='__main__':
    import argparse
    from package_integrity import check_manifest,new_output
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='replay_primary')
    args=ap.parse_args();integrity=check_manifest();dest=new_output(args.output)
    result=run()
    expected=json.loads((ROOT/'results/published_results.json').read_text(encoding='utf-8-sig'))
    assert result==expected,'Exact replay differs from published record'
    (dest/'exact_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    report={'status':'PASS','arithmetic':'Python standard-library integers and Fraction','exact_published_record_match':True,'finite_ring':3,'uniform_edges':[32,128],'uniform_n_min':[5,7],'package_integrity':integrity}
    (dest/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
