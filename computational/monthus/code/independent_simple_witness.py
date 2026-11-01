"""Exact, human-checkable local trial and its positive N=3 realization."""
import argparse
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib
import json
import sys


def sqrt_interval(q):
    scale=10**40
    k=isqrt((q.numerator*scale*scale)//q.denominator)
    lo,hi=Q(k,scale),Q(k+1,scale)
    assert lo*lo<=q<hi*hi
    return lo,hi


def fixed(q,upper=False):
    scale=10**30
    k=-((-q.numerator*scale)//q.denominator) if upper else (q.numerator*scale)//q.denominator
    s=str(abs(k)).rjust(31,'0')
    v=('-' if k<0 else '')+s[:-30]+'.'+s[-30:]
    assert Q(v)>=q if upper else Q(v)<=q
    return v


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    a,b,d=Q(11,100),Q(187,1000),Q(129,250)
    assert a+2*b+d==1 and min(a,b,d)>0
    # Two four-decimal numbers already prove the claimed strict lower bound.
    pab,pbd=Q('0.1434'),Q('0.3106')
    assert pab*pab<a*b and pbd*pbd<b*d
    penalty=a+Q(3,2)*b+d/2
    coarse=pab+pbd/2-penalty
    assert penalty==Q('0.6485') and coarse==Q('-0.3498')
    abl,abu=sqrt_interval(a*b);bdl,bdu=sqrt_interval(b*d)
    lo,hi=abl+bdl/2-penalty,abu+bdu/2-penalty
    assert lo>coarse
    weights={0:Q('0.384'),1:Q('0.132'),2:Q('0.055'),3:Q('0.055')}
    rho=[weights[i.bit_count()] for i in range(8)]
    assert sum(rho)==1 and min(rho)>0
    expected={(1,1):a,(1,0):b,(0,1):b,(0,0):d}
    marginals=[]
    for j in range(3):
        observed={}
        for pair in expected:
            observed[pair]=sum(rho[i] for i in range(8)
                              if ((i>>((j-1)%3))&1,(i>>j)&1)==pair)
        assert observed==expected
        marginals.append({str(k):str(v) for k,v in observed.items()})
    assert a*d-b*b==Q(21791,10**6)>0
    out={'status':'PASS','python':sys.version,'dependencies':'standard library only',
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'trial':{'a':str(a),'b':str(b),'d':str(d),
                  'F':{'rational_lower':str(lo),'rational_upper':str(hi),
                       'decimal_lower':fixed(lo),'decimal_upper':fixed(hi,True)},
                  'strict_coarse_lower':'-0.3498'},
         'human_square_inequalities':{
                  '0.1434_squared':str(pab*pab),'ab':str(a*b),
                  '0.3106_squared':str(pbd*pbd),'bd':str(b*d),
                  'diagonal_penalty':str(penalty)},
         'positive_global_density_binary_integer_order':[str(v) for v in rho],
         'all_3_pair_marginals':marginals,
         'ad_minus_b_squared':str(a*d-b*b)}
    if args.output.exists():raise FileExistsError('Refusing overwrite')
    args.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    if not __debug__: raise RuntimeError('Do not disable certificate assertions')
    main()
