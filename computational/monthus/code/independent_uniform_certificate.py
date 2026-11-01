"""Independently verify rational product-factor graph certificates.

Only input factor/potential strings are read. No discovery/verifier imports and
no floating-point arithmetic. Tuple words and full periodic products are used
as two independent constructions of each N=2R-1 residual.
"""
import argparse
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json
import sys


def decimal_out(q, digits=30, upper=False):
    q=Q(q); scale=10**digits
    k=-((-q.numerator*scale)//q.denominator) if upper else (q.numerator*scale)//q.denominator
    s=str(abs(k)).rjust(digits+1,'0')
    ans=('-' if k<0 else '')+s[:-digits]+'.'+s[-digits:]
    assert Q(ans)>=q if upper else Q(ans)<=q
    return ans


def record(lo,hi):
    return {'lower_rational':str(lo),'upper_rational':str(hi),
            'lower_decimal':decimal_out(lo),'upper_decimal':decimal_out(hi,upper=True)}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--candidate',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--claimed-lower',required=True)
    ap.add_argument('--claimed-upper',required=True)
    ap.add_argument('--local-record',type=Path,required=True)
    args=ap.parse_args()
    source=args.candidate.read_bytes()
    data=json.loads(source.decode('utf-8-sig'))
    R=int(data['R']); ell=2*R-1; center=R-1
    assert R>=2 and data['n_min']==ell
    assert data['bit_order']=='leftmost_is_least_significant'
    assert data['model']=='periodic_soft_east_activity'
    c,kappa,z=Q(data['c']),Q(data['kappa']),Q(data['tilt_factor'])
    assert (c,kappa,z)==(Q(1,2),Q(1),Q(1,2))
    f=list(map(Q,data['factor'])); gu=list(map(Q,data['z_upper'])); gl=list(map(Q,data['z_lower']))
    assert len(f)==2**R and len(gu)==len(gl)==2**(ell-1)
    assert min(f)>0
    # Tuple-keyed factor and vertex dictionaries read the stated little-endian convention.
    def word(i,length):return tuple(int(b) for b in format(i,f'0{length}b')[::-1])
    factor={word(i,R):value for i,value in enumerate(f)}
    upper_potential={word(i,ell-1):value for i,value in enumerate(gu)}
    lower_potential={word(i,ell-1):value for i,value in enumerate(gl)}
    lower_claim,upper_claim=Q(args.claimed_lower),Q(args.claimed_upper)
    assert lower_claim<upper_claim
    edges=[]; table={}; ratios={}
    for w in product((0,1),repeat=ell):
        changed=list(w);changed[center]=1-changed[center];changed=tuple(changed)
        before=after=Q(1)
        for start in range(R):
            before*=factor[w[start:start+R]]
            after*=factor[changed[start:start+R]]
        ratio=after/before
        rate=(kappa+w[center-1])*(c*(1-w[center])+(1-c)*w[center])
        q=rate*(z*ratio-1)
        tail,head=w[:-1],w[1:]
        lower=q+lower_potential[tail]-lower_potential[head]
        upper=q+upper_potential[tail]-upper_potential[head]
        assert lower_claim<=lower and upper<=upper_claim
        table[w]=q;ratios[w]=ratio
        edges.append({'word_x0_to_xlast':list(w),'flip_ratio':str(ratio),'q':str(q),
                      'lower_adjusted':str(lower),'upper_adjusted':str(upper),
                      'lower_slack':str(lower-lower_claim),'upper_slack':str(upper_claim-upper)})
    exact_lo=min(Q(e['lower_adjusted']) for e in edges)
    exact_hi=max(Q(e['upper_adjusted']) for e in edges)
    assert len(edges)==2**ell

    # Recompute h_n directly on the full N-bit ring with integer factor indices.
    # This does not use the tuple-word local product quotient above.
    n=ell; configs=list(product((0,1),repeat=n))
    def full_h(x):
        val=Q(1)
        for site in range(n):
            index=sum((2**offset)*x[(site+offset)%n] for offset in range(R))
            val*=f[index]
        return val
    h={x:full_h(x) for x in configs}
    assert min(h.values())>0
    residuals=[];flip_identity_checks=0
    for x in configs:
        full_residual=Q(0);local_residual=Q(0);upper_sum=lower_sum=Q(0)
        for site in range(n):
            y=list(x);y[site]=1-y[site];y=tuple(y)
            full_ratio=h[y]/h[x]
            rate=(kappa+x[(site-1)%n])*(c*(1-x[site])+(1-c)*x[site])
            full_residual+=rate*(z*full_ratio-1)
            # The centered 2R-1 block starts R-1 sites before this flip.
            w=tuple(x[(site-center+j)%n] for j in range(ell))
            assert full_ratio==ratios[w]
            flip_identity_checks+=1
            local_residual+=table[w]
            tail,head=w[:-1],w[1:]
            upper_sum+=table[w]+upper_potential[tail]-upper_potential[head]
            lower_sum+=table[w]+lower_potential[tail]-lower_potential[head]
        assert full_residual==local_residual==upper_sum==lower_sum
        assert n*lower_claim<=full_residual<=n*upper_claim
        residuals.append(full_residual/n)

    local=json.loads(args.local_record.read_text(encoding='utf-8'))
    Flo=Q(local['trial']['F']['rational_lower'])
    assert Flo>upper_claim
    gap=Flo-upper_claim
    out={
        'status':'PASS','candidate':str(args.candidate.resolve()),
        'candidate_sha256':hashlib.sha256(source).hexdigest(),
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python':sys.version,'dependencies':'Python standard library only',
        'R':R,'ell':ell,'n_min':ell,
        'checked_counts':{'positive_factors':len(f),'upper_potentials':len(gu),
            'lower_potentials':len(gl),'all_graph_edges':len(edges),
            'full_ring_states':len(configs),'individual_full_vs_local_flip_ratios':flip_identity_checks,
            'full_vs_local_and_both_telescoping_residuals':len(configs)},
        'uniform_bounds':{'claimed_lower':args.claimed_lower,'claimed_upper':args.claimed_upper,
                          'exact_graph_extrema':record(exact_lo,exact_hi)},
        'smallest_rational_edge_slacks':{'lower':str(exact_lo-lower_claim),
                                         'upper':str(upper_claim-exact_hi)},
        'N_min_full_product_CW':record(min(residuals),max(residuals)),
        'strict_trial_F_minus_uniform_upper':{'rational_lower':str(gap),
                                              'decimal_lower':decimal_out(gap)},
        'all_n_argument':[
          'Each positive factor defines h_n as the product over all cyclic length-R windows.',
          'For n>=2R-1, flipping site i changes precisely the R factors whose starts are i-R+1,...,i; their quotient is the verified local ratio on 2R-1 consecutive sites.',
          'Thus (L_s h_n)/h_n is the sum of the verified q over n overlapping length-(2R-1) words.',
          'Their length-(2R-2) tails and heads form a closed de Bruijn walk; both potential differences sum to zero.',
          'Every edge satisfies lower<=q+lower_potential(tail)-lower_potential(head), and q+upper_potential(tail)-upper_potential(head)<=upper.',
          'Hence n*lower<=L_s h_n/h_n<=n*upper at every configuration. Positivity and irreducibility give Collatz-Wielandt bounds on lambda_n/n.',
          'Consequently limsup_n(lambda_n/n)<=upper<F(trial), regardless of whether existence of the full limit is invoked.'
        ],
        'does_not_use':['floating discovery density','floating optimization claims','parent verifier output','primary certificate implementation'],
        'edges':edges
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.output.exists():raise FileExistsError('Refusing to overwrite independent result')
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:out[k] for k in ['status','candidate_sha256','checked_counts','uniform_bounds','strict_trial_F_minus_uniform_upper']},indent=2))
    print('Saved:',args.output)


if __name__=='__main__':
    if not __debug__: raise RuntimeError('Do not disable certificate assertions')
    main()
