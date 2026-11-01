"""Independent exact N=3/Monthus local-relaxation audit, standard library only.

Discovery and certification both use integer power iteration plus exact rational
inequalities.  No floating-point optimizer or eigensolver enters the certificate.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from functools import reduce
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json
import platform
import sys


def rat(q):
    return str(Q(q))


def fixed_integer(k, digits):
    sign = '-' if k < 0 else ''
    s = str(abs(k)).rjust(digits + 1, '0')
    return sign + s[:-digits] + '.' + s[-digits:]


def outward(q, digits, upper=False):
    q = Q(q)
    scale = 10 ** digits
    floor = (q.numerator * scale) // q.denominator
    k = -((-q.numerator * scale) // q.denominator) if upper else floor
    return fixed_integer(k, digits)


def interval_record(lo, hi, digits=30):
    assert lo <= hi
    dl, du = outward(lo, digits), outward(hi, digits, True)
    assert Q(dl) <= lo <= hi <= Q(du)
    return {'rational_lower': rat(lo), 'rational_upper': rat(hi),
            'decimal_lower': dl, 'decimal_upper': du,
            'decimal_places': digits}


def sqrt_bounds(q, digits=45):
    q = Q(q)
    assert q >= 0
    scale = 10 ** digits
    root_floor = isqrt((q.numerator * scale * scale) // q.denominator)
    lo = Q(root_floor, scale)
    hi = lo if lo * lo == q else Q(root_floor + 1, scale)
    assert 0 <= lo <= hi and lo * lo <= q <= hi * hi
    return lo, hi


def witness(M, diagonal_shift, iterations=250, scale=10**20):
    n = len(M)
    B = [[M[i][j] + (diagonal_shift if i == j else 0)
          for j in range(n)] for i in range(n)]
    assert all(v >= 0 for row in B for v in row)
    v = [1] * n
    for _ in range(iterations):
        v = [sum(b*x for b,x in zip(row,v)) for row in B]
        common = reduce(gcd, v)
        v = [x // common for x in v]
    mx = max(v)
    v = [scale * x // mx for x in v]
    assert min(v) > 0
    return v


def cw(M, v, denominator):
    ratios = [Q(sum(a*x for a,x in zip(row,v)), denominator*v[i])
              for i,row in enumerate(M)]
    return min(ratios), max(ratios), ratios


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    n = 3
    states = [tuple((x >> i) & 1 for i in range(n)) for x in range(1 << n)]
    # Explicitly assemble the row-action generator from rational flip rates.
    L = [[Q(0) for _ in states] for _ in states]
    L_monthus = [[Q(0) for _ in states] for _ in states]
    untilted = [[Q(0) for _ in states] for _ in states]
    c, kappa, z = Q(1,2), Q(1), Q(1,2)
    eps = kappa/(1+kappa)
    for x,bits in enumerate(states):
        for i in range(n):
            left, own = bits[(i-1)%n], bits[i]
            base = c*(1-own)+(1-c)*own
            r = (kappa+left)*base
            rm = (left+eps*(1-left))*base
            y = x^(1 << i)
            L[x][x] -= r
            L[x][y] += z*r
            L_monthus[x][x] -= rm
            L_monthus[x][y] += z*rm
            untilted[x][x] -= r
            untilted[x][y] += r
    assert all(sum(row) == 0 for row in untilted)
    assert all(L[i][j] == L[j][i] for i in range(8) for j in range(8))
    assert all(L[i][j] == 2*L_monthus[i][j] for i in range(8) for j in range(8))
    M = [[int(4*v) for v in row] for row in L]
    assert all(Q(M[i][j],4) == L[i][j] for i in range(8) for j in range(8))
    h = witness(M, 12)
    lam_lo, lam_hi, ratios = cw(M, h, 12)  # 4L / (4*3)
    assert Q('-0.350706886146') < lam_lo <= lam_hi < Q('-0.350706886143')

    # A rational point chosen independently near the positive local optimizer.
    a, b, d = Q(110309,10**6), Q(186559,10**6), Q(516573,10**6)
    assert min(a,b,d)>0 and a+2*b+d==1 and b<Q(1,3)
    ablo, abhi = sqrt_bounds(a*b)
    bdlo, bdhi = sqrt_bounds(b*d)
    diag = a+Q(3,2)*b+Q(1,2)*d
    F_lo = ablo+bdlo/2-diag
    F_hi = abhi+bdhi/2-diag
    assert F_lo > lam_hi
    gap = F_lo-lam_hi
    assert gap > Q('0.000945858229')

    # Exact local maximum: F is a generalized Rayleigh quotient, D=diag(1,2,1).
    # A=D^{-1}K is similar to symmetric D^{-1/2}KD^{-1/2}.
    A8 = [[-8,4,0],[2,-6,1],[0,2,-4]]
    v = witness(A8, 8)
    local_lo, local_hi, local_ratios = cw(A8, v, 8)
    assert F_hi < local_lo
    # Pure integer discovery of the six-decimal point from this positive vector.
    norm = v[0]**2+2*v[1]**2+v[2]**2
    rounded_a = (2*10**6*v[0]**2+norm)//(2*norm)
    rounded_b = (2*10**6*v[1]**2+norm)//(2*norm)
    assert (rounded_a,rounded_b)==(110309,186559)

    # A full-support three-spin distribution with the EXACT trial pair marginals.
    orbit_mass_per_state = {0:Q(380014,10**6),1:Q(136559,10**6),
                           2:Q(50000,10**6),3:Q(60309,10**6)}
    rho = [orbit_mass_per_state[sum(bits)] for bits in states]
    assert sum(rho)==1 and min(rho)>0
    pair_target = {(1,1):a,(1,0):b,(0,1):b,(0,0):d}
    pair_results = []
    for i in range(n):
        measured = {pair:sum(rho[x] for x,bits in enumerate(states)
                    if (bits[(i-1)%n],bits[i])==pair) for pair in pair_target}
        assert measured==pair_target
        pair_results.append({str(pair):rat(value) for pair,value in measured.items()})

    # The same pair is also an admissible stationary Markov-chain marginal.
    pi = [a+b,b+d]  # state order +,-
    P = [[a/pi[0],b/pi[0]],[b/pi[1],d/pi[1]]]
    assert all(sum(row)==1 for row in P)
    assert all(sum(pi[i]*P[i][j] for i in range(2))==pi[j] for j in range(2))
    assert [[pi[i]*P[i][j] for j in range(2)] for i in range(2)]==[[a,b],[b,d]]

    # Equality in BOTH aggregate Cauchy inequalities would require ratios
    # t+=a/b and t-=b/d for every context.  The 000 -> 110 flip square
    # gives t-*t+ along one route but t-^2 along the other.
    tplus, tminus = a/b,b/d
    path12, path21 = tminus*tplus,tminus*tminus
    assert a*d!=b*b and path12!=path21

    # Direct rational-sqrt enclosure of the Rayleigh value of the explicit rho.
    rlo = rhi = Q(0)
    Dlo, Dhi = {0:Q(0),1:Q(0)}, {0:Q(0),1:Q(0)}
    for x,bits in enumerate(states):
        rlo += L[x][x]*rho[x]/n
        rhi += L[x][x]*rho[x]/n
        for i in range(n):
            y = x^(1 << i)
            if x<y:
                lo, hi = sqrt_bounds(rho[x]*rho[y])
                rlo += 2*L[x][y]*lo/n
                rhi += 2*L[x][y]*hi/n
                left = bits[(i-1)%n]
                Dlo[left] += lo/n
                Dhi[left] += hi/n
    assert rhi < lam_lo
    assert Dhi[1] < ablo and Dhi[0] < bdlo
    assert rlo == Dlo[1]+Dlo[0]/2-diag
    assert rhi == Dhi[1]+Dhi[0]/2-diag

    result = {
        'status':'PASS',
        'scope':'Exact finite N=3 discrepancy and local pair-relaxation compatibility; no infinite-volume spectral conclusion',
        'runtime':{'python':sys.version,'implementation':platform.python_implementation(),
                   'dependencies':'Python standard library only'},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_url':'https://doi.org/10.1140/epjb/s10051-022-00281-5',
        'parameters':{'c':'1/2','kappa':'1','epsilon':'1/2','e_minus_s':'1/2',
                      'Monthus_k':'-s','project_time_factor':'2','N':3},
        'states_little_endian_x1_x2_x3':[list(s) for s in states],
        'matrix_4L':M,
        'cw_N3':{'positive_integer_h':h,'shift_used_only_for_discovery':12,
                 'iterations_used_only_for_discovery':250,
                 'ratios_lambda_per_site':[rat(q) for q in ratios],
                 'lambda_per_site':interval_record(lam_lo,lam_hi),
                 'matches_user_coarse_interval':True},
        'trial':{'a':rat(a),'b':rat(b),'d':rat(d),'normalization':rat(a+2*b+d),
                 'sqrt_ab':interval_record(ablo,abhi,45),
                 'sqrt_bd':interval_record(bdlo,bdhi,45),
                 'F':interval_record(F_lo,F_hi),
                 'F_minus_exact_spectrum_lower_bound':rat(gap),
                 'gap_decimal_lower':outward(gap,30)},
        'local_maximum':{'matrix_8_D_inverse_K':A8,'positive_integer_v':v,
                         'ratios':[rat(q) for q in local_ratios],
                         'maximum':interval_record(local_lo,local_hi),
                         'point_discovered_by_integer_power_iteration_and_rounding':True},
        'global_realization_N3':{'rho_by_state':[rat(q) for q in rho],
                                 'all_3_edge_marginals':pair_results,
                                 'full_support':True,
                                 'D_plus':interval_record(Dlo[1],Dhi[1]),
                                 'D_minus':interval_record(Dlo[0],Dhi[0]),
                                 'rayleigh_per_site':interval_record(rlo,rhi)},
        'stationary_spatial_markov_realization':{'state_order':['+','-'],
                'pi':[rat(x) for x in pi],'transition':[[rat(x) for x in row] for row in P]},
        'equality_obstruction':{'t_plus':rat(tplus),'t_minus':rat(tminus),
                'ad_minus_b_squared':rat(a*d-b*b),
                'ratio_path_000_100_110':rat(path12),
                'ratio_path_000_010_110':rat(path21),
                'same_positive_global_density_cannot_satisfy_both_paths':True},
        'limitations':[
                'Not formal verification or external peer review',
                'No extrapolation from finite N=3 to the thermodynamic SCGF',
                'The pair itself is globally realizable; only the simultaneous equality conditions fail']
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.output.exists():
        raise FileExistsError('Refusing to overwrite existing independent result')
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','parameters','cw_N3','trial','local_maximum','equality_obstruction']},indent=2))
    print('Saved:',args.output)


if __name__=='__main__':
    if not __debug__: raise RuntimeError('Do not disable certificate assertions')
    main()
