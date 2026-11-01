"""Interval verification of the practical-window condition and log residual.

Stored decimal candidates and parameters are exact rationals. Polynomial
transforms use rational arithmetic for the trial and interval arithmetic for
the residual. Exponentials and every weighted norm are outward enclosed.
The mathematical interpretation is stated in PRACTICAL_LOCALITY_PROOF.md.
"""
import argparse
import hashlib
import json
import time
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv
from verify_certificate import endpoint, outward


def I(x):
    x = F(x)
    return iv.mpf(x.numerator) / x.denominator


def upper(x):
    return endpoint(x._mpi_[1])


def lower(x):
    return endpoint(x._mpi_[0])


def record(x):
    lo, hi = lower(x), upper(x)
    return dict(lower=outward(lo, False, 65), upper=outward(hi, True, 65),
                lower_exact=str(lo), upper_exact=str(hi))


def constants(c='0.3', kappa='1', s0='0.02', a='0.02', radius='0.054', dps=80):
    iv.dps = dps
    c, kappa, s0, a, radius = map(F, (c, kappa, s0, a, radius))
    assert 0 < c < 1 and kappa > 0 and s0 > 0 and a > 0 and radius > 0
    q, w, b = max(c, 1-c), c*(1-c), 1-2*c
    gamma, D = kappa+c, 2+abs(b)/q
    ea, er = iv.exp(I(a)), iv.exp(I(D*radius))
    J = I(q*gamma)+I(q*q)*ea
    A = I(abs(q*(2*w+gamma*b)))+I(abs(b)*q*q)*ea
    E = er-1-I(D*radius)
    ell = (I(q)*ea+J*I(D)*(er-1))/I(gamma)
    ell0 = (I(q)+I(q*gamma+q*q)*I(D)*(er-1))/I(gamma)
    margin = (I(gamma)-I(q)*ea)*I(radius)-J*E-(iv.exp(I(s0))-1)*A
    # For c=.3, A and J are affine in gamma; this slope proves the entire
    # half-line kappa >= 1, in addition to the monotonicity of ell.
    slope = I(radius)-I(q)*E-(iv.exp(I(s0))-1)*I(q*abs(b))
    passed = upper(ell) < 1 and lower(margin) > 0
    result = dict(status='CERTIFIED_CONTRACTION' if passed else 'INCONCLUSIVE',
                  c=str(c), kappa=str(kappa), s0=str(s0), a=str(a), radius=str(radius),
                  gamma=str(gamma), q=str(q), D=str(D), interval_dps=dps,
                  ell=record(ell), ell0=record(ell0), self_map_margin=record(margin),
                  kappa_slope=record(slope), J=record(J), A=record(A),
                  half_line_c03_kappa_ge1=(c == F('0.3') and kappa == 1 and passed
                                         and lower(slope) > 0))
    return result


def transform(values, c, q):
    """Tensor basis transform: f0,f1 -> ((1-c)f0+c f1, q(f1-f0))."""
    out = list(values)
    step = 1
    while step < len(out):
        for start in range(0, len(out), 2*step):
            for j in range(start, start+step):
                v0, v1 = out[j], out[j+step]
                out[j] = (1-c)*v0+c*v1
                out[j+step] = q*(v1-v0)
        step *= 2
    return out


def collect(coef):
    out = {}
    for mask in range(1, len(coef)):
        shift = (mask & -mask).bit_length()-1
        key = mask >> shift
        out[key] = out.get(key, 0)+coef[mask]
    return out


def verify_candidate(path, condition=None, dps=80, local_a='0.3', local_cap='0.02'):
    started = time.perf_counter()
    iv.dps = dps
    path = Path(path)
    raw = path.read_bytes()
    data = json.loads(raw)
    assert data['model'] == 'periodic_soft_east_activity'
    assert data['density_convention'] == 'U=sum(density); no implicit s factor'
    assert data['bit_order'] == 'leftmost_is_least_significant'
    c, kappa, s = (F(data[k]) for k in ('c', 'kappa', 's'))
    condition = condition or constants(str(c), str(kappa), dps=dps)
    assert condition['status'] == 'CERTIFIED_CONTRACTION'
    assert c == F(condition['c']) and kappa == F(condition['kappa'])
    assert abs(s) <= F(condition['s0'])
    a, radius = F(local_a), F(local_cap)
    assert a >= F(condition['a']) and radius <= F(condition['radius'])
    q, gamma = max(c, 1-c), kappa+c
    R, v = data['R'], list(map(F, data['density']))
    assert R >= 2 and len(v) == 1 << R
    coeff = collect(transform(v, c, q))
    weights = [iv.exp(I(a*d)) for d in range(2*R-1)]
    trial_norm = sum((I(abs(z)*key.bit_count())*weights[key.bit_length()-1]
                      for key, z in coeff.items()), I(0))
    if upper(trial_norm) > radius:
        return dict(status='OUTSIDE_CERTIFIED_BALL', trial_norm=record(trial_norm))
    es = iv.exp(I(s))
    residual = []
    for word in range(1 << (2*R-1)):
        delta = F(0)
        for start in range(R):
            w = (word >> start) & ((1 << R)-1)
            delta += v[w ^ (1 << (R-1-start))]-v[w]
        bit, prev = (word >> (R-1)) & 1, (word >> (R-2)) & 1
        rate = (kappa+prev)*(c if bit == 0 else 1-c)
        residual.append(I(rate)*(iv.exp(I(delta))-es))
    residual_coeff = collect(transform(residual, I(c), I(q)))
    norm0 = sum((abs(z) for z in residual_coeff.values()), I(0))
    norma = sum((abs(z)*weights[key.bit_length()-1]
                 for key, z in residual_coeff.items()), I(0))
    # These are derivative bounds on the total-norm cap, not a claim that
    # the zero-centered cap is invariant. Invariance is checked about v.
    D = 2+abs(1-2*c)/q
    expcap = iv.exp(I(D*radius))
    J = I(q*gamma)+I(q*q)*iv.exp(I(a))
    local_ella = (I(q)*iv.exp(I(a))+J*I(D)*(expcap-1))/I(gamma)
    local_ell0 = (I(q)+I(q*gamma+q*q)*I(D)*(expcap-1))/I(gamma)
    ella, ell0 = upper(local_ella), upper(local_ell0)
    if ella >= 1:
        return dict(status='LOCAL_DERIVATIVE_INCONCLUSIVE', ell=record(local_ella))
    E0 = I(upper(norm0)/(gamma*(1-ell0)))
    Ea = I(upper(norma)/(gamma*(1-ella)))
    local_margin = I(radius)-trial_norm-Ea
    if lower(local_margin) <= 0:
        return dict(status='LOCAL_BALL_INCONCLUSIVE', margin=record(local_margin))
    flip = E0/I(q)
    # For a target reference density of length R, N >= 2R-1 suffices:
    # the trial flip sees no other copy of the zero-extended target word.
    reference_error = 2*flip
    period_error = 2*flip
    graph_width = 2*iv.exp(-I(s))*norm0
    return dict(status='CERTIFIED_INFINITE_AND_ALL_RING_LOG_ERROR',
                candidate=path.name, candidate_sha256=hashlib.sha256(raw).hexdigest(),
                c=str(c), kappa=str(kappa), s=str(s), R=R, n_min=2*R-1,
                interval_dps=dps, condition=condition,
                local_a=str(a), local_cap=str(radius), local_ella=record(local_ella),
                local_ell0=record(local_ell0), local_ball_margin=record(local_margin),
                trial_Xa_norm=record(trial_norm), residual_F0=record(norm0),
                residual_Fa=record(norma), infinite_X0_error=record(E0),
                infinite_Xa_error=record(Ea), all_ring_X0_error=record(E0),
                flip_log_error=record(flip), periodic_flip_error=record(period_error),
                exact_vacuum_reference_coefficient_error=record(reference_error),
                residual_only_graph_width=record(graph_width), seconds=time.perf_counter()-started,
                scope=('Analytic contraction plus all local residual coefficients; '
                       'certifies an infinite-chain interaction and every n>=n_min. '
                       'Reference coefficient error compares exact finite Perron and '
                       'infinite vacuum densities of length R for N>=n_min. '
                       'The stored trial need not be an exact eigenvector.'))


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('candidate', type=Path)
    p.add_argument('--dps',type=int,default=80)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=verify_candidate(args.candidate,dps=args.dps)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','candidate','R','seconds'] if k in result},indent=2))
