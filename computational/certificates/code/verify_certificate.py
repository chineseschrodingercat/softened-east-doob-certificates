"""Independent all-size graph verifier: stdlib plus mpmath only.

Decimal inputs and graph potentials are exact rationals. Only exp is
enclosed by mpmath.iv. All subsequent endpoint operations are exact rational
operations; the delivered decimal endpoints are rounded outwards.
"""
import argparse
import hashlib
import json
import time
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv


def endpoint(t):
    sign, mantissa, exponent, bits = t
    if bits < 0:
        raise ArithmeticError("Nonfinite interval endpoint")
    result = F((-1 if sign else 1)*mantissa)
    return result*(2**exponent) if exponent >= 0 else result/F(2**(-exponent))


def outward(value, upper, places=35):
    scale = 10**places
    x = value*scale
    q = -((-x.numerator)//x.denominator) if upper else x.numerator//x.denominator
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q//scale}.{q%scale:0{places}d}"


def verify(path, dps=60, write_result=True, places=35):
    if not isinstance(places, int) or places < 1 or not isinstance(dps, int) or dps < 15:
        raise ValueError('Positive endpoint places and interval dps >= 15 are required')
    started = time.perf_counter()
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    assert data["model"] == "periodic_soft_east_activity"
    assert data["density_convention"] == "U=sum(density); no implicit s factor"
    assert data["bit_order"] == "leftmost_is_least_significant"
    R = data["R"]
    assert R >= 2 and data["n_min"] == 2*R-1
    c, kappa, s = [F(data[key]) for key in ["c", "kappa", "s"]]
    assert 0 < c < 1 and kappa > 0
    v, zl, zu = [[F(x) for x in data[key]] for key in ["density", "z_lower", "z_upper"]]
    V = 1 << (2*R-2)
    assert len(v) == 1 << R and len(zl) == len(zu) == V
    iv.dps = dps
    lo, hi = None, None
    max_exp_radius = F(0)
    for word in range(2*V):
        delta = F(0)
        for start in range(R):
            window = (word >> start) & ((1 << R)-1)
            delta += v[window ^ (1 << (R-1-start))]-v[window]
        bit, prev = (word >> (R-1)) & 1, (word >> (R-2)) & 1
        rate = (kappa+prev)*(c if bit == 0 else 1-c)
        exponent = delta-s
        exp_iv = iv.exp(iv.mpf(exponent.numerator)/exponent.denominator)
        a = rate*(endpoint(exp_iv._mpi_[0])-1)
        b = rate*(endpoint(exp_iv._mpi_[1])-1)
        max_exp_radius = max(max_exp_radius, (b-a)/2)
        tail, head = word & (V-1), word >> 1
        a += zl[tail]-zl[head]
        b += zu[tail]-zu[head]
        lo = a if lo is None else min(lo, a)
        hi = b if hi is None else max(hi, b)
    lower, upper = outward(lo, False, places), outward(hi, True, places)
    assert F(lower) <= lo <= hi <= F(upper)
    result = dict(status="VERIFIED_ALL_EDGES", candidate_sha256=hashlib.sha256(raw).hexdigest(),
        c=data["c"], kappa=data["kappa"], s=data["s"], R=R,
        all_integer_n_at_least=2*R-1, lower_density=lower, upper_density=upper,
        width=float(F(upper)-F(lower)), exact_width=str(F(upper)-F(lower)),
        endpoint_places=places, vertices=V, edges=2*V,
        max_edge_enclosure_radius=float(max_exp_radius),
        max_edge_enclosure_radius_exact=str(max_exp_radius), interval_dps=dps,
        seconds=time.perf_counter()-started,
        scope="Exact decimal candidate and duals; outward exponential enclosures. No Perron accuracy or optimizer optimality assumption.")
    if write_result:
        target = Path(path).with_name(Path(path).stem+"_verified.json")
        target.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--places", type=int, default=35)
    args = parser.parse_args()
    print(json.dumps(verify(args.candidate, args.dps, places=args.places), indent=2))
