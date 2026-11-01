"""Exact-decimal, outward-rounded h=1 certificate, all periodic n >= 3."""
import hashlib
import json
import time
from fractions import Fraction as F
from mpmath import iv
from verify_certificate import endpoint, outward


def constant_certificate(c="0.3", kappa="0.01", s="0.001", dps=80):
    started = time.perf_counter()
    cf, kf, sf = map(F, (c, kappa, s))
    if not (0 < cf < 1 and kf > 0):
        raise ValueError("Require 0<c<1 and kappa>0")
    # The only simple cycles in the two-state directed activity graph.
    means = [kf*cf, (kf+1)*(1-cf), (kf+cf)/2]
    iv.dps = dps
    exponential = iv.exp(-iv.mpf(sf.numerator)/sf.denominator)
    a, b = endpoint(exponential._mpi_[0])-1, endpoint(exponential._mpi_[1])-1
    lower = outward(min(a*x for x in means), False)
    upper = outward(max(b*x for x in means), True)
    inputs = dict(model="periodic_soft_east_activity", c=c, kappa=kappa, s=s,
                  method="constant", trial="h=1", all_integer_n_at_least=3)
    return dict(**inputs, status="VERIFIED_CONSTANT_CYCLES", interval_dps=dps,
                cycle_activity_means=[str(x) for x in means],
                lower_density=lower, upper_density=upper,
                exact_width=str(F(upper)-F(lower)), width=float(F(upper)-F(lower)),
                input_sha256=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest(),
                seconds=time.perf_counter()-started,
                proof="Every closed binary walk decomposes into loops 0,1 and cycle 01; Collatz-Wielandt applies to h=1.")


if __name__ == "__main__":
    print(json.dumps(constant_certificate(), indent=2))
