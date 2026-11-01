"""A posteriori log-Perron accuracy from a positive power, without SciPy.

All supplied decimals denote exact rationals. This conservative verifier is
intended for small reference volumes; it does not form a dense matrix.
"""
from fractions import Fraction as F
from mpmath import iv
from verify_certificate import endpoint, outward


def fi(x):
    x = F(x)
    return iv.mpf(x.numerator)/x.denominator


def upper(x):
    return endpoint(x._mpi_[1])


def lower(x):
    return endpoint(x._mpi_[0])


def log_error_bound(N, c, kappa, s, candidate, dps=90, places=65):
    """Return a rigorous upper bound on inf_b ||log(candidate)-log(h)-b||inf."""
    if not isinstance(places, int) or places < 1:
        raise ValueError('Positive output places required')
    c, kappa, s = map(F, (c, kappa, s))
    h = list(map(F, candidate))
    if N < 3 or len(h) != 2**N or min(h) <= 0 or not (0 < c < 1 and kappa > 0):
        raise ValueError("Invalid model or positive reference vector")
    iv.dps = dps
    J = max(c, 1-c)*(1+kappa)
    shift = 2*N*J
    es = iv.exp(-fi(s))
    rows = []
    for x in range(2**N):
        rates = [(kappa+((x >> ((i-1) % N)) & 1)) *
                 (1-c if ((x >> i) & 1) else c) for i in range(N)]
        rows.append((fi(1-sum(rates)/shift), [fi(r/shift)*es for r in rates]))
    qh = [fi(x) for x in h]
    for _ in range(N):
        qh = [rows[x][0]*qh[x]+sum(rows[x][1][i]*qh[x ^ (1 << i)]
                                  for i in range(N)) for x in range(2**N)]
    ratios = [q/fi(x) for q, x in zip(qh, h)]
    ratio_lo = min(lower(x) for x in ratios)
    ratio_hi = max(upper(x) for x in ratios)
    if ratio_lo <= 0:
        raise ArithmeticError("Positive lower endpoint lost; increase precision")
    residual_hi = upper(iv.log(fi(ratio_hi)/fi(ratio_lo)))
    beta_lo = min(F(1, 2), lower(fi(kappa*min(c, 1-c)/shift)*iv.exp(-fi(abs(s)))))
    B_hi = upper(1+(iv.exp(fi(abs(s)))-1)/2)
    chi_hi = (B_hi/beta_lo)**N
    # 1/(1-tau) <= (chi+1)/2, with no subtractive cancellation.
    eta_hi = residual_hi*(chi_hi+1)/4
    return dict(N=N, interval_dps=dps,
                projective_residual_upper=outward(residual_hi, True, 75),
                log_error_upper=outward(eta_hi, True, places),
                log_error_upper_exact=str(eta_hi), output_places=places,
                contraction_gap_lower=outward(2/(chi_hi+1), False, places),
                contraction_gap_lower_exact=str(2/(chi_hi+1)),
                scope="Certified up to one additive constant in log h; finite N only.")
