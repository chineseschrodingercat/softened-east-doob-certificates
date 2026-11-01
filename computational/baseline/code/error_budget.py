"""Symbolic size planning and interval checks: never allocate 2**N states."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from mpmath import mp, iv
from certified_reference import fi, upper, lower


def approx(x):
    return mp.mpf(x.numerator)/x.denominator


def plan(c="0.3", kappa="0.01", s0="1e-3150", block_length=33554432, tolerance="1e-6"):
    mp.dps = 100
    iv.dps = 100
    s0_input = str(s0)
    c, kappa, s0, tol = map(F, (c, kappa, s0, tolerance))
    if not (0<c<1 and kappa>0 and s0>=0 and tol>0 and block_length>=1):
        raise ValueError("Invalid planning parameters")
    ell = int(block_length)
    a = F(1, 2*ell)
    J = max(c, 1-c)*(1+kappa)
    A = fi(J)*iv.exp(fi(s0+1))
    theta = iv.log(1+fi(tol)/(4*A))
    theta_lo = lower(theta)
    Hstar = 8*(1+2/fi(a))/(3*(1-iv.exp(-fi(a)))**2)
    R = max(2, int(mp.ceil(approx(upper(2/fi(a)*iv.log(4*Hstar/fi(theta_lo))))))+2)
    N_periodic = 16*ell*iv.log(32*iv.exp(2)*R/fi(theta_lo))
    N_exterior = R-1+iv.log(128*R/(3*fi(theta_lo)*(1-iv.exp(-fi(a)))))/fi(a)
    N = max(2*R-1, int(mp.ceil(approx(upper(N_periodic))))+2,
            int(mp.ceil(approx(upper(N_exterior))))+2)
    eta, qv = theta_lo/(16*R), theta_lo/(8*R)
    H = 8*iv.exp(-fi(a)*R)*((R+1)-R*iv.exp(-fi(a)))/(3*(1-iv.exp(-fi(a)))**2)
    delta = 2*iv.exp(2)*iv.exp(-fi(F(N,16*ell))) + 8*iv.exp(-fi(a)*(N-R+1))/(3*(1-iv.exp(-fi(a))))
    terms = [H, 2*R*delta, fi(4*R*eta), fi(2*R*qv)]
    if not all(upper(x) <= theta_lo/4 for x in terms):
        raise ArithmeticError("Planning inequalities did not verify")
    E = sum(terms)
    ideal_width = 2*A*(iv.exp(E)-1)
    arith = tol/2
    total = upper(ideal_width)+arith
    if total > tol:
        raise ArithmeticError("Total budget failed")
    return dict(status="INTERVAL_BUDGET_VERIFIED", c=str(c), kappa=str(kappa),
        s0_input=s0_input, block_length=ell, tolerance=str(tol), R=R,N=N,
        theta_lower=mp.nstr(approx(theta_lo),22), eta_max=mp.nstr(approx(eta),22),
        density_rounding_max=mp.nstr(approx(qv),22),
        log10_reference_dimension=mp.nstr(N*mp.log10(2),22),
        log10_graph_vertices=mp.nstr((2*R-2)*mp.log10(2),22),
        four_log_error_terms_upper=[mp.nstr(approx(upper(x)),22) for x in terms],
        width_bound_upper_display=mp.nstr(approx(total),22),
        check_precision_digits=100,
        scope="Verifies the size/error algebra, conditional on the supplied A1 block length and tilt window. No state space allocated.")


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/budget_plan.json"))
    args=parser.parse_args()
    result=plan()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
