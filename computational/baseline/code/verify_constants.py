"""Direct interval recheck, independent of the v7 log-inequality implementation."""
from pathlib import Path
import json
from mpmath import iv
from mpmath.libmp import mpf_cmp

ROOT=Path(__file__).resolve().parents[1]
iv.dps=100


def leq(a,b):
    # Compare binary endpoints without expanding huge exponents into
    # fractions with multi-gigabit denominators.
    return mpf_cmp(a._mpi_[1],b._mpi_[0])<=0


def main():
    old=json.loads((ROOT/"inputs"/"parameter_verification.json").read_text())
    cases=[dict(c=row["c"],kappa=row["kappa"],T=row["T"],l=row["l"],
                beta=row["beta"],s0=row["s0"],source="v7_saved_parameters")
           for row in old["cases"]]
    for c,kappa,T in [("0.01","0.000001",2**48),("0.99","0.000001",2**48),
                       ("0.3","0.00000001",2**52),("0.001","0.01",2**36)]:
        cases.append(dict(c=c,kappa=kappa,T=T,l=32*T,beta="1e-3199",s0="1e-3200",
                          source="new_extreme_parameter_case"))
    verified=[]
    for row in cases:
        c,kappa=iv.mpf(row["c"]),iv.mpf(row["kappa"])
        # c values are exact input decimals; branch selection uses their
        # strings as rational decimals and is not a numerical interval test.
        from fractions import Fraction
        m=c if Fraction(row["c"])<=Fraction(1,2) else 1-c
        q=1-m
        alpha=1/(1+kappa) if Fraction(row["kappa"])<=1 else iv.mpf("0.5")
        g=2*(1-alpha)
        T,l=row["T"],row["l"]
        tau=iv.ln(2*l-1)/(4*g*c*(1-c))
        t=T-tau
        eps=iv.mpf("1e-32")
        C0=2/m**2
        assert leq(iv.mpf(1),t)
        clock=(2*iv.exp(1)*T/l)**l
        assert leq(clock,(eps/(2*iv.exp(1)))**iv.mpf("2.5"))
        prefactor=128*(1+t)/(1-alpha)**2
        target=(eps/(iv.exp(1)*C0))**126
        assert leq(prefactor*alpha**(l//2),target)
        assert leq(prefactor*iv.exp(-2*(1-iv.sqrt(alpha))*t),target)
        beta=iv.mpf(row["beta"])
        s0=iv.mpf(row["s0"])
        assert leq(4*l*T*beta,iv.exp(-1)*(eps/2)**97)
        n0=16*l
        assert leq(n0*beta*(1+iv.mpf(n0)*iv.ln(1/m)/2)/g,iv.mpf("0.25"))
        assert leq(2*q*s0*iv.exp(s0),beta)
        ratio=128*iv.exp(2)*iv.mpf("0.0001")
        assert leq(ratio/(1-ratio),iv.mpf(1)/3)
        verified.append(dict(**row,status="DIRECT_INTERVAL_VERIFIED"))
    # Excluding only the high-energy half-line is insufficient for the
    # resolvent norm claim. This is a scalar domain diagnostic, not an
    # East-model or A1 counterexample.
    alpha=iv.mpf("0.5")
    delta=iv.mpf("0.1")
    actual=2*alpha/delta**2
    overbroad_bound=iv.sqrt(alpha)/delta
    assert leq(2*overbroad_bound,actual)
    result=dict(status="ALL_DIRECT_INTERVAL_CHECKS_PASSED",
        precision=iv.dps,parameter_cases=len(verified),cases=verified,
        scalar_domain_diagnostic=dict(HB="2",R="-1",z="1.9",actual_norm="100",
                                     invalid_bound="sqrt(0.5)/0.1",
                                     scope="Rejects an overly broad domain interpretation; not a counterexample to A1."),
        scope="Direct inequalities without the v7 logarithmic comparison code.")
    (ROOT/"results"/"independent_constants.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(f"Verified {len(verified)} direct interval parameter cases at {iv.dps} digits.")


if __name__=="__main__":
    main()
