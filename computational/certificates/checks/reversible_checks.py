"""Independent mathematical-review diagnostics; no output files are written.

Run: python -B certificates/checks/reversible_checks.py
The random finite-state experiment is floating diagnostic evidence, not a
proof. The closed-form transfer comparison uses interval arithmetic but is
only a low-dimensional implementation crosscheck of the general validator.
"""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh
from mpmath import iv

EXTENSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXTENSION_ROOT / "code"))
from certify_infinite import verify, seed
from phase0 import I, bounds


def random_reversible_checks():
    rng = np.random.default_rng(928341)
    records = []
    for n in [2, 3, 4, 5]:
        count = 1 << n
        states = np.arange(count)
        bits = (states[:, None] >> np.arange(n)) & 1
        for rep in range(12):
            fields = rng.normal(0, .5, n)
            couplings = rng.normal(0, .03, (n, n))
            couplings = np.triu(couplings, 1)
            U = bits @ fields + np.einsum(
                "si,ij,sj->s", bits, couplings, bits
            )
            nu = np.exp(U - U.max())
            nu /= nu.sum()
            p = np.empty((count, n))
            b = np.empty((count, n))
            for i in range(n):
                zero = states & ~(1 << i)
                one = zero | (1 << i)
                p[:, i] = nu[one] / (nu[zero] + nu[one])
                k = np.exp(.1 * rng.normal(size=count))
                k = k[zero]
                b[:, i] = k * np.where(
                    bits[:, i], 1 - p[:, i], p[:, i]
                )
            C = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    if i != j:
                        C[i, j] = np.max(np.abs(
                            p[:, i] - p[states ^ (1 << j), i]
                        ))
            alpha = max(np.max(C.sum(0)), np.max(C.sum(1)))
            beta = min(
                np.min(b[:, i] + b[states ^ (1 << i), i])
                for i in range(n)
            )
            V = (
                .15 * bits @ rng.normal(size=n)
                + .07 * np.einsum(
                    "si,ij,sj->s", bits,
                    rng.normal(size=(n, n)), bits
                )
            )
            d = max(
                np.max(np.abs(V - V[states ^ (1 << i)]))
                for i in range(n)
            )
            A = np.diag(V - b.sum(axis=1))
            for i in range(n):
                flip = states ^ (1 << i)
                A[states, flip] = np.sqrt(b[:, i] * b[flip, i])
            actual = (
                eigvalsh(
                    A, subset_by_index=[count - 1, count - 1]
                )[0] - nu @ V
            ) / n
            bound = d*d / (4*beta*(1-alpha)**2)
            assert alpha < 1
            assert actual >= -1e-12
            assert actual <= bound + 1e-12
            records.append((n, alpha, actual / bound))
    return {
        "random_seed": 928341,
        "random_nontranslation_invariant_cases": len(records),
        "max_ratio_actual_gap_over_bound": max(r[2] for r in records),
        "max_alpha": max(r[1] for r in records),
        "status": "PASS_DIAGNOSTIC_ONLY"
    }


def closed_form_transfer_checks():
    iv.dps = 100
    results = []
    cases = [
        ("0.04", ["0.01", "-0.02", "0.05", "0.08"]),
        ("-0.03", ["-0.04", "0.03", "0.06", "-0.01"]),
        ("0", ["0.03", "-0.05", "-0.02", "0.04"]),
    ]
    for tilt, v in cases:
        data = {
            "R": 2, "c": "0.3", "kappa": "0.17",
            "s": tilt, "density": v
        }
        c = I(F(3, 10))
        k = I(F(17, 100))
        s = I(F(tilt))
        val = [I(F(x)) for x in v]
        tm = [[None]*2 for _ in range(2)]
        for a in range(2):
            for b in range(2):
                tm[a][b] = (
                    (c if b else 1-c) * iv.exp(2*val[a+2*b])
                )
        a, b = tm[0]
        cc, d = tm[1]
        rho = (a+d+iv.sqrt((a-d)**2+4*b*cc))/2
        r = [b, rho-a]
        l = [cc, rho-a]
        den = l[0]*r[0]+l[1]*r[1]
        marginal = [l[j]*r[j]/den for j in range(2)]
        P = [
            [tm[i][j]*r[j]/(rho*r[i]) for j in range(2)]
            for i in range(2)
        ]
        exact = I(0)
        for prev in range(2):
            for bit in range(2):
                for nxt in range(2):
                    delta = (
                        val[prev+2*(1-bit)] + val[(1-bit)+2*nxt]
                        - val[prev+2*bit] - val[bit+2*nxt]
                    )
                    phi = (
                        (k+prev)*(c if bit == 0 else 1-c)
                        * (iv.exp(delta-s)-1)
                    )
                    exact += (
                        marginal[prev]*P[prev][bit]*P[bit][nxt]*phi
                    )
        lo, hi = bounds(exact)
        precise = verify(data, seed(data), 80)["Q_infinity"]
        coarse = verify(
            data, {"right": ["1", "1"], "left": ["1", "1"]}, 80
        )["Q_infinity"]
        assert F(precise["lower"]) <= hi
        assert F(precise["upper"]) >= lo
        assert F(coarse["lower"]) <= lo
        assert F(coarse["upper"]) >= hi
        results.append({
            "s": tilt,
            "density": v,
            "exact_mid": float((lo+hi)/2),
            "precise_width": float(
                F(precise["upper"])-F(precise["lower"])
            ),
            "coarse_seed_width": float(
                F(coarse["upper"])-F(coarse["lower"])
            )
        })
    return {
        "independent_closed_form_2x2_checks": results,
        "status": "PASS"
    }


if __name__ == "__main__":
    print(json.dumps({
        "random_reversible_checks": random_reversible_checks(),
        "closed_form_transfer_checks": closed_form_transfer_checks()
    }, indent=2))
