"""Independent R=2 analytic Perron checks for the Q_infinity validator.

Run with python -B certificates/checks/closed_form_checks.py. This script writes no
files. The analytic reference uses a direct two-state matrix, explicit Perron
data, and a direct three-spin residual; it does not call edges, apply, or
local_tables from the implementation under review.
"""
import sys
import json
from pathlib import Path
from fractions import Fraction as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from mpmath import iv
from certify_infinite import verify, seed
from phase0 import I, bounds


def independent_r2(data, dps):
    iv.dps = dps
    c, k, s = [I(F(data[q])) for q in ['c', 'kappa', 's']]
    v = [I(F(x)) for x in data['density']]
    # Construct by spin pairs, independently of the implementation's shifts.
    A = [[(c if j else 1-c)*iv.exp(2*v[i+2*j]) for j in (0, 1)]
         for i in (0, 1)]
    a, b = A[0]
    cc, d = A[1]
    rho = (a+d+iv.sqrt((a-d)**2+4*b*cc))/2
    right = [b, rho-a]
    left = [cc, rho-a]
    num = I(0)
    den = I(0)
    for x in (0, 1):
        for y in (0, 1):
            for z in (0, 1):
                weight = left[x]*A[x][y]*A[y][z]*right[z]
                delta = (v[x+2*(1-y)]+v[(1-y)+2*z]
                         -v[x+2*y]-v[y+2*z])
                rate = (k+x)*(c if y == 0 else 1-c)
                phi = rate*(iv.exp(delta-s)-1)
                den += weight
                num += weight*phi
    return bounds(num/den), bounds(rho)


def main():
    cases = [
        {'R': 2, 'c': '0.3', 'kappa': '0.1', 's': '-0.03',
         'density': ['0.06', '-0.02', '0.09', '-0.03']},
        {'R': 2, 'c': '0.7', 'kappa': '0.01', 's': '0.04',
         'density': ['-0.11', '0.02', '0.03', '0.08']},
        {'R': 2, 'c': '0.17', 'kappa': '1', 's': '0',
         'density': ['0.3', '-0.24', '0.16', '0.11']},
        {'R': 2, 'c': '0.9', 'kappa': '0.3', 's': '-0.2',
         'density': ['-0.9', '-0.12', '0.17', '0.19']},
    ]
    results = []
    for data in cases:
        actual, rho = independent_r2(data, 100)
        for mode in ('accurate', 'poor'):
            trial = (seed(data) if mode == 'accurate'
                     else {'right': ['1', '2.5'], 'left': ['1', '0.4']})
            cert = verify(data, trial, 80)
            q = cert['Q_infinity']
            r = cert['perron_root']
            overlap = max(F(q['lower']), actual[0]) <= min(F(q['upper']), actual[1])
            root_overlap = max(F(r['lower']), rho[0]) <= min(F(r['upper']), rho[1])
            assert overlap and root_overlap
            results.append({
                'case': data,
                'mode': mode,
                'Q_center_approx': float(sum(actual, F(0))/2),
                'Q_width': float(F(q['upper'])-F(q['lower'])),
                'rho_overlap': root_overlap,
                'Hilbert_right_bound': cert['right_Hilbert_error_upper_float'],
            })
    print(json.dumps({'status': 'PASS_INDEPENDENT_ANALYTIC_R2_CHECKS',
                      'checks': len(results), 'results': results}, indent=2))


if __name__ == '__main__':
    main()
