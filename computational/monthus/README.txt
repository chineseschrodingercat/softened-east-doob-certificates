PUBLISHED TWO-SPIN CLOSURE COUNTEREXAMPLE (SM11)
Reference: C. Monthus, European Physical Journal B 95, article 32 (2022).
DOI: 10.1140/epjb/s10051-022-00281-5. Formula numbering is from the version of record.
At c=1/2, kappa=1, exp(-s)=1/2 in the article's units, the Monthus softening
is epsilon=1/2, her tilt is k=-s, and her generator is multiplied by 2.
The exact-SCGF identification from published equations (109)-(111) is ruled out
at this point. This is not an assessment of every result or the hard limit.

From the package root (Python standard library only):
  python -B monthus/code/replay_monthus.py --output replay_primary
  python -B monthus/code/replay_independent.py --output replay_independent

Outputs are fresh subdirectories under monthus/. Both commands verify the full
package manifest. Do not use python -O. No NumPy, SciPy, optimization or network
is needed for these two commands. Primary replay recomputes the frozen expected
record exactly. Independent replay uses separately written matrix/word code.

The simple feasible pair (.110,.187,.516) satisfies F>-0.3498 and is globally
realizable. The 32-edge R3 witness proves lambda_n/n<-0.3499 for every n>=5.
The 128-edge R4 witness gives [-0.350871819626929,-0.350726876487370] for n>=7.
All decimals in the factor and potential tables are exact rational numbers.
Full n=5 and n=7 checks cross-check implementation; the all-n guarantee follows
from the cyclic telescoping identity proved in SM11, not finite-size extrapolation.
The finite n=3 matrix enclosure is independently reconstructed as well.

results/published_results.json: exact expected primary output and rational endpoints.
results/uniform_R3_simple_candidate.json: 8 factors and 16 upper/lower potentials.
results/uniform_R4_candidate.json: 16 factors and 64 upper/lower potentials.
results/R3_UPPER_EDGES.csv: all 32 simple upper inequalities, with exact slacks.

These are posterior certificates. They do not enlarge the prior locality window.
Publisher PDF, archived drafts, optimizer output and internal review reports are
not included. The main Supplementary Information supplies the mathematical proof.
