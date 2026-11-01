# External comparison: exact endpoints, frozen inputs, and generation costs

This directory compares the same periodic softened-East activity generator at
`c=0.3` and `(kappa,s)=(1,0.02),(0.1,0.01),(0.01,0.001)`. It separately records
the finite target `lambda_18/18` and the thermodynamic limit. The operator and
source-to-model correspondences are in `THEORY_AUDIT.md` and article SM12.

## Verify without optimization

From the package root, after installing `requirements.txt`:

```
python -B external/code/replay_external.py --output replay_external
```

Use a new output name each time. The default checks the entire package manifest,
then reconstructs every retained certificate and compares exact decimal endpoints.
No SDP solver, optimizer, author runtime, network, reference-vector cache, or GPU
is needed to replay. Do not use Python `-O`. `--skip-manifest` is a development
option and was not used for the final fresh-package replay.

For each successful result:

- `candidate.json` supplies rational Anderson vectors or SDP multipliers, maps,
  trace bounds and factor-residual PSD witnesses.
- Block/variance trials instead supply `trial.json`, `trial_graph.json`,
  `seed.json` and separately certified finite and infinite centers.
- `result.json` contains exact outward endpoints, generation and verification
  times, memory, solver/optimizer status, and coefficient-error corrections.
- `failure.json` and `process.log` retain unsuccessful solver attempts.

`tables/all_paired_results.csv` contains every parameter, level, solver and target;
`best_by_method.json` selects the narrowest certified interval within each method's
retained grid. Every external upper bound is paired with the same fresh range-four
Jack--Sollich lower bound, charging its preparation and validation. Its transfer
verification timer includes both targets, so charging the complete timer to each
target is conservative. `campaign_costs.csv` also charges all retained repairs,
levels and failures. These are single-run measurements, not timing distributions.
Interpreter/import startup is excluded; certificate-factor construction is part
of generation. Peak memory includes imports and is the fresh worker's peak working
set. Sequentially paired workers use the maximum, not the sum, of their peaks.
Failed-worker peak memory was not captured and is explicitly unavailable.

For cost at a matched tolerance, `matched_tolerances.csv` allows every external
upper bound to choose a cheaper verified range-two/three/four, variance, or constant
trial lower bound. It is not forced to pay for range four. `all_cost_pairings.json`
records the exact endpoints and sum of measured component costs for every pairing.
The campaign ledger charges the complete eligible lower pool once per external
method/case, in addition to all upper-bound attempts and any map optimization.
The constant-trial formula and its measured validation cost are in `lower_sources/`.

Every reported method is also intersected with the elementary upper bound zero:
all cases have s >= 0, and L_s applied to the constant function is nonpositive.
`raw_upper` preserves the un-intersected certificate endpoint. The replay validates
the pairing, this uniform intersection, and every exact width in the main ledger.

`trial_certificate_overhead.json` isolates the extra graph-dual and upper-bound
verification cost beyond generating and validating each one-sided trial value.
Historical main-table times remain distinct and exclude candidate generation.

## Optional candidate generation

Install `requirements-generation.txt` in addition to the root requirements.
Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` for comparable single-worker
measurements. The recorded environment is in `ENVIRONMENT.json`.

The recorded Clarabel runs retained its internal `max_threads=0` automatic
default; actual solver-internal thread counts were not measured. The one-thread
BLAS/OMP environment must not be interpreted as a verified one-core constraint
on every solver.

For example:

```
python -B external/code/benchmark.py --method js --case strong --m 4 --output external/regenerated/js_strong_R4
python -B external/code/benchmark.py --method lti --case moderate --m 6 --solver SCS --output external/regenerated/lti_moderate_m6
python -B external/code/benchmark.py --method kull_rg --case weak --last-k 5 --solver SCS --output external/regenerated/kull_weak_k5
python -B external/code/benchmark.py --method variance --case weak --m 3 --output external/regenerated/variance_weak_R3
```

Outputs must be fresh and remain inside `external/`. Floating optimization need
not reproduce identical candidates across library versions or machines; every new
candidate is verified independently before it receives an endpoint. The RG command
uses the shipped range-two trial to build its MPS map, and charges that trial's
recorded generation time. The separate `kull_native` rows run the unmodified author MATLAB dual functions under R2026b, YALMIP R20250626_fix2 and SDPT3 4.0-20240410. The same exact verifier checks their rationalized output. Julia/MOSEK execution is not claimed.

The ordered `*_JOBS.json` lists show the entire experimental grid and the added
solver repairs. The suite launchers skip existing output directories; their purpose
is to document the actual run order, not to overwrite the frozen release.

## Scope

Anderson, Eisert's marginal program, Kull's uncompressed LTI hierarchy and a separate
bond-two renormalized implementation, a small zero-temperature FFS relaxation,
Jack--Sollich block trials, and variance-minimized trials are actually executed.
The upper endpoints use exact rational slack reconstruction plus interval-enclosed
transcendental coefficients. Raw primal objectives are diagnostics only.

The comparison does not reproduce optimized large-bond author benchmarks, the full
FFS hierarchy, positive-temperature equilibrium algorithms, or Wang's structured
NPA implementation. MPS/DMRG and stochastic simulations without matched deterministic
error bounds are not credited as equivalent certification tasks. No general method
dominance, new generic Gibbsianity result, or enlarged prior locality window follows.

Source metadata and author-code commits are in `SOURCE_RECORDS.json` and
`CODE_PDF_RECORDS.json`. Third-party article PDFs, repository snapshots and the local
Python runtime are kept in the independent extension workspace, not redistributed
inside this submission package.


## Native MATLAB regeneration

`code/native_kull.m` calls the pinned author's `Relaxation_dual_SDP`,
`declareDualVars` and `setDualConstraints` without changes. Only paths, model data,
solver options and output serialization are supplied by our wrapper. A work root
with `author_code/RDM_Constraints_Renormalization`, `matlab_runtime/YALMIP_release`
and `matlab_runtime/sdpt3_release` is needed for generation. Download records and
hashes are in `MATLAB_RELEASES.json`; the author commit is in `CODE_PDF_RECORDS.json`.
Use the exact sources to reconstruct these directories. Set `--work-root` to that
root and run `python external/code/run_native_kull.py --work-root YOUR_ROOT` in a
fresh copy with the existing `kull_native_*` result folders relocated. All three
parameters and all twelve native jobs are generated sequentially; no directory is
overwritten. Frozen certificate replay needs neither MATLAB nor YALMIP/SDPT3.

Core generation excludes process startup consistently with the Python comparisons;
native startup and I/O are separately recorded in each result. Native memory uses
the maximum aggregate resident set sampled every 50 ms for the Python worker and
MATLAB process tree, also taking the Python and MPS preparation peaks. This includes
startup but can miss transient peaks; it is not identical to an operating-system
exact aggregate peak. MATLAB was launched with `-singleCompThread`. No global path,
installation files or user settings were modified by the wrappers.

The explicit native startup ledger is `tables/native_startup_costs.csv`: MATLAB startup and I/O add 12.21--12.83 seconds to each fresh native run. Full Python cold-start times were not retained in structured form, so complete cold-start rankings are not claimed.

Native YALMIP return codes are preserved as `native_problem_0`, `native_problem_5`, and `native_problem_-1` (7, 4, and 1 runs). Nonzero native statuses are not labelled approximately optimal. Every exported point is checked independently; certificate success is distinct from solver optimality. The campaign ledger separates CVXPY inaccurate flags from native nonzero codes.
