COMPUTATIONAL SUPPLEMENT
Journal of Statistical Physics
Pointwise Perron locality and certified Doob dynamics for the softened East model
Minhao Liu
minhliu@ttu.edu

This archive contains the frozen computational evidence, certificate verifiers,
and numerical diagnostics accompanying the article and Supplementary Material.
Manuscript sources and figures are supplied separately in source_files.zip.
INDEX.txt maps the reported calculations to their data and executable checks.

SETUP AND INTEGRITY
Run the commands below from the extracted package root, after installing the
dependencies in requirements.txt. That file retains the recorded replay versions;
external/ENVIRONMENT.json and the stored runtime records identify the original
benchmark environments. A later validation environment must be reported separately.

  python -m pip install -r requirements.txt
  python -B code/check_manifest.py

MANIFEST.json lists the byte count and SHA-256 of every shipped file except itself.
The main replay commands validate this manifest before computing. Python assertions
are required: do not use -O or PYTHONOPTIMIZE. Use a new output name for every run.
Numerical inputs, values, endpoints and timing arrays are preserved. Obsolete local
paths in provenance metadata and logs have been made portable; PACKAGING_NOTES.txt
describes those changes. The current manifest covers this archive's exact bytes.

QUICKEST REPRESENTATIVE CERTIFICATE (SM9.2--SM9.3)
  python -B certificates/code/replay_submission.py --scope example --output replay_example

ALL SEVEN ROWS OF THE FULL SM9.4 TABLE, INCLUDING MAIN TABLE 1
  python -B certificates/code/replay_submission.py --scope table1 --output replay_table1

ALL 22 GRAPH AND 20 QUADRATIC CERTIFICATES
  python -B certificates/code/replay_submission.py --scope all --output replay_all

These commands create new subdirectories of certificates/. They recompute graph
endpoints, local constants, quadratic corrections and strict infinite-volume centers
from the frozen candidates and positive seeds. Published endpoints are compared
using exact rational arithmetic. The worked 30-decimal certificate is regenerated
from canonical endpoints. Inspect REPLAY_REPORT.json and endpoints.csv. Range-nine
verification takes substantially longer than the representative example. Replay
elapsed times do not replace the original three-run timing measurements.

MONTHUS COUNTEREXAMPLE (SM11; PYTHON STANDARD LIBRARY)
  python -B monthus/code/replay_monthus.py
  python -B monthus/code/replay_independent.py

The independent implementations reconstruct the finite-ring discrepancy and check
the rational all-volume witnesses. Their default fresh output directories are
monthus/replay_primary/ and monthus/replay_independent/. Use --output with a new
name to repeat. See monthus/README.txt for the exact model correspondence and scope.

EXTERNAL COMPARISON (SM12)
  python -B external/code/replay_external.py --output replay_external

This creates external/replay_external/ and checks the retained Anderson, marginal,
LTI, renormalized SDP, zero-temperature equilibrium and block/variance certificates.
Candidate optimization is not repeated. Exact inputs, raw solver statuses, repairs,
failed attempts, timings and native MATLAB exports are retained. See external/README.md
for regeneration dependencies and external/THEORY_AUDIT.md for the model mappings.
Frozen certificate replay requires no MATLAB, SDP solver, network or vector cache.

BASELINE GRAPH, LOCALITY AND PLOTTED-BOUND CHECKS (SM4--SM8)
  python -B baseline/code/reproduce_v14.py --output baseline/replay_baseline

The stable script name identifies the original computational implementation. It
copies evidence into a fresh directory, then checks 92 graph certificates, six
constant trials, 32 locality certificates, 24 intersections, 74 tail-curve points,
seven reference-bound points, 16 parameter cases and output-precision checks.
Original inputs and published results stay unchanged.

Seventeen reference-vector diagnostic summaries are retained with their provenance.
Their large reference vectors were not part of the original compact archive. Replay
reports reference_vectors_rechecked=0 and stored_reference_records_retained=17.
Graph/locality/quadratic certificate replay uses frozen decimal candidates directly.

COST/ACCURACY AND LOW-DIMENSIONAL CHECKS
  python -B certificates/code/check_cost_accuracy.py
  python -B certificates/checks/closed_form_checks.py
  python -B certificates/checks/reversible_checks.py

The cost check recomputes widths and medians from exact endpoints and stored timing
arrays, including the common-target comparison; it performs no new benchmark.
Two-state and finite-state implementation diagnostics are identified separately
from the rigorously enclosed spectral endpoints.

ONE-COMMAND REPLAY OF THE CORE FROZEN CERTIFICATES
  python -B code/replay_all.py --tag replay_complete

This runs manifest verification, the external replay, all graph/quadratic cases,
both Monthus implementations, cost and low-dimensional checks, and the baseline
replay. It creates replay_logs/replay_complete/ and fresh per-module output folders.
The additional diagnostics below are separate from this core replay.

LOCAL-OBSERVABLE AND SOFTENING-SCALE DIAGNOSTICS
The recorded diagnostic dependency versions are in
validation/local_observables/requirements-diagnostics.txt. With those dependencies:
  python -B validation/local_observables/validate_new_results.py --output replay_diagnostics

This writes a new subdirectory of validation/local_observables/. It checks finite-state
coupling/semigroup inequalities and sample sufficient softening-scale inequalities.
These checks support implementation consistency and are not a proof of an
infinite-volume locality theorem.

DOOB AND ALGEBRAIC DIAGNOSTICS
Install validation/doob_diagnostics/requirements-diagnostics.txt separately if needed.
  python -B validation/doob_diagnostics/check_new_math.py --output replay_new_math.json
  python -B validation/doob_diagnostics/check_doob_identity.py --output replay_doob_identity.json

The exact finite algebra and the floating finite-state diagnostics are distinguished
in their output. DOOB_CERTIFICATE.json retains the strict representative entropy
bound, with the assumptions and initialization specified in the article.

INDEPENDENT IMPLEMENTATION CHECKS (OPTIONAL; MAY BE EXPENSIVE)
  python -B validation/independent_math/independent_checks.py --output replay_independent_math
  python -B validation/independent_math/support_metadata_checks.py --output replay_support_checks
  python -B validation/independent_certificates/independent_audit.py --output validation/independent_certificates/replay_exact
  python -B validation/independent_certificates/local_center_audit.py --output validation/independent_certificates/replay_centers
  python -B validation/independent_certificates/cost_provenance_audit.py --output validation/independent_certificates/replay_costs

These are implementation-level cross-checks with retained numerical records, not
external peer-review reports. All output directories must be new. Recorded numerical
outputs retain their original script hashes and environments. Relocation and
output-safety edits to the scripts are described in PACKAGING_NOTES.txt.

OPTIONAL FINITE CENTERS AND FIGURES
  python -B certificates/code/finite_variational.py --dps 120 --constants certificates/published/finite_constants --output certificates/replay_finite

For figure redrawing, add matplotlib (originally tested with 3.10.0):
  python -B baseline/replay_baseline/code/figures_v14.py
  python -B certificates/code/plot_ranges.py --data certificates/published/range_study.json --output certificates/replay_figures/quadratic_range

The first figure command follows a successful baseline replay. Font substitution
changes rendering only. Each finite-center record applies to its specified ring.
Graph intervals cover every n >= 2R-1; quadratic infinite-center intervals certify
the thermodynamic spectral density. These scopes remain distinct throughout.
