# Pointwise Perron locality and certified Doob dynamics for the softened East model

**Computational materials by Minhao Liu, Texas Tech University**  
Contact: [minhliu@ttu.edu](mailto:minhliu@ttu.edu)

This repository accompanies the research manuscript *Pointwise Perron locality and certified Doob dynamics for the softened East model*. It provides the fixed numerical inputs, certificate verifiers, results, and supporting implementation checks used in the study. The work connects pointwise locality of Perron eigenfunction ratios to rigorous spectral and dynamical approximation of the activity-biased softened East model.

## Download the fixed release

[**Download computational_supplement.zip — fixed release v1.0.0**](https://github.com/chineseschrodingercat/softened-east-doob-certificates/releases/download/v1.0.0/computational_supplement.zip)

[Browse the files at this release](https://github.com/chineseschrodingercat/softened-east-doob-certificates/tree/v1.0.0/computational) · [Release notes](https://github.com/chineseschrodingercat/softened-east-doob-certificates/releases/tag/v1.0.0)

Use the named release above to obtain the computational archive accompanying the manuscript. Its SHA-256 is:

```text
097d164e2cc05166564d12cecfc33b8bf4e77079721337528d033838ba54cd9b
```

The release attachment preserves the original ZIP bytes. The repository's [`computational/`](computational/) directory contains the unpacked contents of that archive; its [`MANIFEST.json`](computational/MANIFEST.json) records the byte count and SHA-256 of each supplied file except the manifest itself. The archive checksum identifies the downloadable ZIP, while the manifest identifies its individual contents.

## What is here

The package contains independently executable checks of fixed computational evidence:

- Graph certificates enclosing the spectral density on all rings above a stated minimum size.
- Quadratic residual certificates, with rigorously enclosed variational centers and verified local constants.
- Locality certificates, spatial-tail bounds, and range/reference-volume studies.
- Exact-rational witnesses for the two-spin closure counterexample.
- Comparisons with finite-block, marginal, and semidefinite methods, including unsuccessful runs and solver-status records.
- Finite-state and algebraic diagnostics supporting the implementation of the Doob, local-observable, and softening-scale calculations.

The full instructions are in [`computational/README.txt`](computational/README.txt). [`computational/INDEX.txt`](computational/INDEX.txt) maps manuscript and Supplementary Material results to data files and checks. A `PASS` report refers to the verification steps and inputs named in that report.

## Quick start in three steps

### 1. Download and enter the package directory

Download `computational_supplement.zip` from the fixed release above and extract it into a new directory. Open a terminal in the extracted directory containing `README.txt`, `requirements.txt`, and `MANIFEST.json`.

If you obtained the repository instead, enter its package directory:

```bash
cd computational
```

All commands below run from this package directory.

### 2. Install the recorded dependencies and check the files

Use a Python environment with the dependencies recorded in the package:

```bash
python -m pip install -r requirements.txt
python -B code/check_manifest.py
```

Python assertions must remain enabled: run without `-O` and without `PYTHONOPTIMIZE`. The recorded benchmark environments are described separately in the package; a replay on another environment does not replace those original measurements.

### 3. Recompute the representative certificate

```bash
python -B certificates/code/replay_submission.py --scope example --output replay_example
```

Inspect `certificates/replay_example/REPLAY_REPORT.json` and `endpoints.csv`. This example recomputes the graph endpoints, local constants, quadratic correction, and enclosed infinite-volume center from the supplied candidate and positive seeds. It also checks the endpoints using exact rational arithmetic. Choose a new output name for every repeat; existing output directories are rejected.

After dependency installation, replaying the fixed certificates requires no network access, MATLAB, SDP solver, candidate optimization, or reference-vector cache.

## Directory guide

| Location | Contents and purpose |
|---|---|
| [`computational/README.txt`](computational/README.txt) | Setup, replay commands, optional checks, and dependency guidance. |
| [`computational/INDEX.txt`](computational/INDEX.txt) | Mapping from the manuscript and Supplementary Material to supporting files. |
| [`computational/MANIFEST.json`](computational/MANIFEST.json) | File sizes and SHA-256 checksums for the supplied package. |
| [`computational/certificates/`](computational/certificates/) | Main-table and SM9 graph/quadratic cases, positive seeds, transfer centers, exact endpoints, and recorded timing arrays. |
| [`computational/baseline/`](computational/baseline/) | SM4–SM8 graph and locality certificates, constant trials, intersections, spatial tails, and range/reference-volume studies. |
| [`computational/monthus/`](computational/monthus/) | SM11 exact-rational counterexample, factor tables, and two replay implementations. |
| [`computational/external/`](computational/external/) | SM12 external-method inputs, solver outputs, feasibility corrections, costs, failures, and certificate replay. |
| [`computational/validation/`](computational/validation/) | Algebraic and finite-state diagnostics and additional implementation cross-checks. |
| [`computational/validation/doob_diagnostics/DOOB_CERTIFICATE.json`](computational/validation/doob_diagnostics/DOOB_CERTIFICATE.json) | The representative path-entropy certificate and its stated assumptions. |

## Additional replay entry points

Run these from the package directory. Each run needs fresh output names.

**All 22 graph and 20 quadratic certificates:**

```bash
python -B certificates/code/replay_submission.py --scope all --output replay_all
```

**The exact-rational closure counterexample, using two implementations:**

```bash
python -B monthus/code/replay_monthus.py
python -B monthus/code/replay_independent.py
```

**The external-method certificates:**

```bash
python -B external/code/replay_external.py --output replay_external
```

**The core certificate suite, including the baseline checks:**

```bash
python -B code/replay_all.py --tag replay_complete
```

The core suite covers the manifest, external certificates, all graph/quadratic cases, both closure-counterexample implementations, cost/accuracy checks, low-dimensional checks, and baseline replay. Additional diagnostics and candidate regeneration have separate entry points and dependencies documented in [`README.txt`](computational/README.txt). Range-nine cases and the complete suite take longer than the representative example.

## What the certificates establish

The manuscript supplies the mathematical hypotheses and proofs. The computational records distinguish the following targets:

| Result | Scope |
|---|---|
| **Graph interval** | A single interval encloses `lambda_n(s)/n` for every `n >= 2R - 1` and therefore also its thermodynamic limit. The candidate range is `R`. |
| **Quadratic interval using an infinite-chain center** | Encloses the thermodynamic spectral density, with the local residual, rate, and influence hypotheses checked for the candidate. |
| **Finite-ring center** | Applies to the particular ring size specified in its record. Finite-ring path-entropy statements use this center. |
| **Path relative entropy** | The spectral error gives candidate-to-exact Doob path relative entropy per site and unit time when both processes start from the candidate stationary law, in the original physical clock. |
| **Local-observable dynamics** | The analytical finite-time bound requires the stated locality and flip-gradient hypotheses. Finite-state diagnostics check implementation consistency; their success does not prove the infinite-volume theorem or establish its hypotheses at every posterior numerical example. |

Decimal inputs treated as exact rationals, rational witness checks, and outward interval enclosures underpin the certificates. Floating eigensolves and optimization generate candidates or provide diagnostics; raw floating objectives and solver status flags do not serve as certified endpoints. The package retains failures and corrections alongside successful cases.

Recorded timing arrays belong to the original runs. Recomputing their widths and medians checks those records; the elapsed time of a new replay is a separate measurement.

## Scope of the 17 historical reference-vector records

The archive retains **17 historical summaries of reference-vector diagnostics**. The large reference vectors underlying those summaries were not included in the original compact archive and are absent from this release. Baseline replay therefore reports:

```text
reference_vectors_rechecked = 0
stored_reference_records_retained = 17
```

These summaries document the earlier calculations, but their vector-based computations are not re-executed by this release. The graph, locality, and quadratic certificate replays instead use the supplied fixed decimal candidates and positive seeds directly, so they can be recomputed without those large vectors. [`INDEX.txt`](computational/INDEX.txt) identifies the retained records and separates them from the replayable certificates.

## Identifying this version

When referring to these computational materials, include the manuscript title, author, repository URL, fixed release tag, and ZIP checksum. The repository is a companion computational artifact for the research manuscript; its description and verification records do not assert journal acceptance or completion of peer review.

```text
Minhao Liu
Pointwise Perron locality and certified Doob dynamics for the softened East model
Computational materials, release v1.0.0
https://github.com/chineseschrodingercat/softened-east-doob-certificates
Archive SHA-256: 097d164e2cc05166564d12cecfc33b8bf4e77079721337528d033838ba54cd9b
```
