# All-method comparison

Each external upper bound is paired with the same fresh R=4 Jack-Sollich lower bound. Generation and verification of that lower bound are charged. For every method, the universal bound lambda/n <= 0 for s >= 0 is intersected with the raw endpoint. Selected rows minimize certified width within the retained grid; all runs, both targets and failures remain in all_paired_results.csv. Times are one fresh run; no asymptotic cost claim.

| Case | Method | Level | Strict lower (outward) | Strict upper (outward) | Width | Gen s | Verify s | Total s | Peak MiB |
|---|---|---|---|---|---|---|---|---|---|---|
| strong | anderson | m=8 | -0.010754146351791 | -0.010403616037372 | 0.00035053 | 0.5296 | 0.07527 | 0.6049 | 82.4 |
| strong | eisert | m=4,t=2 | -0.010754146351791 | -0.010751669627685 | 2.47672e-06 | 0.5534 | 0.04364 | 0.597 | 277.0 |
| strong | lti | m=4 | -0.010754146351791 | -0.010754144477712 | 1.87408e-09 | 10.03 | 0.04217 | 10.07 | 270.0 |
| strong | kull_rg | D=2,k=5, block MPS | -0.010754146351791 | -0.010754140409101 | 5.94269e-09 | 51.81 | 0.08924 | 51.9 | 274.1 |
| strong | kull_native | D=2,k=5, optimized MPS | -0.010754146351791 | -0.010754118387092 | 2.79647e-08 | 2.285 | 0.09992 | 2.385 | 779.6 |
| strong | ffs_t0 | m=4 | -0.010754146351791 | -0.010754143156549 | 3.19524e-09 | 2.641 | 1.34 | 3.981 | 280.1 |
| strong | js_block_plus_certificate | R=4 | -0.010754146351791 | -0.010754144938594 | 1.4132e-09 | 0.5045 | 0.07178 | 0.5763 | 81.0 |
| strong | variance_plus_certificate | R=3 | -0.010754146433951 | -0.010754053525320 | 9.29086e-08 | 0.9296 | 0.02222 | 0.9519 | 194.5 |
| moderate | anderson | m=8 | -0.001644826501742 | -0.001466272281080 | 0.000178554 | 0.687 | 0.07557 | 0.7626 | 83.1 |
| moderate | eisert | m=4,t=3 | -0.001644826501742 | -0.001641552670000 | 3.27383e-06 | 2.025 | 0.07607 | 2.101 | 508.7 |
| moderate | lti | m=6 | -0.001644826501742 | -0.001644750969834 | 7.55319e-08 | 125.3 | 0.07874 | 125.4 | 274.6 |
| moderate | kull_rg | D=2,k=5, optimized MPS | -0.001644826501742 | -0.001644755521735 | 7.098e-08 | 63.42 | 0.1029 | 63.52 | 275.7 |
| moderate | kull_native | D=2,k=5, optimized MPS | -0.001644826501742 | -0.001644777414363 | 4.90874e-08 | 2.409 | 0.1014 | 2.51 | 774.7 |
| moderate | ffs_t0 | m=4 | -0.001644826501742 | -0.001644108942247 | 7.17559e-07 | 2.633 | 1.068 | 3.7 | 280.7 |
| moderate | js_block_plus_certificate | R=4 | -0.001644826501742 | -0.001596875287091 | 4.79512e-05 | 0.6614 | 0.07321 | 0.7346 | 81.7 |
| moderate | variance_plus_certificate | R=3 | -0.001645266565584 | -0.001432781248772 | 0.000212485 | 0.9219 | 0.02391 | 0.9458 | 194.5 |
| weak | anderson | m=8 | -0.000129383182438 | -0.000110697181153 | 1.8686e-05 | 0.7286 | 0.0761 | 0.8047 | 83.2 |
| weak | eisert | m=3,t=3 | -0.000129383182438 | -0.000128982002038 | 4.0118e-07 | 2.47 | 0.0781 | 2.548 | 503.0 |
| weak | lti | m=6 | -0.000129383182438 | -0.000129127145495 | 2.56037e-07 | 23.41 | 0.07634 | 23.49 | 275.3 |
| weak | kull_rg | D=2,k=5, optimized MPS | -0.000129383182438 | -0.000129109273655 | 2.73909e-07 | 64.5 | 0.1013 | 64.61 | 275.5 |
| weak | kull_native | D=2,k=5, optimized MPS | -0.000129383182438 | -0.000129196732556 | 1.8645e-07 | 3.256 | 0.1009 | 3.357 | 778.4 |
| weak | ffs_t0 | m=4 | -0.000129383182438 | -0.000128661801235 | 7.21381e-07 | 2.738 | 1.084 | 3.823 | 279.9 |
| weak | js_block_plus_certificate | R=4 | -0.000129383182438 | 0.000000000000000 | 0.000129383 | 0.7031 | 0.07238 | 0.7755 | 81.0 |
| weak | variance_plus_certificate | R=3 | -0.000129716242299 | -0.000081735024654 | 4.79812e-05 | 1.832 | 0.02298 | 1.854 | 195.5 |

Failed runs: 10. Paired rows: 228. Finite n=18 and limit rows are separate in the machine-readable ledger.

For extremely small width, the printed rounded endpoint difference can exceed the width column; the width uses the full 35-decimal endpoints.

## Cost at a common width tolerance

Cheapest eligible pair of measured upper/lower computations at each tolerance, selected after the grid for descriptive comparison. External methods may choose any fresh range-2/3/4 block, variance, or constant-trial lower endpoint; they are not forced to pay for range 4. This is not cumulative search cost; see the following ledger. NA means not attained within the retained grid. All three tolerances and both targets are in matched_tolerances.csv.

| Case | Method | Total s at width 1e-6 | Total s at width 1e-8 |
|---|---|---|---|
| strong | anderson | NA | NA |
| strong | eisert | NA | NA |
| strong | lti | 0.04085 | 9.626 |
| strong | kull_rg | 21.96 | 51.46 |
| strong | kull_native | 0.8086 | NA |
| strong | ffs_t0 | 0.4593 | 3.539 |
| strong | js_block_plus_certificate | 0.1113 | 0.5763 |
| strong | variance_plus_certificate | 0.9519 | NA |
| moderate | anderson | NA | NA |
| moderate | eisert | NA | NA |
| moderate | lti | 0.1636 | NA |
| moderate | kull_rg | 23.56 | NA |
| moderate | kull_native | 0.9796 | NA |
| moderate | ffs_t0 | 3.137 | NA |
| moderate | js_block_plus_certificate | NA | NA |
| moderate | variance_plus_certificate | NA | NA |
| weak | anderson | NA | NA |
| weak | eisert | 0.1802 | NA |
| weak | lti | 0.1464 | NA |
| weak | kull_rg | 23.74 | NA |
| weak | kull_native | 0.7817 | NA |
| weak | ffs_t0 | 3.2 | NA |
| weak | js_block_plus_certificate | NA | NA |
| weak | variance_plus_certificate | NA | NA |

## Cumulative campaign costs

This ledger includes every retained solver level, repair and failure; the eligible lower pool is charged once per external method/case, including map optimization where applicable. Timings exclude interpreter/import startup. Failed-worker peak memory was not recorded and is not inferred from successful runs.

| Case | Method | Completed certificates | No certificate | Solver warnings | Campaign s |
|---|---|---|---|---|---|
| strong | anderson | 4 | 0 | 0 | 1.681 |
| strong | eisert | 10 | 1 | 9 | 426.7 |
| strong | lti | 9 | 0 | 7 | 172.7 |
| strong | kull_rg | 5 | 3 | 5 | 222.7 |
| strong | kull_native | 4 | 0 | 2 | 6.289 |
| strong | ffs_t0 | 2 | 0 | 1 | 5.473 |
| strong | js_block_plus_certificate | 3 | 0 | 0 | 0.732 |
| strong | variance_plus_certificate | 1 | 0 | 0 | 0.9519 |
| moderate | anderson | 4 | 0 | 0 | 1.874 |
| moderate | eisert | 10 | 1 | 7 | 437.9 |
| moderate | lti | 9 | 0 | 4 | 131.1 |
| moderate | kull_rg | 5 | 2 | 5 | 225.6 |
| moderate | kull_native | 4 | 0 | 2 | 6.607 |
| moderate | ffs_t0 | 2 | 0 | 2 | 5.22 |
| moderate | js_block_plus_certificate | 3 | 0 | 0 | 0.9335 |
| moderate | variance_plus_certificate | 1 | 0 | 0 | 0.9458 |
| weak | anderson | 4 | 0 | 0 | 2.805 |
| weak | eisert | 10 | 1 | 7 | 463.1 |
| weak | lti | 9 | 0 | 3 | 30.84 |
| weak | kull_rg | 5 | 2 | 5 | 233.3 |
| weak | kull_native | 4 | 0 | 1 | 8.172 |
| weak | ffs_t0 | 2 | 0 | 2 | 6.255 |
| weak | js_block_plus_certificate | 3 | 0 | 0 | 0.9514 |
| weak | variance_plus_certificate | 1 | 0 | 0 | 1.854 |

Solver warnings count CVXPY optimal_inaccurate flags and native nonzero YALMIP problem codes. Native problem 5 or -1 is not called optimal_inaccurate; each raw code/info is retained independently of exact certificate success.
