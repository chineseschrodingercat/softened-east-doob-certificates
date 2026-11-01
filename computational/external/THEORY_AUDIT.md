# External methods: definitions, assumptions and certificate derivations

Status: working audit, 12 October 2026. Original v19 proofs are preserved. All new experiments retain c=0.3 and the three existing (kappa,s) pairs. No numerical experiment extends the prior locality window.

## Sources and version lock

| Work | Version checked | Publication | Code and execution scope |
|---|---|---|---|
| Eisert | arXiv:2301.06142v3, 14 Nov 2025; formal PDF read through APS/web | PRA 113, 022214, 17 Feb 2026; DOI 10.1103/cz6k-y46r | No author implementation identified in the paper; equation-based implementation |
| Kull, Schuch, Dive, Navascues | arXiv:2212.03014v3, 10 Apr 2024, full PDF | PRX 14, 021008, 9 Apr 2024 | IlyaKull/RDM_Constraints_Renormalization, commit 2e9015fff5d9bc5b170cdc6cee98fbbb928decda; Python implementation of the published primal constraints plus unmodified MATLAB author dual routines, using the same MPS-product maps and exact invertible gauge |
| Jack and Sollich | arXiv:1308.6124v2, 12 Dec 2013, full PDF | J. Phys. A 47, 015003 (2014) | Section 6.1.1 and Appendix A.1 implemented for the present generator; no author repository identified |
| Fawzi, Fawzi, Scalet | arXiv:2311.18706v3, 14 Aug 2025, full PDF | Nature Communications 15, 7394 (2024) | Actual zero-temperature, single-site perturbation SDP on 3/4-site marginals; not the complete hierarchy or positive-temperature implementation |
| Wang et al. | arXiv:2310.05844v5, 18 Jul 2024, full PDF | PRX 14, 031006 (2024) | wangjie212/ManyBodySOS commit 0b152dfcd9bcde3403aa4d967642ce5261c80712; README now names QMBCertify, Julia/JuMP/MOSEK and Heisenberg models; no native execution claimed |
| De Roeck, Maes, Netocny, Schutz (DMNS) | arXiv:1312.4782v2, 8 Feb 2015, full PDF | J. Math. Phys. 56, 023301 (2015) | Theoretical comparison, Assumptions 1/2 and Theorems 4.2/7.6 |
| Aragones-Soria, Aberg, Park, Kastoryano | arXiv:2010.11643v2, 17 Sep 2021, full PDF | J. Math. Phys. 62, 093511 (2021) | Proposition 1, Definition 3, Theorem 4; theoretical comparison |

Latest arXiv HTML/PDF, metadata and hashes are under sources/. Author repositories are archived locally under author_code/ and sources/, not redistributed in the release. The arXiv version of Eisert equation (28) prints `2+1`; the formal PRA PDF corrects this to `s+1`. The implementation uses the formal equation. Publication pages for PRA and both PRX articles were checked. The latest FFS preprint postdates the journal article; theorem numbering here explicitly follows v3. Failed command-line publisher downloads are retained in PUBLISHER_RECORDS.json; access to the formal Eisert PDF through the web reader and full article in Chrome was successful. No claim is made that all publisher PDFs were saved locally.

## Common operator and targets

With row-action generator L, pi=Ber(c)^n, D=diag(sqrt(pi)), define H=-D L D^-1. Its local two-site term is

    h = diag(kappa,kappa+1) tensor [[c,-a],[-a,1-c]],
    a = exp(-s) sqrt(c(1-c)).

Thus E0(H_n)/n=-lambda_n/n. An energy lower bound is a spectral upper bound. There is no factor 2, e^s, 1+kappa or n hidden in the reported units. Periodic n=18 and the thermodynamic limit are distinct columns. Translation averaging of a ground-state density matrix supplies feasible translation-invariant marginals even if a finite-volume ground state is degenerate. The common upper certificates are valid for every ring at least as long as their largest physical support; each certificate records n_min. Infinite variational centers use verified transfer eigenvectors; finite centers use the exact trace of the 18th transfer power.

Jack-Sollich originally considers the hard East model and a bias nu of the time-integrated escape rate. For the off-diagonal activity bias used here, L_s=exp(-s) W_R(1-exp(s)), with the same softening inserted in the rates. Its variational free energy F therefore gives lambda_s/n >= -exp(-s) F. We implement the equivalent Rayleigh expression directly in the original L_s units. The overlapping-block ansatz, transfer matrix, and one-sided variational principle are prior work, not new contributions. We do not transplant the original hard-model numerical conclusions to positive kappa.

## Validated external bounds

All input c,kappa,s and dual entries are exact rationals. An interval evaluation at 100 decimal digits encloses a. A rational approximation a0 is used by the solver; delta=(kappa+1) max(|a0-a_lower|,|a0-a_upper|) bounds the local Hamiltonian error. Every energy-density lower bound loses at least delta. This correction also covers square-root and exponential coefficient errors.

For any symmetric rational matrix S, a stored rational factor B and rational shift mu give

    lambda_min(S) >= mu - max_i sum_j |S_ij - mu I_ij - (BB^T)_ij|.

The row sums are evaluated with integer matrix products and exact fractions. Eigenvalues and Cholesky decompositions generate B only; their floating accuracy is not assumed. This is an elementary PSD residual correction, not a new theorem. A solver can fail or return optimal_inaccurate while a supplied dual still yields a valid, possibly weaker, bound.

**Anderson.** Let H_m be the sum of m-1 open bonds. Averaging all translates of H_m on any ring n>=m gives (m-1)H_n. Hence E0(H_n)/n >= E0(H_m)/(m-1). For the stoquastic open block, a positive rational vector v gives E0(H_m)>=min_i (H_m v)_i/v_i; both invariant sectors of the leftmost spin are solved. No unverified Ritz value enters the certificate.

**Eisert improved marginal program.** For overlap t, use trace-one PSD omega on m sites and sigma on 2t sites. Match the first/last t-site marginals of sigma to the last/first marginals of omega, and minimize 2 tr(H_m omega)+tr(h_central sigma). Any translated true state is feasible and the objective equals (2m-1)e. This directly proves the required lower bound, without relying on equality to a finite open 2m-site optimum. Dual matrices Y_L,Y_R give two slacks

    S_omega=2H_m-embed_last(Y_L)-embed_first(Y_R),
    S_sigma=h_central+Y_L tensor I+I tensor Y_R.

Certified minimum eigenvalues gamma_omega,gamma_sigma yield e >= (gamma_omega+gamma_sigma)/(2m-1)-delta. We use the displayed published program on TI states; we do not require its strongest identification with the unconstrained open-patch spectrum. The t hierarchy need not numerically improve monotonically after different solver residuals; all levels are retained.

**LTI and Kull renormalization.** LTI minimizes tr(h rho) over an m-site PSD trace-one marginal with equal left/right (m-1)-site reductions. Its dual slack h+I tensor Y-Y tensor I telescopes on a ring. This is the uncompressed hierarchy of Kull equation (5), not the RG algorithm. The separate RG experiments implement the author primal constraints: a four-site rho, compressed omega_k on physical-left/bond-pair/physical-right, V0 rho V0* matching both partial traces of omega_3, then L/R omega_k L/R* matching omega_(k+1). The MPS bond dimension is D=2 and k=3,5,7. This modest bond budget is not a benchmark of the authors' optimized larger-D implementation.

The map V_k sends a physical word to the vectorized product of its rational MPS matrices. L and R multiply one matrix on the left and right, so the composition identities hold algebraically. A rational invertible change of basis W improves conditioning; its inverse is computed exactly, and V_k becomes W V_k. No rounded approximate isometry is assumed. Compressed states need not have trace one: for a physical density matrix, tr(omega_k)<=||W V_k||_F^2=:B_k. If reconstructed dual slacks have certified minima gamma_0,gamma_k, then

    e >= gamma_0 + sum_k min(gamma_k,0) B_k - delta.

This repairs every dual PSD defect, including intermediate compressed constraints. The author code also provides FeasibleDualPt.m; this experiment replaces its numerical eigenvalue correction by an exact rational factor residual. A failure of our solver or map choice does not refute RG convergence or performance with better maps.

The initial maps use the optimized range-two block trial. To test this restriction, a further k=5 run at every parameter uses general real D=2 MPS tensors optimized through their four-dimensional transfer operator from three deterministic starts, retaining the original block tensors as a fallback. Their floating energy values are only map-generation diagnostics. Rounded rational tensors define the actual maps; the same exact dual verification supplies the endpoint. General MPS-map optimization and the initial block trial are included in generation costs. Larger D remains untested. After MATLAB R2026b became available, the unmodified author dual functions were run with YALMIP R20250626_fix2 and SDPT3 4.0-20240410 on all three cases, k=3,5,7 for block maps and k=5 for general optimized D=2 maps. These native runs are separately labelled kull_native. Exported duals receive the same exact rational correction; solver values are diagnostic only. The default automatic HSD mode failed a 2-by-2 environment probe; standard mode 1 passed and was frozen for every native run. No global MATLAB path was saved. Native startup time is recorded separately, and peak memory is the maximum sampled aggregate Python worker plus MATLAB process tree (50 ms sampling).

**FFS at zero temperature.** On a marginal containing the support and both neighbors of a one-site perturbation, enforce stationarity and C_ij=tr(rho a_i* [H,a_j]) PSD for the four one-site matrix units, in addition to LTI. Minimize energy. This is a small, valid restriction of their ground-state conditions, not their full hierarchy. Its PSD dual is replaced by an exact rational Gram matrix, and the slack is certified as above. Since the ground-state conditions contain H, coefficient error must also be propagated through the dual: each C_ij changes by at most 4 delta, each one-site stationarity expression by at most 4 delta, and each Hermiticity difference by at most 8 delta. The implementation explicitly adds these dual-weighted errors. The finite-temperature matrix relative-entropy constraints are irrelevant to the present zero-temperature spectral target and are not implemented.

## What is and is not new in the convergence comparison

The elementary open-block interval already proves computable termination and O(1/m) boundary budgets without a spectral gap. Eisert gives an Anderson performance budget; LTI/SDP hierarchies can improve far faster in individual models. O(1/m) is not an impossibility theorem for those methods. Kull's paper demonstrates much stronger behavior after compression, depending on the maps; its empirical rates for other spin models are not universal East-model rates.

FFS Theorems 1.1/1.2 give convergent two-sided observable bounds for general, including noncommuting, equilibrium states. The quantitative exponential theorem 1.5 assumes commuting interactions and strictly positive temperature in one dimension (beta finite); it does not give a zero-temperature rate for the noncommuting softened East Hamiltonian. Our model-specific exponential flip-error and range budget is a different conclusion, on its stated small-tilt window. General certified convergence is already known. Wang et al. combine a variational energy upper bound and an SDP lower bound to constrain other ground-state observables. For O=energy their energy-window step does not add a new target beyond the energy interval already being compared. Their native code would need a new East Hamiltonian/model builder and a Julia/MOSEK environment; this run does not claim to reproduce it. LTI is not labelled a reproduction of Wang's structured NPA numerics.

DMNS supplies exponentially summable classical potentials under a product-vacuum Peierls condition, excitation-sector commutation, small perturbation and non-null measurement overlaps. Such a potential already gives configuration-uniform flip control by summing its tails. Thus neither Gibbsianity nor this implication is new here. For softened East, the relative term remains O(1) at zero tilt, and absorbing it into their reference breaks the required excitation-pattern commutation; main Section 3.1 computes the nonzero commutator. Our contribution is the specific physical supremum-norm estimate, periodic matching and certificate construction under the paper's hypotheses, not a claim that all alternative applications of DMNS are impossible.

The MPS quasi-local Gibbs result assumes an injective finite-bond MPS, fixed measurement basis and purity, and obtains exponential decay of (classical/post-measurement) conditional mutual information. Proposition 1 controls relative entropy to a local Gibbs approximation by a sum of CMI terms. This is distribution-weighted control, not a bound on every configuration's log-amplitude flip ratio. A small Hilbert-space, trace-norm or relative-entropy error can hide a large relative error on rare configurations; no volume-independent lower bound on all probabilities is supplied here. It also does not establish that the exact East Perron vector is a fixed-bond injective MPS satisfying purity. These are distinct assumptions and norms, not a dismissal of the MPS theorem.

## Additional candidate methods and exclusion record

- Positive finite-range trials are a structured MPS subclass. General real D=2 MPS tensors are also optimized to generate RG maps, but their floating energy is not promoted to a certified variational endpoint. Unrestricted larger-D MPS/DMRG certification is not run; the exact coefficient/positivity and certificate interfaces would need additional validated contractions. No general MPS speed or accuracy claim follows.
- Variance minimization is run on the same n=18 target with the same range-three positive ansatz, using exact state enumeration for a floating candidate objective. The resulting candidate is independently certified by the same graph and quadratic verifiers. Small variance alone is never accepted as a rigorous eigenvalue error bound.
- The same positive candidate defines approximate Doob rates exp(-s) r_i h(x^i)/h(x). Their diagonal compensator is exactly the residual being certified. Thus the graph bounds add a deterministic compensator enclosure to this candidate family; this is not a new independent SCGF estimator or a claim of certified path-space accuracy.
- Importance sampling/cloning estimates are stochastic and have finite-time/population errors. Without a proved confidence/finite-time-to-eigenvalue bound at the same target they do not supply the requested deterministic endpoints. No simulation is billed as a competing rigorous interval. Implementing that additional theory would exceed this round's fixed theoretical scope.

The Monthus counterexample remains independent: only the stated exact two-spin activity-closure consequence at c=1/2,kappa=1,exp(-s)=1/2 is excluded. It is neither a general refutation of that paper nor evidence enlarging the prior locality window.
