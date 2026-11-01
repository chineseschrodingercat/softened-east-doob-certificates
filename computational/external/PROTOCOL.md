# Frozen external comparison, 2026-10-12

This extension preserves v19 and does not extend the model or the proved tilt window.

- Generator: equation (A.1) of v19, c=0.3; (kappa,s)=(1,0.02),(0.1,0.01),(0.01,0.001), labelled strong, moderate, weak softening. Time is the original rate unit, not a rescaled Hamiltonian unit.
- Targets: largest eigenvalue per site on a periodic ring (primary n=18) and thermodynamic limit, kept separate. Volume-uniform bounds record the minimum n. No finite-ring variational value is silently used as an infinite-volume lower endpoint.
- Symmetric Hamiltonian H=-diag(sqrt(pi)) L diag(1/sqrt(pi)); local term diag(kappa,kappa+1) tensor [[c,-exp(-s)sqrt(c(1-c))],[-exp(-s)sqrt(c(1-c)),1-c]]. Energy bounds change sign and reverse order to give spectral bounds.
- Primary baselines: Anderson open blocks; Eisert improved Anderson marginal SDP (with any correction to the printed formulation explicitly documented); Kull author-code investigation followed by runnable RG/LTI relaxation; Jack-Sollich overlapping-block trial with transfer expectation, adapted to the specified softening and activity tilt. Compare the identical trial before/after graph certification.
- Screen Fawzi-Fawzi-Scalet, Wang et al., MPS, approximate Doob, variance minimization and importance sampling. Implement comparable affordable deterministic cases; give explicit reasons otherwise. No stochastic confidence interval is relabelled a deterministic certificate.
- Candidate generation/optimization and deterministic verification are separately timed with perf_counter. Total includes both; fresh subprocess peak working-set memory is reported. Frozen-candidate replay is distinguished from end-to-end generation. One measured run is not a benchmark distribution.
- All three parameter points appear in every primary comparison, with budget failures retained. Start with block sizes 2,4,6,8; SDP matrices up to 64 unless a pilot supports expansion; local trials R=2,3,4 and archived R=3,5,7 where useful. Pilot-driven budget changes must be logged, not selected by outcome.
- Floating eigensolver/SDP output is a candidate only. Reconstruct dual identities and bound all coefficient/PSD errors with interval or exact rational arithmetic. Negative certificate margins trigger repair or failure, never an unqualified endpoint.
- Keep old theoretical results attributed. Compare generic O(1/m) block budgets and the present model-specific conditional exponential range bound; the posterior examples lie outside the proved prior window.
- After results stabilize, two independent AI reviewers inspect mathematics and certificates, independently reconstruct selected checks, and replay the actual release. Repair and rerun before updating root deliverables.

Native-runtime addition after the user announced MATLAB availability: run the unmodified author dual code at every original parameter with D=2 block maps k=3,5,7 and optimized general D=2 maps k=5. No parameter/model/window changed. Record startup separately; retain every native problem code and certify exported points independently.
