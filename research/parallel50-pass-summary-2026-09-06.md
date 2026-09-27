# Pass 50: endpoint energy and the upper tower

The Paley conjecture and Proximity Prize remain unproved. This pass
improves the reduction to the missing energy estimate; no uniform
exponent improves.

Let Y_s be the actual primitive triple-fiber mass at subgroup order s,
Ytot_N=sum_s Y_s, and S_N=sum_s sqrt(Y_s), with sums over dyadic levels.
Subtracting the two-representation baseline before Cauchy gives

    3N+6N Ytot_N <= E_N,
    sqrt(E_N) <= sqrt(14N^2-25N)+2sqrt(2N) S_N,
    E_N <= 28N^2+16N(log_2 N-1)Ytot_N.

Thus total triple mass is necessary and sufficient up to a logarithm
for the corresponding energy power. The radical bound improves the
previous sufficient weighted criterion, particularly for collisions at
smaller subgroup sizes. No estimate for this mass has been proved.

The recovery also works from any intermediate subgroup M. Reusing the
audited 49/20 energy bound lets us choose M as the largest dyadic order
below N^(80/87)/(1+log N)^(4/29). To reach the intermediate energy7/3
and triangle17/3 bounds, it now suffices to control

    sum_(M<s<=N) sqrt(Y_s(p)) << N^(2/3).

The smaller levels are already handled. This leaves about
(7/87)log_2 N levels, plus a logarithmic correction, and remains a
growing-order problem. The sufficient estimate above is still open,
and even the intermediate triangle target would not prove the full goal.

The [proof](parallel50-tower-mass-2026-09-06.md) gives the constants,
the source hypotheses, and abstract profiles showing the limitation of
the logarithmic recovery step. These profiles are not field examples.
The direct mass criterion also avoids the higher-precision and
cancellation costs of the stronger critical-content certificate.

The [checker](../experiments/parallel50_tower_mass.py) independently
counted 1,689,552 pairs in 70 previously certified prime-field towers.
All 390 level checks and 460 intermediate-cutoff checks passed, with
exact rational bounds for radicals. The
[results](../results/parallel50_tower_mass_2026_09_06.json) retain every
level and reproduce the three earlier quartic examples. This is finite
verification, not a uniform mass estimate or a new prime classification.

Root completed ordinary proofs and checks. No separate-agent, Lean,
external review or novelty claim is made. Previously usage-limited
agents were not restarted. No Lean process was polled or changed,
and no Prove2Me submission was made. The
[audit](../results/parallel50_pass_audit_2026_09_06.json) records
preserved dependencies, sources and artifact hashes.

The next task is a uniform estimate on actual triple mass in this upper
tower interval, or another sufficient input. Uniform W86/15,
energy49/20, period71/72, absolute prime exception5/3 and the full goals
remain unchanged. The goal remains active.
