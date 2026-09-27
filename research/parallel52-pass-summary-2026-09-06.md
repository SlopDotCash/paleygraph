# Pass 52: charge new primitive collisions to cubic increments

The Paley conjecture and Proximity Prize remain unproved. This pass
quantifies the cost of primitive collisions at each tower step, but
does not prove a new uniform bound.

Let X_s be the nontrivial shifted multiplicative energy at subgroup
order s, and Delta_s=X_s-X_(s/2). Four fine incidence cells merge into
one parent, so Delta_s counts exactly the ordered triples that first
join at this step. The primitive fibers occupy two equally sized
children of certain parent cells. Their exact contribution proves

    Delta_s >= 36Y_s,
    Y_s^2 <= s Delta_s/24.

The [full proof](parallel52-cubic-increments-2026-09-06.md) also bounds
the new normalized additive-energy excess by the same increment and
characterizes when the increment vanishes.

Summing by Holder gives

    sum_(M<s<=N) sqrt(Y_s) << N^(1/4)(X_N-X_M)^(1/4).

Thus an upper-interval increment bound X_N-X_M<<N^(5/3), at the
cutoff from pass50, would suffice for its intermediate target. This
condition remains unproved. The existing X_N<<N^2 log N theorem only
returns energy power5/2 through this argument, weaker than49/20.
Subtracting the smaller excess supplies no automatic power saving.

In the existing quartic example p2144280833,N256, the exact values are
X128=X256=1260. Rich cells persist without creating new cubic triples.
This refutes automatic fixed-factor growth at every eligible step; a
finite example does not exclude an eventual or additive-error variant.

The [checker](../experiments/parallel52_cubic_increments.py) independently
reconstructs repeated incidence cells and counts shifted products.
All368 distinct steps,71 towers and397 cutoff checks pass. Detailed
[results](../results/parallel52_cubic_increments_2026_09_06.json) and the
[audit](../results/parallel52_pass_audit_2026_09_06.json) record the
exact counts and source scope. No new prime classification is claimed.

The next task remains a quantitative saving for actual upper-tower
mass or its increment budget. Uniform energy49/20, triangle86/15,
period71/72, absolute exception5/3 and the full goals remain unchanged.
Root completed ordinary proofs and checks. No independent-agent,
Lean, external review or novelty claim is made. Previously usage-limited
agents were not restarted. No Lean process was polled or changed,
and no Prove2Me submission was made. The goal remains active.
