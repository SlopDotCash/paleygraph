# Pass 48: three quartic triple fibers and an aggregate criterion

The full Paley and Proximity Prize goals remain unproved. This pass
refutes a proposed pointwise shortcut on actual quartic fields and
completely classifies it at fixed orders through128. No uniform exponent
improves.

## What changed

The primitive polynomial R_n was already constructed in the repository.
Root now uses

    G_n=gcd(|Res(R_n,R_n')|,|Res(R_n,R_n'')|)

to obtain a finite list containing every prime with a primitive fiber
of size at least three. Complete factorizations and checks of every
eligible factor give the exact lists through n128. This avoids treating
a finite prime scan as a complete classification.

At order128, precisely three eligible primes in [128^4/4,128^4] have
fibers larger than two:

| Prime | Fiber histogram | Balanced energy |
|---:|:---|---:|
| 77,796,353 | 27 singletons, one doubleton, one tripleton | 5,120 |
| 118,593,281 | 27 singletons, one doubleton, one tripleton | 5,120 |
| 181,312,129 | 27 singletons, one doubleton, one tripleton | 5,120 |

All other eligible primes in that interval have maximum fiber at most
two. The three exceptions therefore refute the proposed uniform
maximum-fiber-two condition, while the weaker energy bound B<=8192
still holds throughout this entire fixed-order interval. All smaller
dyadic levels also satisfy their maximum-fiber-two condition there.
This is a finite tower result, not a bound as its order grows.

For a general dyadic n=2k, let

    Y_n(p)=sum_(c>=3) c(c-2)

over the primitive root multiplicities. Root proves

    B_n<=2k^2+nY_n(p),   Y_n(p)<=v_p(G_n).

For a fixed prime splitting at every dyadic level up to N, this yields

    E_N<=44N^2+22N^2 V_N(p),
    V_N(p)=sum_(s=4,8,...,N)
             (3/4)^(log_2(N/s)) Y_s(p)/s.

A uniform V_N<<N^(1/3) in the quartic prime range would reach the
intermediate energy7/3 and triangle17/3 targets. That input is unproved,
as is the stronger sufficient condition with v_p(G_s) in place of Y_s.
The elementary resultant height bound remains too weak. The full goals
need further estimates even after the intermediate triangle target.
Proofs and an explicit triple representation are in the
[mathematical note](parallel48-collision-eliminants-2026-09-06.md).

## Verification and scope

The [checker](../experiments/parallel48_collision_eliminants.py) and
[certificate](../results/parallel48_collision_eliminants_2026_09_06.json)
cover five exact resultant pairs, five complete factorizations,66
distinct trial-division primality certificates and70 splitting-prime
cases. Polynomial resultants agree with a separate multiplication-matrix
determinant algorithm in every order; standard-library Bareiss agrees
through order64. Direct balanced and parent pair enumeration covers
1,056,320 pairs, in addition to the child-energy checks.

Root completed the ordinary derivations and exact checks. No Lean,
separate-agent, external review or novelty claim is made. The existing
parallel lanes were last observed at their usage limit; they were not
restarted. The compiled arithmetic backend is isolated in the project's
ignored work directory with its version pinned. The superseded standard-library
exploration finished; the exploratory order256 factorization was
terminated, with no order256 certificate claimed. No Lean process or
Prove2Me submission was changed. The
[audit](../results/parallel48_pass_audit_2026_09_06.json) records the
source, artifact and prior-state checks.

The next mathematical task is a uniform estimate on V_N or a different
sufficient aggregate. The pointwise maximum-fiber-two condition can no
longer be used as an unproved shortcut expected to hold throughout the
quartic range. Uniform W86/15, energy49/20, period71/72, absolute prime
exception5/3 and the full targets remain unchanged. The goal is active.
