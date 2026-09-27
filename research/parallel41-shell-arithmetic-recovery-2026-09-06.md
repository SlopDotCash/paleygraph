# Real-subfield recovery removes the scalar reciprocal exponent

An additional real-subfield certificate recovers a norm cofactor
linear in the coset's norm defect. If `F` is an actual coset polynomial
with `Norm(F)=p^(N-1)*tau`, and a real element `h` in the appropriate
evaluation kernel has real norm `lambda*p`, then

    f=h*bar(F)/p in R,       Norm(f)=p*lambda^2*tau.

The hypotheses and evaluation root are made precise below. A real
norm-`p` element therefore replaces the scalar-reciprocal cofactor
`tau^(N-1)` from
[pass40](parallel40-shell-converse-2026-09-06.md) by `tau`.

For `p=6700417,n=64`, a concrete real generator gives cofactor 1217
from the known short actual coset, replacing `1217^31`. It does not
improve the already known cofactor 641. A separate exact example
shows why the associated lattice correspondence does not solve
cofactor minimization by centering: centering lowers the squared
radius while raising the norm defect from 641 to 178771841.

The uniform thin-annulus estimate and the full goals remain open.
This is an ordinary mathematical argument and bounded exact
verification, with no build, proof-status change, or novelty claim.

## 1. The quantities that must remain separate

Let `n=2N>=4` be a power of two, let `H=<g>` have order `n` modulo
the prime `p`, and let

    R=Z[X]/(X^N+1),           bar(X)=X^(-1),
    F=sum_(j=0)^(N-1) r_j X^j,

where `r_j` is centered and congruent to `a g^j`, with `a!=0`.
Set `u=g^(-1)`. Modulo `p`, the only nonzero split-root evaluation
of `F` is at `u`. As before,

    V=sum_j r_j^2/p^2,
    tau=Norm(F)/p^(N-1) is a positive integer.

The squared radius `V`, the multiplicative norm defect `tau`, and
the smallest cofactor of **any** element of the evaluation kernel
are different quantities. A flat identity `F bar(F)=p^2` is yet
another condition and is not assumed here.

Searches of the research notes for real-kernel/cofactor recovery
found the preceding scalar reciprocal result and its flat-Gram
discussion, but no transfer through a real norm generator.

## 2. Recovery through a real norm certificate

Put `D=N/2` and let `R^+` be the subring fixed by conjugation.
It has integer basis

    1, s_j=X^j+X^(-j), 1<=j<D.

For `h in R^+`, let `Norm_+(h)` be the determinant of its integer
multiplication matrix in this basis, equivalently its norm in the
real cyclotomic subfield. Thus `Norm(h)=Norm_+(h)^2`.

Suppose a nonzero `h in R^+` satisfies

    h(u)=0 mod p,
    |Norm_+(h)|=lambda*p,
    lambda>=1 and p does not divide lambda.                         (1)

Then

    f=h*bar(F)/p belongs to R,
    f(u)=0 mod p,
    Norm(f)=p*lambda^2*tau.                                        (2)

If `p` does not divide `tau`, this is a certificate with `p`
occurring to exponent exactly one, of the kind used in
[pass39](parallel39-shell-structural-input-2026-09-06.md).
The condition `V<1` is sufficient for `p` not dividing `tau`,
by the previously established inequality `tau<=p V^(N/2)`.
No condition `V<1` is needed for the integral formula (2) itself.

**Proof.** Since `h` is real, it vanishes at both `u` and `u^(-1)`.
Its complex multiplication determinant has `p`-adic valuation two.
Consequently these are exactly its zero split-root evaluations;
any additional zero would force its determinant to be divisible
by `p^3`.

The reduction of `bar(F)` is supported only at `u^(-1)`. Hence
`h bar(F)` vanishes at every split root. It is zero in `R/pR`,
proving coefficientwise divisibility by `p` and integrality of `f`.

To identify the zero of `f` after this division, define

    c=lambda*p/h in R^+.

This is integral by the real multiplication adjugate identity;
the possible negative sign of `Norm_+(h)` has no effect. Its
complex norm is `(lambda*p)^(N-2)`. The identity `hc=lambda*p`
shows that `c` vanishes at every split root other than `u` and
`u^(-1)`. It is nonzero at both of those two roots: an additional
zero would force a higher `p`-adic valuation than `N-2` in its
norm. Now

    fc=lambda*bar(F).

At `u`, the right side is zero and `c(u)!=0`. Thus `f(u)=0`.
Finally,

    Norm(f)=Norm(h) Norm(F)/p^N
           =(lambda*p)^2*p^(N-1)*tau/p^N
           =p*lambda^2*tau.

This proves (2) without a class-group assumption or a unit-norm
assumption. The additional arithmetic information is the explicit
real certificate (1).

For example, define `k_min(u)` as the minimum norm cofactor over
all nonzero `f in R` satisfying `f(u)=0` and `v_p(Norm(f))=1`.
For every admissible `h` and every `F` with `p` not dividing `tau`,

    k_min(u)<=lambda^2*tau.                                        (3)

If the real evaluation kernel is principal, a generator has
`|Norm_+(h)|=p`, so one may take `lambda=1`. This concerns the
real kernel; it does not require the full complex kernel to be
principal. No uniform bound on the best `lambda` is proved here.

## 3. When lambda is one, the lattice correspondence is exact

Assume now that (1) holds with `lambda=1`, so `c=p/h` is integral.
Define the two full-rank integer lattices

    I_u={F in R:F(v)=0 mod p for every split root v!=u},
    P_u={f in R:f(u)=0 mod p}.

The maps

    Phi(F)=h*bar(F)/p,
    Psi(f)=c*bar(f)=p*bar(f)/h                                      (4)

are mutually inverse additive bijections `I_u <-> P_u`.
The first integrality and image assertion are the argument above.
For the second, `c` is zero modulo `p` outside `u,u^(-1)`, and
`bar(f)` is zero at `u^(-1)`; their product belongs to `I_u`.
Composing the formulas gives the identity because `h` and `c`
are real and `hc=p`.

The norm relation is exactly

    Norm(Phi(F))/p = Norm(F)/p^(N-1).                               (5)

Thus minimization of the two normalized norms agrees over nonzero
vectors of these **unrestricted integral lattices**, including if one restricts to
norm cofactors coprime to `p`. This supplies a much stronger
converse than integer scalar reciprocals when the real certificate
is available.

The actual coset representatives form a particular section of
`I_u/pR`, obtained by centering the coefficients. The map (4) sends
`pR` to `hR`, so it also identifies `I_u/pR` with `P_u/hR`.
However, coefficient centering is a minimization of Euclidean
length, not of the cyclotomic norm. Equation (5) does not assert
that a norm-minimizing lattice vector lies in that centered
section. The explicit failure below is relevant to any attempt
to turn (4) into a cofactor minimization algorithm.

## 4. A real generator for the known order-64 example

Take `p=6700417,n=64,N=32,D=16,g=2`. In the real basis of Section 2,
use the coefficient vector

    h: [1,0,1,0,-1,-1,1,0,-1,0,-1,0,-1,1,1,0].

The exact coefficient certificate verifies

    bar(h)=h,
    Norm_+(h)=p,
    Norm(h)=p^2,
    h(2)=h(2^(-1))=0 mod p.

This vector was discovered in one bounded reduction of the
16-dimensional real evaluation-kernel lattice. Its claimed
properties are checked independently of that numerical discovery
step: a separate integer 16-by-16 multiplication determinant gives
the real norm, and integer quadratic norm descent in degree 32
gives the complex norm. No shortest-vector claim is made.

For the actual `a=1` polynomial from pass40,

    V=6120237/p<1,       Norm(F)=p^31*1217,

the recovered `f=h bar(F)/p` has coefficient vector

    [-1,0,0,-1,0,0,1,0,1,0,1,0,0,0,0,-1,
     0,0,-1,1,0,1,0,0,0,-1,-1,0,0,-1,1,0].

It satisfies `f(2^(-1))=0` and `Norm(f)=1217*p`. This replaces
the scalar-reciprocal cofactor `1217^31` by 1217. It does not
improve the already known element `2-X^(-1)` of norm `641*p`,
and no minimum cofactor across all elements is claimed.

A second finite check uses `h_large=5-2s_1`. Its real norm is
`641*p`, and (2) gives exactly the cofactor `641^2*1217`.
Thus the general `lambda^2` factor is checked as well as the
principal-real-kernel case.

## 5. Centering can greatly increase the norm cofactor

Start instead from the known smaller-cofactor element

    f_0=2-X^(-1),       Norm(f_0)=641*p.

Its inverse image under (4) is the integral vector

    F_0=c*bar(f_0),       Norm(F_0)=p^31*641.

Let `F_1` be its coefficientwise centered reduction modulo `p`.
It is an actual coset polynomial, with `a=2922708` and generator 2.
Exact computation gives

    sum_j (F_0)_j^2/p^2 =56196475/p,
    V(a)=sum_j (F_1)_j^2/p^2=20588347/p,
    Norm(F_1)=p^31*178771841.

Centering strictly decreases the squared radius, as it must, but
increases the norm defect from 641 to 178771841. The centered
radius here is greater than one; this example tests norm behavior
under the canonical-section operation, not a new `V<1` assertion.

The decrease in radius and increase in norm are both exact integer
or rational facts. They disprove the specific proposed shortcut
that centering in (4) preserves or lowers the cofactor being
minimized. They do not prove that the global minimum over all
centered cosets exceeds the unrestricted lattice minimum.

## 6. Verification and the remaining uniform obstruction

The stdlib-only
[verifier](../experiments/parallel41_shell_real_recovery.py)
checks the displayed coefficients, real and complex determinants,
all split-root supports, both recovery formulas, norm transfers,
the `lambda=641` case, and the centering counterexample. It writes
the [exact certificate](../results/parallel41_shell_real_recovery_2026_09_06.json).
Run `python3 experiments/parallel41_shell_real_recovery.py`.
The finite check makes no coset scan or numerical trigonometric
estimate.

This gives a quantitative arithmetic recovery unavailable to the
scalar reciprocal method: a real norm cofactor `lambda` turns
the defect `tau` into an actual cofactor `lambda^2*tau`, and
`lambda=1` gives an exact correspondence of unrestricted norm
minimization problems. But neither a uniformly small `lambda`
nor a strong enough relation between `V` and `tau` has been
established here. The old bound `tau<=p V^(N/2)` still permits
an exponential defect even well below the dimension-scale mean.
The explicit centering failure also prevents replacing norm
minimization by shortest-coordinate representatives. These are
the additional arithmetic losses that remain before this route
could control the required uniform shell deviations.
