# A short shell need not come from a principal evaluation kernel

There is an actual quartic-window counterexample to the proposed
converse `V(a)<=1` implies a norm-`p` generator. At

    p=6700417, n=64, H=<2>, a=1,

the shell is `V(1)=6120237/6700417<1`, but the corresponding
evaluation kernel in `Z[X]/(X^32+1)` is nonprincipal. The proof uses
the factorization `2^32+1=641*p` and an elementary Gauss-sum upper
bound at the auxiliary prime 641. It needs no class-group calculation.

An exact converse with a larger cofactor is available: if the actual
coset polynomial has norm `p^(N-1)*tau` and `V<1`, then it produces
a norm-`p*tau^(N-1)` element in the evaluation kernel. In the example,
`tau=1217`, and the exponent `N-1=31` is unavoidable within this
scalar-reciprocal construction. Thus a short shell does not give the
small cofactor needed to reverse
[pass39](parallel39-shell-structural-input-2026-09-06.md).

The uniform thin-annulus estimate and the full goals remain open. These
are ordinary mathematical arguments and exact finite checks; no build,
proof-status change, or literature novelty claim is made.

## 1. What an arbitrary actual coset does provide

Use the notation of
[pass36](parallel36-shell-inversion-2026-09-05.md): `n=2N>=4` is a
power of two, `H=<g>` has order `n` modulo the prime `p`, and

    F(X)=sum_(j=0)^(N-1) r_j X^j in R=Z[X]/(X^N+1),

where `r_j` is the centered integer congruent to `a g^j`, with `a!=0`.
Its reduction has a unique nonzero evaluation at `u=g^(-1)`:
`F(u)=aN`, while all other split-root evaluations vanish. Let

    P_u={h in R:h(u)=0 mod p},
    tau=Norm(F)/p^(N-1).

The earlier rank-one argument shows that `tau` is a positive integer.
We use its already proved norm inequality in the form

    tau <= p V(a)^(N/2).                                           (1)

The following construction is an algebraic converse; it does not
improve that arithmetic-geometric-mean inequality.

**Reciprocal converse.** For every such `F`,

    T=tau*p/F belongs to R,
    T(u)=0 mod p,
    Norm(T)=p*tau^(N-1).                                           (2)

If `V(a)<1`, then `tau<p`, so `p` occurs in `Norm(T)` to exponent
exactly one. Thus (2) is a certificate of the type required in pass39,
but with cofactor `tau^(N-1)`.

**Proof.** Let `M_F` be the integer multiplication matrix of `F`.
Modulo `p` it has rank one. Every `(N-1)`-row minor of `M_F` is
divisible by `p^(N-2)`: integer row operations can leave all but one
row divisible by `p`, and each such minor contains at least `N-2`
of those rows. Inverting the unimodular row operations preserves
this divisibility of minors. Consequently

    adj(M_F)/p^(N-2)

is an integer matrix. Applied to the coefficient vector of one it
represents `Norm(F)/(p^(N-2)F)=tau*p/F`, proving integrality.
Evaluating `FT=tau*p` at `u` gives `T(u)=0`, because `F(u)!=0`.
Multiplicativity of the norm gives

    Norm(T)=(tau*p)^N/(p^(N-1)*tau)=p*tau^(N-1).

Finally (1) and `V<1` give `tau<p`. This proves the assertions.

In particular `tau=1` forces `P_u` to be principal: (2) then gives
an element `T in P_u` of norm `p`, so the two subgroups `TR` and
`P_u` of `R` have the same index `p` and are equal. A useful precise
hypothesis ensuring this case is

    V(a)<(2/p)^(2/N),

because then (1) implies `tau<2`. Conversely, if `P_u` is
nonprincipal, the known norm lower bound sharpens to

    V(a)>=(2/p)^(2/N).                                             (3)

The gain in (3) is only a factor `2^(2/N)` over the previous norm
bound. It is not a dimension-scale lower bound and does not approach
the required thin-annulus estimate.

## 2. A short actual coset whose kernel is nonprincipal

Take `p=6700417`, `n=64`, `N=32`, and `g=2`. The exact centered
vector at `a=1` is

    [1,2,4,8,16,32,64,128,
     256,512,1024,2048,4096,8192,16384,32768,
     65536,131072,262144,524288,1048576,2097152,-2506113,1688191,
     -3324035,52347,104694,209388,418776,837552,1675104,3350208].

It satisfies

    sum r_j^2=41008140038829,
    V(1)=6120237/6700417<1,
    Norm(F)=p^31*1217.                                             (4)

The prime is in the working window `64^4/4<=p<=64^4`.
This is the pre-existing order-64 subgroup from
[the finite spectral obstruction](finite-spectral-obstruction.md);
the small-shell calculation and the nonprincipality argument below
are the relevant new consequences. Searches for this shell value,
the norm defect 1217, and the proposed reciprocal converse found no
earlier research note establishing them.

We first rule out a norm-641 generator at the auxiliary prime

    q=641,   H_q=<2>,   |H_q|=64,   [F_q^*:H_q]=10.

Suppose `v in R` has `Norm(v)=641` and `v(2)=0 mod 641`.
Pass39 would give a nonprincipal-frequency period with

    eta >=64-4 pi^2 >64-1936/49=1200/49.                            (5)

On the other hand, an elementary character calculation gives

    M(H_q)<=(1+9 sqrt(641))/10 <(1+9*26)/10=47/2.                   (6)

The exact margin is `1200/49-47/2=97/98>0`, contradicting (5).
Thus no such `v` exists.

For completeness, (6) uses only the usual finite character
orthogonalities. If a subgroup has index `m`, its indicator on
`F_q^*` is `m^(-1) sum_chi chi`, summed over the `m` multiplicative
characters trivial on the subgroup. The principal character
contributes `-1` to a nonzero additive frequency. Every other
character contributes a Gauss sum of magnitude `sqrt(q)`, so
the triangle inequality gives `(1+(m-1)sqrt(q))/m`. To check the
magnitude directly, for a nontrivial character `chi` substitute
`x=ty` in the squared Gauss sum. The inner additive sum over
`y!=0` is `q-1` when `t=1` and `-1` otherwise. Since
`sum_(t!=0)chi(t)=0`, the squared magnitude is `q`.

Now suppose the kernel `P_2` at the original prime `p` is principal,
say `P_2=fR`. Its generator has norm `p`. The elementary identity

    Norm(2-X)=2^32+1=641*p

and `2-X in P_2=fR` show that

    v=(2-X)/f in R,      Norm(v)=641.

Since the multiplication determinant of `f` is `p`, it is invertible
modulo 641. In particular its evaluation at 2 modulo 641 is
nonzero. Evaluating `fv=2-X` there forces `v(2)=0`, contradicting
the preceding paragraph. Therefore `P_2` is nonprincipal.

Conjugation `X -> X^(-1)` sends `P_2` to `P_(2^(-1))`, so the
kernel `P_u` corresponding to the actual polynomial `F` in (4) is
also nonprincipal. Indeed all evaluation kernels above `p` are
conjugate under `X -> X^j`, `j` odd, and none can be principal.

This is a counterexample to `V<=1` implying a norm-`p` generator
in the actual quartic family. It is not a counterexample to pass39,
which proved the opposite implication.

## 3. The reciprocal loss is exact in this example

In (4), `tau=1217` is prime and different from `p`. The converse
(2) gives `T=1217*p/F` with the coefficient vector

    [109,-467,-556,-331,-11,-232,-357,547,
     805,-613,237,-312,-606,274,579,-69,
     -17,84,-773,-261,359,159,145,71,
     -494,-351,231,4,195,683,133,-514].

Direct integer multiplication gives `FT=1217*p`, and

    Norm(T)=p*1217^31.

This is the least possible cofactor among scalar reciprocals
`T_m=m*p/F`, with `m` a positive integer. If `T_m` is integral,
its integer norm must be

    Norm(T_m)=p*m^32/1217.

Primality of 1217 and `1217!=p` therefore imply `1217|m`.
The choice `m=1217` is integral by (2), so it is the exact minimum.
Its norm cofactor is the 96-digit integer `1217^31`.

This is a limitation of scalar reciprocal recovery, not a lower
bound for every possible norm-cofactor certificate: the already
available `2-X^(-1)` has norm `641*p` and vanishes at the same
root `u`. It illustrates how much algebraic information is lost by
recovering a kernel element through `m*p/F` alone.

For an additional exact consistency check, the sparse carry
polynomial `C` defined by

    (2-X^(-1))F=pC

has nonzero coefficients `C_21=1`, `C_22=-1`, `C_23=1`,
`C_24=-1`, and `C_31=1`. Its norm is

    Norm(C)=641*1217=780097,

which recovers the norm factor in (4) independently through this
multiplicative identity. No minimal-cofactor assertion is made for
the value 641.

## 4. Flat Gram data require a separate converse

The example above is not flat: `F bar(F)!=p^2`, and its full
negacyclic Gram product has 30 nonzero nonconstant coefficients.
The scalar inequality `V<1` must not be substituted for a flat
Gram identity.

Here is the precise ideal-factorization interpretation when an
integral `A` really does satisfy `A bar(A)=p^2` and its reduction
has a unique nonzero evaluation at `u`. This paragraph uses the
usual prime-ideal factorization in the cyclotomic integer ring;
none of it is needed for the elementary counterexample above.
Writing `P_bar=P_(u^(-1))`, one has

    (A)=(p) P_bar/P_u.                                             (7)

Indeed the valuation at `P_u` is zero. The identity `A bar(A)=p^2`
forces valuation two at `P_bar` and one at every remaining split
prime. There are no prime factors above any other rational prime.
Thus (7) implies the class equality `[P_u]=[P_bar]`. It does not
by itself identify their common class with the principal class.
Any such additional assertion needs a result about the classes
fixed by conjugation.

Conversely, if `P_bar/P_u` is principal, choose a generator `beta`.
Then `epsilon=beta bar(beta)` is a totally positive unit fixed by
conjugation. An integral flat representative with the ideal (7)
exists exactly when

    epsilon^(-1)=e bar(e) for some unit e in R.

In that case take `A=p beta e`. This is an actual norm-of-unit
condition; it cannot be replaced by an assumption that every
totally positive unit is a square. This formulation makes no
claim that nontrivial fixed classes or an obstructing unit occur
in any particular dyadic field.

Even the field-level conjugate-quotient equation alone does not
recover a norm-`p` generator. For a flat `A`, an explicit integral
solution is

    f_0=p+bar(A),       bar(f_0)/f_0=A/p.                            (8)

It is nonzero because `A=-p` would contradict the nonzero reduction
hypothesis. But modulo `p`, `f_0=bar(A)` is nonzero only at
`u^(-1)`: it vanishes at the other `N-1` split roots. Therefore
`p^(N-1)` divides its norm. All field solutions of the quotient
equation differ by a nonzero scalar fixed by conjugation. Finding
such a scalar that removes the extra integral factors while
retaining a small norm is an additional arithmetic problem, not
an automatic consequence of the quotient equation.

## 5. Verification and the resulting restriction

The stdlib-only
[verifier](../experiments/parallel40_shell_converse.py) checks the
three relevant primes, both generator orders, the quartic window,
the exact vector and split-root values, norm descent with every
integer inverse identity, `FT=1217*p`, the sparse carry identity,
the nonflat Gram product, and all rational inequalities at 641.
Its [certificate](../results/parallel40_shell_converse_2026_09_06.json)
contains the exact coefficient arrays and norm factors. Run
`python3 experiments/parallel40_shell_converse.py` from the root.
No coset census or numerical trigonometry is used.

The next shell argument cannot infer a small principal norm
cofactor from `V<=1`: the actual order-64 example rules that out.
Reciprocal recovery is available with a precisely quantified but
potentially enormous cofactor, and the flat Gram case retains
class and unit-norm conditions. Controlling those losses or finding
a different ingredient remains necessary before this converse
could contribute to the uniform thin-annulus estimate.
