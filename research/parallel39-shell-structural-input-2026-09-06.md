# Principal evaluation kernels force a small shell coset

The uniform thin-annulus estimate remains open. This note gives an
additional structural constraint: an element of norm `kp` in a split
evaluation kernel, with `p` occurring to exponent one, forces an actual
nonzero coset with `V<=k^2`. In particular a principal evaluation kernel
forces `V<=1`. The known norm-prime trinomial at `p=215535361, n=128`
supplies an exact witness with `V=1` and refutes the literal amplitude
constant `C=2` at this finite quartic-window parameter.

This does not refute a bound with an unspecified absolute constant or an
asymptotic bound allowing finitely many exceptions. Existing sublinear
period bounds already prevent norm-prime templates from persisting to
arbitrarily large orders in the quartic window. The new implication is
proved below by ordinary algebra; no Lean build or proof-status change
was made, and no literature novelty is claimed.

## 1. Scope and prior-work check

The notation and target are those of
[pass36](parallel36-shell-inversion-2026-09-05.md) and
[the subgroup target](subgroup-target.md). Thus `n=2N` is a power of two,
`n>=4`, `p` is prime with `n|(p-1)`, and `H=<g>` has order `n`.
For centered representatives of `aH/{+1,-1}`,

    V(a) = sum_(j=0)^(N-1) r_j(a)^2/p^2,
    vbar = n(p+1)/(24p),
    eta(a) = sum_(h in H) exp(2 pi i ah/p).

The earlier shell argument bounds `V(a)` below by `p^(-2/N)` using a
rank-one cyclotomic norm. Searches of the research notes for shell,
negacyclic, autocorrelation, principal-kernel, conjugate-quotient and
the finite witness terms did not locate the implication below.

The trinomial and its norm-prime property are **previous results**:
[pass30](parallel30-triangle-remainder-2026-09-05.md),
[pass31](parallel31-short-multiples-2026-09-05.md), and
[pass32](parallel32-linear-triangle-length-2026-09-05.md) use them to
analyze short relations and triangle derivations. This note uses their
conjugate quotient to obtain a shell and period witness.

## 2. A norm cofactor controls an actual shell coset

Let

    R = Z[X]/(X^N+1),        bar(F)(X) = F(X^(-1)).

Suppose a nonzero `f in R` satisfies

    f(g)=0 mod p,     Norm(f)=kp,     k>=1,     p does not divide k.    (1)

Here `Norm` is the determinant of integer multiplication by `f` in
the coefficient basis. Then there is a nonzero `a in F_p` such that

    V(a) <= k^2,
    eta(a) >= n-4 pi^2 k^2.                                         (2)

**Proof.** The polynomial `X^N+1` is irreducible over `Q` and splits
over `F_p` into the distinct roots `g^u`, `u` odd modulo `n`.
Multiplication by `f` modulo `p` is diagonal in the evaluation basis,
with those evaluation values as eigenvalues. Its determinant has
`p`-adic valuation one. Consequently exactly one eigenvalue vanishes:
two zero eigenvalues would make the integer determinant divisible by
`p^2`, by integer elimination or Smith normal form. The unique root
where `f` vanishes is `g`.

Put `d=kp` and `B=d/f in R`. Integrality follows from the integer
adjugate identity; equivalently,

    B(X) = product_(u odd mod n, u != 1) f(X^u) in R.

The full product including `u=1` is `Norm(f)=d`. At `X=g`, every factor
in the displayed product is nonzero modulo `p`, so `B(g)!=0`.
At every other split root `t`, the identity `f(t)B(t)=d` and
`f(t)!=0` give `B(t)=0` modulo `p`.

Now define the integral polynomial

    A = bar(f) B = d bar(f)/f.                                      (3)

Because `g!=g^(-1)` when `n>=4`, `bar(f)(g)=f(g^(-1))!=0`. Hence
`A(g)!=0`, and `A(t)=0` at all other split roots. Finite Fourier
inversion gives, coefficient by coefficient,

    A_j = a g^(-j) mod p,     a=A(g)/N != 0.                        (4)

For a direct check of this inversion, the polynomial
`(1/N) sum_(j=0)^(N-1) g^(-j)X^j` takes value one at `g` and zero
at every other split root, by a finite geometric sum.

Conjugating (3) gives the exact identity

    A bar(A) = d^2 in R.                                            (5)

The constant coefficient of `A bar(A)` is `sum_j A_j^2`.
Thus `sum_j A_j^2=d^2`. Center each coefficient independently modulo
`p`, writing the result as `r_j`. Centering cannot increase its absolute
value. Equation (4) shows that these are the representatives of
`a g^(-j)`, and `g^(-1)` generates the same `H`. Therefore

    V(a)=sum_j r_j^2/p^2 <= d^2/p^2 = k^2.

Finally `cos(t)>=1-t^2/2` implies

    eta(a)=2 sum_j cos(2 pi r_j/p)
          >= n-4 pi^2 V(a) >= n-4 pi^2 k^2.

This proves (2). In particular the construction produces a genuine
subgroup coset; it is not a synthetic vector with selected moments.

If the index-`p` evaluation kernel `P_g={F in R:F(g)=0}` is principal,
say `P_g=fR`, its generator has norm `p`: multiplication by `f` has
image of index `p`, and its determinant is positive because the
complex embeddings occur in conjugate pairs. Hence (1) applies with
`k=1`.

## 3. What the higher norm and trace identities do, and do not, enforce

Every complex embedding of `A` in (3) has absolute value `d`. Thus
`Norm(A)=d^N`; its conjugate norms are perfectly flat. If `M_A` is its
negacyclic multiplication matrix, (5) says

    M_A M_A^T = d^2 I_N.

All normalized positive Gram traces equal `N`, at every positive
integer power. For `k=1`, these especially rigid norm and trace
identities coexist with a shell of squared radius at most one,
whereas the exact average shell radius is about `N/12`.

This supplies a concrete warning for the proposed additional trace
approach: flat conjugate magnitudes and vanishing nontrivial
negacyclic autocorrelations do not themselves force a coset near its
dimension-scale mean. The arithmetic availability of `f` matters.
This finite example does not establish that those constraints have
such realizations in an unbounded quartic-window family.

The implication has a useful quantitative reverse consequence. If a
valid uniform period upper bound is `M<=U(n,p)`, then every `f`
satisfying (1) obeys

    k^2 >= max(0, (n-U(n,p))/(4 pi^2)).                              (6)

Wherever `U(n,p)=o(n)`, this gives `k>=(1-o(1))sqrt(n)/(2 pi)`.
For example, the already documented ordinary argument in
[pass28](parallel28-moment-recurrence-2026-09-05.md) gives
`U << n^(71/72)(log n)^(1/18)=o(n)` in `n^4/4<=p<=n^4`.
That use imports the literature inputs and verification scope of
pass28; the elementary implication (6) does not need that particular
upper bound. Its absolute constant is not converted here into a
finite size threshold.

Consequently norm-prime or, more generally, `k=o(sqrt(n))` templates
cannot persist to arbitrarily large `n` in that window. Searching for
an unbounded principal quartic-window family would conflict with
already available sublinear cancellation, so it is not a viable
route to refuting the unspecified-constant shell target.

The desired stronger shell estimate `max|V-vbar|<=D` would further
force `k^2>=vbar-D` for every such certificate. That is only a
necessary consequence, not a converse or a new proof of the target.

## 4. Exact finite witness and the literal constant two

Use the previously identified parameters

    p=215535361, n=128, N=64, g=25525303,
    f=1+X+X^19,       Norm(f)=p.

They satisfy `n^4/4=67108864<=p<=268435456=n^4`. Reconstructing
the exact `B` with `fB=p` and setting `A=bar(f)B` gives

    a=38468180,
    A_j = centered_residue(a g^(-j) mod p),
    sum_j A_j^2 = p^2 = 46455491841400321,
    A bar(A)=p^2,
    V(a)=1,
    vbar=3448565792/646606083,
    max|V-vbar| >= 2801959709/646606083.                            (7)

Every coefficient of `A` is already centered, so the bound in (2)
is an equality for `V`. The exact arrays `f,B,A`, all root evaluations,
and rational comparisons are in the
[coefficient certificate](../results/parallel39_principal_shell_2026_09_06.json).
The scalar is `a=A(g)/N=A_0 mod p`, not `B(g)/N`.

No numerical trigonometry is needed for the following finite
obstruction. Since `pi<22/7`, (7) gives

    eta(a) >= 128-4 pi^2 > 4336/49.                                (8)

For a self-contained way to obtain the strict bound on `pi`, integrate
`x^4(1-x)^4/(1+x^2)` from zero to one. Polynomial division gives
its value `22/7-pi`, and the integrand is strictly positive inside
the interval.

Also `e>1+1+1/2+1/6=8/3`, and the exact integer comparison is

    128*8^15 - 215535361*3^15 = 1410902777170069 > 0.

Hence `log(p/n)<15`. Finally,

    (4336/49)^2 - 4*128*15 = 361216/2401 > 0.

Together these prove

    eta(a) > 2 sqrt(n log(p/n)).                                    (9)

Thus the literal `C=2` version fails at this finite point. The
pre-existing order-64 example refuted the literal `C=sqrt(2)`;
(9) is a stronger finite constant obstruction. It supplies no
uniform lower growth rate and does not refute an unspecified
absolute constant.

## 5. Verification and the remaining action

Run `python3 experiments/parallel39_principal_shell.py` from the
repository root. The stdlib-only
[verifier](../experiments/parallel39_principal_shell.py) reconstructs
the adjugate by exact quadratic norm descent, checks every lifted
identity, and checks primality, generator order, all 64 centered
coordinates, the full negacyclic product, all 64 split-root values,
and the integer/rational comparisons in (7)--(9). It performs no
scan over `F_p` or its quotient cosets and evaluates no approximate
cosines. The finite verification is separate from the ordinary
uniform implication (2).

The useful new restriction is (6): a small cofactor of a principal
multiple forces a highly coherent coset, so any sublinear period
bound already forbids cofactors below a constant times `sqrt(n)`
asymptotically. Future norm/trace work must retain that arithmetic
cofactor obstruction instead of assuming that stronger Gram
regularity automatically gives the mean shell. Neither (6), the
finite constant obstruction, nor the existing norm lower bound
controls the required deviations of all other cosets at the
`sqrt(n log(p/n))` scale.
