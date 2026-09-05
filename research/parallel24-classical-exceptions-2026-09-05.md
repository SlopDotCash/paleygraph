# Exceptional sets: exact local drift and the remaining quantitative loss

**Status: no upper bound for every exceptional Sidon set was obtained.**
This bounded investigation proves exact local identities, measures the
family forced by one natural persistence argument, and quantifies the
loss in removing the signs of a good unsigned inverse. It does not claim
another typical-set theorem or a new worst-case Paley exponent. No
abstract countermodel is introduced.

The [pass-23 theorem](parallel23-classical-upper-2026-09-05.md) and its
[independent review](parallel23-classical-independent-review-2026-09-05.md)
leave a proportion O(n²/p) of the Sidon slice uncontrolled. The
[pass-20 transfer](parallel20-inversion-moments-2026-09-05.md) still needs
a bound on every input of that slice. The analysis below addresses
local persistence and the actual transported signs, without reproducing
the separate root lane's unsigned inverse-orbit averaging calculation.

Throughout, χ(0)=0 is the quadratic character of an odd prime field,
F_C(x)=Σ_(c∈C)χ(x−c), and M_j(C)=Σ_xF_C(x)^j for even j.
For the quantitative comparison, n=|C|≥6 and n⁴≤p<(n+1)⁴.
Normalize an input by R=M_6(C)/(pn³).

## 1. Exact one-swap drift, including the zero entries

Choose a uniformly from C and b uniformly outside C, and set
C′=C\{a}∪{b}. Write d=n(p−n). For a fixed field row x, put

    F=F_C(x),       N=n−1_(x∈C),
    S=n(p−1)+(p−2n)N.

Then the exact row formula is

    d E[F_(C′)(x)^6]
      =[d−6(p−5)]F^6
       +(15S−80p+180)F^4
       +(15S+90N(p−1−N)−96p+122)F²
       +S+30N(p−1−N).                                  (1)

In particular, the coefficient of M_6 in the complete swap average is

    λ_6=1−6(p−5)/[n(p−n)].                             (2)

Here is a direct derivation. Let u=χ(x−a), v=χ(x−b), and Δ=v−u.
The two selections are independent, conditional on C. The population
inside C has sum F and N nonzero entries; the outside population has
sum −F and p−1−N nonzero entries. Binomial expansion therefore gives

    d EΔ   =−pF,
    d EΔ²  =S+2F²,
    d EΔ³  =−(4p−3)F,
    d EΔ⁴  =S+8F²+6N(p−1−N),
    d EΔ⁵  =−(16p−15)F,
    d EΔ⁶  =S+32F²+30N(p−1−N).

Substitution into (F+Δ)^6 proves (1). This retains the difference
between N=n and N=n−1 rather than treating every row as fully nonzero.
It works in both congruence classes of odd primes.

For n≥6 and p≥n⁴, all three lower-degree coefficients in (1) are
nonnegative. For example, using N≥n−1 and p≥2n, the fourth-power
coefficient is at least

    p(30n−95)−30n²+15n+180
      ≥30n²−175n+180>0.

The second-power coefficient is at least

    p(120n−201)−120n²+105n+122
      ≥120n²−297n+122>0.

The constant coefficient is nonnegative directly. Thus

    E M_6(C′)≥λ_6 M_6(C).                              (3)

For n≥7, λ_6≥0 and this lower-persistence bound can be iterated.
At n=6 the exact coefficient is slightly negative; no useful positive
iteration is asserted there. As n grows, λ_6=1−6/n+O(1/p).
This is a local mean-persistence statement, not a uniform upper bound.

## 2. A full-character local variance inequality still contains local high moments

For a fixed (n−1)-set A, let F=F_A and

    H=6F⁵+20F³+6F,       Z=F⁴+F²,
    X_A=Σ_xH(x)²,       Y_A=Σ_xZ(x)².

The actual transform Qh(b)=Σ_xχ(x−b)h(x) satisfies

    ||Qh||₂²=p||h||₂²−(Σ_xh(x))².                      (4)

The exact insertion formula from pass 23 is

    M_6(A∪{b})=M_6(A)+15M_4(A)+15M_2(A)+p−1
                          +(QH)(b)−15Z(b),  b∉A.

Consequently, with m=p−n+1,

    Var_(b∉A) M_6(A∪{b})
       ≤ [√(pX_A)+15√Y_A]²/m.                         (5)

This follows from (4), the triangle inequality, and bounding variance
by the complete squared sum divided by m. It applies to an individual
A, with no random-input hypothesis. But X_A contains 36M_10(A), and
Y_A contains M_8(A). Pass 23 controlled their averages over all A;
that average cannot be applied to an adversarial A obtained by deleting
a point of an exceptional C. Replacing them by elementary inequalities
such as M_10(A)≤|A|⁴M_6(A) supplies no new upper bound for that input.
Thus the localization is exact, but the missing local estimate remains
explicit in (5).

## 3. The forced local family is much too small for the variance contradiction

Choose D uniformly among n-sets satisfying |D∩C|=t. Every point of C
is selected with probability t/n and every point outside C with
probability (n−t)/(p−n). Since the complete character row sums to zero,

    E[F_D(x)]=λ_tF_C(x),
    λ_t=(pt−n²)/[n(p−n)].                              (6)

Convexity, now on each actual character row, gives

    E[M_6(D)]≥λ_t^6 M_6(C).                            (7)

The ordinary individual Weil bound also gives, for every n-set,

    M_6(D)≤15pn³+5√p n⁶≤(15+5n)pn³.                   (8)

Indeed at most 15n³ index words have all multiplicities even and
contribute at most p. Each remaining monic degree-six polynomial is
not a square and its complete character sum has absolute value at
most 5√p. This is the standard Weil input already used in the earlier
squarefree notes, not a new estimate for its signed aggregate.

Let θ_t be the proportion in this shell with M_6(D)>15pn³. If
λ_t^6R>15, (7) and (8) imply

    θ_t ≥ (λ_t^6R−15)/(5n).                            (9)

However, 0≤λ_t≤t/n when t≥1, and R≤15+5n. For the right side to
be positive, it is therefore necessary that

    t>n[15/(15+5n)]^(1/6)
       =n[3/(n+3)]^(1/6)>n^(5/6),     n≥6.             (10)

The shell t=0 has negative λ_0, so (7) uses its absolute sixth power.
It cannot force this threshold: p−n≥n³ gives
|λ_0|^6R≤(15+5n)/n^12<15. Thus even the largest exception
allowed by the existing uniform bound must retain more than n^(5/6)
of its original labels before this mean argument forces any bad sets.

For a globally uniform n-set D, write I=|C∩D|. The exact factorial
moment and a union bound give

    P[I≥t]≤E binom(I,t)
      =binom(n,t)(n)_t/(p)_t
      ≤[n²/(p−n)]^t≤(2n²/p)^t≤(2/n²)^t.              (11)

Therefore **even declaring every set retaining the necessary points
bad** would produce a family of relative size at most

    (2n²/p)^(n^(5/6))≤(2/n²)^(n^(5/6)),              (12)

which is smaller than every fixed power of 1/n. At this slice the
pass-23 exceptional allowance is of order n²/p≈n^(−2). The local
lower bound (9), multiplied by its shell volume, is far below that
allowance. The p-dependent form in (12) should be retained if the
upper restriction p<(n+1)⁴ is dropped: the coarse last expression
alone does not give a uniform comparison with n²/p when n is fixed
and p grows. This particular propagation-plus-variance argument cannot
make the exceptional probability zero.

This calculation already uses the unconditional all-set variance bound.
Requiring the propagated family also to be Sidon cannot enlarge its
global measure. No claim that local shells consist of Sidon sets has
been made. The loss above describes the evaluated argument; it is not
a theorem forbidding a stronger arithmetic propagation mechanism.

## 4. Transported signs: balanced splitting returns the old power

Suppose z∉C and D={(c−z)^(-1):c∈C}. The actual transported weights
are

    w_d=χ(c−z)=χ(d),
    M_6(D,w)=M_6(C)−|F_C(z)|⁶+n⁶≥M_6(C).              (13)

This is the already-proved pass-20 identity. Suppose, as an additional
input, this particular unsigned inverse satisfies

    M_6(D)≤Kpn³                                         (14)

for some fixed K. This note does not establish (14) or repeat the root
lane's orbit-wise averaging proof.

Choose a majority sign η and let E={d∈D:w_d≠η}, with a=|E|≤n/2.
The pointwise identity is exactly

    F_(D,w)=η(F_D−2F_E),
    a=(n−|F_C(z)|)/2.                                  (15)

If R>K, the triangle inequality and (13)–(14) force a>0 and

    M_6(E)/(pa³)
       ≥(n/a)³ [R^(1/6)−K^(1/6)]⁶/64.                (16)

At a balanced split a=n/2, the lower bound is asymptotically R/8
when R/K tends to infinity. It does not force a stronger normalized
exception at the smaller size. A gain would require enough imbalance,
for example a<n/4 in that asymptotic comparison; no such pole selection
has been proved.

One can see the lack of a power saving directly. Applying the same
individual Weil bound to E gives

    R^(1/6)
      ≤K^(1/6)+2[15(a/n)³+5a⁶/(√p n³)]^(1/6)
      ≤K^(1/6)+2[15/8+5n/64]^(1/6).                  (17)

For fixed K, the sixth power of the last expression has leading term
5n, exactly the old leading power in (8). This particular use of a
good unsigned inverse, a sign split, and the termwise bound supplies
no improvement in that power.

Strongly unbalanced poles are sparse by the actual second moment:

    #{z∉C: a<n/4}
       ≤4M_2(C)/n²=4(p−n)/n.                          (18)

Indeed a<n/4 means |F_C(z)|>n/2, and M_2(C)=pn−n². A positive-density
supply of poles satisfying an unsigned moment bound therefore need not
meet this O(p/n)-sized set. This is the missing joint-selection issue;
cardinality alone supplies no intersection. The sign constraint is
actual χ(d), not a freely chosen or synthetic weighting.

## 5. Checks and the remaining obligation

The [exact verifier](../experiments/parallel24_classical_exceptions_2026_09_05.py)
checks complete swap averages, shell means, local full-transform
variance bounds, the actual inverse weights, both zero-row corrections
in (13) and (15), and the volume comparison. Its
[results](../results/parallel24_classical_exceptions_2026_09_05.json)
record:

- 1,745 complete swap identities and 20,839 row-coefficient checks;
- 86 exact shell-mean and sixth-moment inequalities;
- 24 local full-transform variance checks;
- 144 actual transported-sign identities and 24 unbalanced-pole checks;
- 7 thin-slice volume and coefficient-positivity certificates.

The actual character calculations use the small primes 5,7,11,13.
They verify exact identities, including cases outside the thin slice;
the nonnegative-coefficient and volume assertions are checked separately
on actual primes with n=6,8,16,32,64,128,256. Those larger cases compute
integer combinatorial certificates only, not their character matrices.
No floating-point threshold is used.

The calculations did not force useful new additive structure on the
high-value rows of an exceptional Sidon set. Nor did they select an
unsigned-good inverse with sufficiently unbalanced transported signs.
The exact missing steps are: local control of the higher moments in
(5), a propagation family much larger than (12), or a joint pole
estimate that overcomes (18). No one of these inputs is asserted.
The uniform classical upper bound and full goal remain unproved.
