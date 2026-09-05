# Positive sixth-aggregate review: an exact quartic-star interface

**Status: no new uniform upper bound for T₆ is proved.** The original
Sidon sixth-moment input remains open. This review obtains an exact
character-sensitive reformulation using only squarefree quartic sums
away from C. It also identifies why applying the ordinary Weil bound
to its off-diagonal terms reproduces the original obstruction.

The [pass-21 squarefree reduction](parallel21-squarefree-moments-2026-09-05.md)
and [pass-20 transfer](parallel20-inversion-moments-2026-09-05.md) were
read directly. No prior file or central assessment was changed by this
review. No novelty claim is made for the elementary identities below.

## 1. The remaining estimate can be stated entirely with quartic sums

Let p be an odd prime, χ the quadratic character with χ(0)=0, and
C⊂Fₚ have n elements. Write

    e_j(x)=Σ_(Q⊂C, |Q|=j) ∏_(c∈Q)χ(x−c),
    T_j=Σ_x e_j(x),
    K(A)=Σ_x ∏_(a∈A)χ(x−a),

where A can be a multiset and K(∅)=p. Define the quartic-star sum

    U_C(t)=Σ_(Q⊂C, |Q|=3) K(Q ⊔ {t})
          =Σ_x χ(x−t)e_3(x),
    L_out(C)=Σ_(t∉C) U_C(t)².                          (1)

Every summand K(Q⊔{t}) appearing in L_out has four distinct roots.
The affine curve y²=∏_(a∈Q⊔{t})(x−a) has a smooth projective
model of genus one; the estimate below only uses its ordinary
multiplicative-character Weil bound and does not require a trace
normalization at infinity.

**Proposition.** Suppose n≥6 and n⁴≤p. Then

    L_out(C) ≤ 20p T_6(C)+p²n³,                        (2)
    20p T_6(C) ≤ L_out(C)+3p²n³.                       (3)

Thus, on this size range,

    T_6(C)≤A pn³   ⇒   L_out(C)≤(20A+1)p²n³,
    L_out(C)≤D p²n³   ⇒   T_6(C)≤((D+3)/20)pn³.       (4)

The constants are deliberately loose. The proposition holds for all
C; it does not assume Sidon. In particular, on the precise requested
slice n=⌊p^(1/4)⌋, a uniform quartic-star energy bound for Sidon C
would give the missing positive T₆ bound. **The new statement (1) is
a reformulation of the missing cancellation, not its proof.**

### Proof of the transform identity

For x,y∈Fₚ,

    Σ_t χ(x−t)χ(y−t)=p·1_(x=y)−1.                      (5)

The diagonal case counts p−1 nonzero squares. If x≠y, the substitution
u=(x−t)/(y−t), excluding t=y, runs through Fₚ\{1}; the summand is
χ(u), whose sum is −1. Consequently the actual, full p-column
translate kernel gives

    L_all:=Σ_t U_C(t)²
       =p Σ_x e_3(x)²−T_3².                           (6)

This is valid for both congruence classes of odd primes. The transpose
in the definition (1) avoids any χ(−1) convention issue. Also
Σ_t U_C(t)=0. Identity (5) concerns every field coordinate x and every
translate t. The Gram matrix of just n abstract columns, even if it
is pI−J, is not a replacement for this full kernel identity.

### Exact overlap ledger and the zero rows

For N nonzero signs, counting pairs of unordered triples by their
intersection size gives

    e_3²=20e_6+6(N−4)e_4+(N−2)(N−3)e_2+binom(N,3).     (7)

For disjoint triples the coefficient is binom(6,3)=20. For intersection
size one, choose its element outside the four remaining indices and
then split these four into an ordered pair of two-element sets, giving
6(N−4). Intersection size two gives (N−2)(N−3); identical triples
give binom(N,3).

At a character row x=c∈C there are n−1 nonzero signs. An especially
simple correction follows by replacing its zero coordinate by an
independent sign σ: the new cubic coefficient is e_3+σe_2, so its
squared average is e_3²+e_2². The average of each elementary symmetric
coefficient on the right of (7) is unchanged. Thus, with

    Z=Σ_(c∈C)e_2(c)² ≥0,
    R=6(n−4)T_4−(n−2)(n−3)binom(n,2)+p binom(n,3)−Z,

we obtain the exact identity

    Σ_x e_3(x)²=20T_6+R.                               (8)

The pair correlation T₂=−binom(n,2) was used here. No exceptional row
has been discarded or treated as having n nonzero signs.

For four distinct roots, |K|≤3√p. Therefore

    |6(n−4)T_4| ≤ (3/4)√p n⁵ ≤ (3/4)pn³,
    0≤(n−2)(n−3)binom(n,2)+Z
       ≤ n⁴/2+n⁵/4 ≤ pn³/4,
    0≤p binom(n,3)≤pn³/6.                              (9)

For the middle inequality, |e₂(c)|≤binom(n−1,2)≤n²/2;
the last step follows from n²+2n≤n⁴≤p for n≥2. The first uses
n²≤√p. Keeping track of the signs of the three terms of R now gives

    −pn³≤R≤(11/12)pn³≤pn³.                             (10)

For the degree-three sums, the ordinary Weil estimate also gives

    T_3²≤4p binom(n,3)²≤pn⁶/9≤p²n³/9.                 (11)

The last step uses n³≤p. The only imported character-sum estimate is
the usual distinct-root bound |K(A)|≤(|A|−1)√p, stated as
[Lemma 2.1 in McDonald–Sahay–Wyman](https://arxiv.org/html/2210.03789v2#S2).
That primary HTML statement was checked directly; no theorem from its
VC-dimension analysis is imported.

### The n omitted parameters already have the required energy scale

If t∈C and Q contains t, write Q={t,a,b}. The repeated root gives

    K(Q⊔{t})=−1−χ((t−a)(t−b)),

since the complete distinct-pair correlation is −1 and the row x=t
has to be subtracted. It follows exactly that

    U_C(t)=Σ_(Q⊂C\{t}, |Q|=3) K(Q⊔{t})
                 −binom(n−1,2)−e_2(t).                (12)

In particular,

    |U_C(t)|≤3√p binom(n−1,3)+2binom(n−1,2)
              ≤(1/2)√p n³+n²≤√p n³,
    L_in:=Σ_(t∈C)U_C(t)²≤pn⁷≤p²n³.                   (13)

The last inequality is exactly where n⁴≤p controls the omitted
parameter set. All quartic sums in the remaining energy L_out are
squarefree. Summing (12) also supplies a useful independent check:

    Σ_(t∈C)U_C(t)=4T_4−n binom(n−1,2)−Σ_(t∈C)e_2(t).   (14)

Combine (6) and (8) to obtain

    L_out+L_in=20pT_6+pR−T_3².                         (15)

Using R≤pn³ and the nonnegativity of L_in,T₃² proves (2). Solving
for 20pT₆ and using (10)–(13) gives

    20pT_6≤L_out+(19/9)p²n³≤L_out+3p²n³,

which is (3). The same identities also give T₆≥−pn³/20; this is a
lower bound and does not address the required positive upper side.

## 2. Why a naive elliptic large-sieve step would be circular

For a triple Q⊂C write f_Q(x)=∏_(q∈Q)χ(x−q) and
U_Q(t)=K(Q⊔{t}). The full translate identity gives, exactly,

    Σ_t U_Q(t)U_R(t)
      =p[ K(Q△R)−Σ_(a∈Q∩R)∏_(d∈Q△R)χ(a−d) ]
           −K(Q)K(R).                                 (16)

For Q=R the bracket is p−3. For distinct triples sharing two
indices it contains a degree-two sum and two explicitly subtracted
values. If they share one index it contains a degree-four sum and
one subtracted value. **For disjoint triples it is exactly**

    Σ_t U_Q(t)U_R(t)=p K(Q∪R)−K(Q)K(R).                (17)

Thus, despite each U_Q(t) being a quartic character sum in x, its
pairwise correlation over t contains the original degree-six sum.
An asserted near-orthogonality bound for this family cannot be
justified by applying the individual quartic Weil bound. Applying
Weil to (17) and taking absolute values gives the old scale
√p n⁶ for the aggregate T₆. At n≈p^(1/4), this is a factor
p^(1/4) larger than pn³. It produces no new exponent.

The Sidon hypothesis is not used in (5)–(17). A future proof must
supply an actual estimate on the signed aggregate of these
character-sensitive correlations, using more than the bare fact
that equal two-term sums in C are trivial. Pass 20 already explains
why this restricted hypothesis, if strong enough, transfers to all
sets at the same size. The weighted sign obstructions in pass 21
are not actual instances of (5), so they neither establish nor
refute the missing estimate (1).

## 3. Exact audit and scope

The [verifier](../experiments/parallel21_positive_upper_review_2026_09_05.py)
uses exact integer character tables, direct distinct-subset products,
and a complete integer translate transform. It checks (6), (8),
(12), (14), and every overlap type in (16), together with the explicit
constants of (2)–(3) at several genuine Sidon sets on the requested
size slice. It explicitly checks integer overflow bounds; no
floating-point threshold decides acceptance.

The exhaustive part covers every C of size at least three in F₅
and F₇. Additional small fields test varied sizes and both prime
congruence classes. The larger examples are actual prime-field
kernels for p=1297, 2411, 4099, 6563, 10007, with
n=6,7,8,9,10 respectively. These examples validate finite identities
and constants, not a uniform upper bound as p tends to infinity.
The exact counts, set labels, moment values, and input hashes are
saved in the [results](../results/parallel21_positive_upper_review_2026_09_05.json).

The bounded review did not prove T₆=O(pn³), the Paley conjecture,
the subgroup square-root conjecture, or a bridge to the official
Reed–Solomon prize. Its concrete output is the identity (15), the
controlled on-C contribution, and the precise unproved off-C
quartic-star energy estimate (1).
