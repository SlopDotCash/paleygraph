# Single-size moment tests and an additive-energy obstruction

**Status: a conditional reduction and an unconditional counterexample to
a stronger proposed bridge. Neither (LM) nor Paley is proved or refuted.**
The reduction shows precisely how moment information at one cardinality
controls larger, possibly very unequal rectangles. The construction shows
that even the exact minimum of every lower additive energy does not imply
a Gaussian high character moment. Its sets have polynomial size for
every sufficiently large prime; there is no least-prime assumption.
No novelty or formalization claim is made.

Let p be an odd prime, chi(0)=0,
F_B(x)=sum_(b in B) chi(x-b), and M_(2r)(B)=sum_x F_B(x)^(2r).
Write S(A,B)=sum_(a in A) F_B(a).

## 1. An exact transfer from one cardinality

For integers 1<=k<=n=|B| and r>=1, let C be a uniformly chosen
k-element subset of B. Pointwise, E F_C=(k/n)F_B. Convexity therefore
gives the exact inequality

\[
 \boxed{M_{2r}(B)\le(n/k)^{2r}
       \frac1{\binom nk}\sum_{C\in\binom Bk}M_{2r}(C).}       \tag{1}
\]

In particular, suppose the following input is available at **only the
single size k**, uniformly in C:

\[
 M_{2r}(C)\le L\qquad(|C|=k).                              \tag{2}
\]

For every A and B with |B|>=k, Holder and (1) imply

\[
 \boxed{\frac{|S(A,B)|}{|A||B|}
       \le\left(\frac{L}{|A|k^{2r}}\right)^{1/(2r)}.}       \tag{3}
\]

The exact average in (1) can replace L in (3). There is no partition
remainder or rounding loss in these statements. Since
S(B,A)=chi(-1)S(A,B), either set can be the one being sampled.
Consequently (2) controls every rectangle with max(|A|,|B|)>=k,
using min(|A|,|B|) in the denominator of (3).

This elementary reduction is useful here because (LM) demands estimates
at all cardinalities. The sufficient inputs below have a smaller stated
domain and tolerate an additional power of p. They remain unproved;
the reduction is not an independent bound on an arbitrary rectangle.

## 2. A sparse sequence of single-size estimates would suffice

Let r_j>=2 be an unbounded sequence of integers. For each j fix
beta_j>=0, with beta_j tending to zero, and put
k_j(p)=floor(p^(1/(r_j+1))). Suppose, for every j, some C_j satisfies

\[
 \boxed{M_{2r_j}(C)\le C_j p^{1+\beta_j}k_j(p)^{r_j}
          \quad\hbox{whenever }|C|=k_j(p)}                 \tag{SS}
\]

for all sufficiently large odd primes p. Then the full classical Paley
two-set conjecture follows.

Indeed, fix epsilon>0 and choose a single j for which
theta_j+beta_j<epsilon/2, where theta_j=1/(r_j+1).
For m,n>p^epsilon, one has n>=k_j and, eventually,
k_j>=p^theta_j/2. Equation (3) gives

\[
 \left(\frac{|S(A,B)|}{mn}\right)^{2r_j}
 \le C_j2^{r_j}p^{\theta_j+\beta_j-\varepsilon}.
\]

Thus delta=epsilon/(8r_j)>0 works for sufficiently large p. All constants
and the prime threshold may depend on the chosen j and epsilon.
No uniformity in a growing r=r(p) is assumed or inferred.

For a single fixed r, (SS) gives the more precise rectangular range

\[
 \max(|A|,|B|)\ge k(p),\qquad
 \min(|A|,|B|)\ge p^{1/(r+1)+\beta+\eta}
 \quad\Longrightarrow\quad
 \frac{|S(A,B)|}{|A||B|}\ll_{r,C}p^{-\eta/(2r)}.            \tag{4}
\]

For example, the unproved estimate M_4(C)<=C p|C|^2 only at
|C|=floor(p^(1/3)) would establish Paley for every epsilon>1/3.
It would not establish the remaining smaller exponents. In the exact
fourth-moment decomposition M_4=R_C+D_C from the preceding passes,
0<=R_C<=3p|C|^2 at these sizes. Hence a one-sided upper bound
D_C<=C p|C|^2 on that slice would suffice for this particular range.
This last observation identifies the required signed estimate; it
does not prove it.

The standard termwise Weil bound gives (SS) with a fixed-r loss

\[
 \beta_r=\frac{r-1}{2(r+1)},\qquad
 \frac1{r+1}+\beta_r=\frac12.                              \tag{5}
\]

This computes exactly why the known termwise input still stops at
one half after this reduction. To cover every epsilon by (SS), the
allowed loss must tend to zero along the selected ranks. The required
Gaussian estimate with beta=0 is not the false all-cardinality Gaussian
claim: its sole set size lies strictly below p^(1/r), where the forced
row obstruction overwhelms a constant Gaussian bound.

There is also a more permissive alternative slice. If for an unbounded
sequence of fixed r one can prove

\[
 M_{2r}(C)\le C_r p^{2+\alpha_r},\qquad
 |C|=\lfloor p^{1/r}\rfloor,\qquad \alpha_r>0,
 \quad\alpha_r\longrightarrow0,                           \tag{CS}
\]

then (3) again proves full Paley: choose r with both 1/r<epsilon and
alpha_r<epsilon/2. Fixed positive alpha_r can absorb a logarithm; the
choice alpha_r=0 with a constant C_r is false, as shown by the existing
forced-row argument and strengthened within structured sets below.
Neither (SS) nor (CS) has been shown equivalent to Paley. No converse
or strict logical separation is claimed. (LM) would imply these
sufficient inputs in their indicated ranges.

## 3. Polynomial-size sets with minimal lower additive energies

A set Q in F_p is a B_h set if equality of two h-term sums from Q,
with repetitions allowed, forces equality of the two multisets.
It is also B_j for every 1<=j<=h, by padding with a fixed element.
Such constructions are classical; see
[Bose and Chowla, Theorems in the additive theory of numbers](https://doi.org/10.1007/BF02566968).
The following digit construction is included in full, so no size claim
about an imported construction is needed.

Fix h>=2. For every sufficiently large p, choose a prime q>h with
q=p^(1/h+o(1)) and

\[
 2h^{h+1}q^h\le p.
\]

Put X=(p/(2h^(h+1)))^(1/h) and m=floor(X/2). Bertrand's postulate
gives a prime m<q<2m<=X. For sufficiently large p, this q is greater
than h and satisfies q>X/2-1, establishing the stated size.
There is no congruence condition on q and no relation between q and p
other than the displayed size bound. Set b=hq and

\[
 a(t)=\sum_{j=1}^h (t^j\bmod q)b^{j-1},\qquad
 Q=\{a(t):t\in\mathbb F_q\}\subseteq\mathbb F_p.           \tag{6}
\]

The first digit is t, so these q integers are distinct. Every a(t)<b^h
and h b^h<=p/2. Equality of two h-term sums modulo p is therefore an
integer equality. No base-b carries occur, since each digit sum is at
most h(q-1)<b. The two multisets consequently have equal power sums
modulo q through degree h. Newton's identities, with 1,...,h invertible
in F_q, give the same monic polynomial with the selected t values as
roots. Thus their multisets agree, and Q is B_h in F_p.

For any subset B of Q, n=|B|, define the usual additive energy

\[
 E_j^+(B)=\#\{(b_1,\ldots,b_j,c_1,\ldots,c_j)\in B^{2j}:
                   \sum b_i=\sum c_i\}.
\]

For 1<=j<=h it equals the exact permutation-only minimum

\[
 E_j^+(B)=(j!)^2[z^j]
          \left(\sum_{a=0}^j\frac{z^a}{(a!)^2}\right)^n.
                                                               \tag{7}
\]

In particular E_2^+(B)=2n^2-n, and E_3^+(B)=6n^3-9n^2+4n
when h>=3. These are additive exponential-sum energies, not the
quadratic-character moments M_(2r).

## 4. The structured obstruction, with all prime-size quantifiers

**Theorem.** Fix integers r>h>=2 and a real 0<alpha<1/h. For every
sufficiently large odd prime p there is B subset F_p with
n=floor(p^alpha), satisfying (7) for every j<=h, but

\[
 \boxed{M_{2r}(B)\ge
       \left(\frac{1/h-\alpha}{\log2}+o(1)\right)
                   n^{2r}\log p.}                         \tag{8}
\]

Here the o(1) depends only on fixed h and alpha. The assertion also
holds for any positive integer moment index in place of r; the
condition r>h selects the range where it obstructs a Gaussian bound.

Use Q from (6), of size q=p^(1/h+o(1)), and put d=(p-1)/2. For a
k-element row set U define
Q(U)={b in Q: chi(a-b)=1 for every a in U}. Each column b has exactly
d positive rows, irrespective of its membership in Q. Double counting
therefore gives

\[
 \sum_{U\in\binom{\mathbb F_p}k}|Q(U)|=q\binom dk,
 \qquad \max_U|Q(U)|\ge q\frac{\binom dk}{\binom pk}.       \tag{9}
\]

Take k=floor(log_2(q/n))-1. For sufficiently large p it is positive,
k=O(log p), and

\[
 \frac{\binom dk}{\binom pk}
 \ge 2^{-k}\left(1-\frac{k^2}{p-k+1}\right)\ge2^{-k-1}.
\]

The right side of (9) is at least n. Select n elements B of such Q(U).
Then B remains B_h, and F_B(a)=n for all a in U. In particular no
chosen row lies in B: chi(0)=0 has been excluded by the definition.
There are
k=((1/h-alpha)/log2+o(1))log p such rows, which proves (8).

Taking 1/r<alpha<1/h gives

\[
 \frac{M_{2r}(B)}{p n^r}\ge
    \left(\frac{1/h-\alpha}{\log2}+o(1)\right)
                    p^{\alpha r-1}\log p\longrightarrow\infty.
\]

Taking alpha=1/r still makes this ratio diverge logarithmically.
Thus the implication

“all additive energies through h are permutation-only, therefore
M_(2r)(B)<=C_(r,h) p|B|^r”

is false for every r>h>=2, including h=r-1. This is the new scope of
the obstruction: it survives exact minimality of every lower additive
energy, on polynomial-size sets, for every sufficiently large prime.
It is not another assertion of prime existence with no polynomial
control, nor an additional lower-bound constant for arbitrary sets.

It does **not** refute (LM), (SS), (CS) with positive alpha_r, or the
size-restricted Fourier-norm candidate. In (8), the logarithmic term
allowed by (LM) accommodates the forced rows. The hypothesis (SS)
lies on a smaller size slice; at that slice, n^r log(p)/p tends to zero.
An upper-bound argument using additional character-sensitive information,
or allowing the necessary logarithmic spike term, is not excluded.

## Verification and remaining work

The companion verifier checks (1) by exact subset averaging and checks
(3) using actual maxima on small complete fields. It independently
checks the digit containers and the permutation-only energies, constructs
forced-row witnesses inside them, and checks the moment lower bounds
without floating-point arithmetic. Rational exponent calculations check
the two sufficient slices and the exact half-threshold in (5).
These finite checks do not establish (SS) or (CS). The uniform
construction and reduction are proved above; the missing upper bound
has not been replaced with numerical plausibility.

The completed bounded run contains 12,384 subset-average checks,
92,410 rectangle-transfer checks, eight structured witnesses, 38
additive-energy checks, 24 moment lower-bound checks, and 99 fixed-r
exponent audits. The eight witnesses include B_2 and B_3 containers;
they are finite checks of the general construction, not its asymptotic
justification. The exact counts and final note/script hashes are in
[the result](../results/parallel5_classical_2026_09_04.json).
