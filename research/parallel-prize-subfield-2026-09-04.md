# Subfield-valued words restrict MCA challenges

**Status: an ordinary mathematical reduction, not a Paley proof, a prize
threshold, or a Lean certificate.** This pass makes one obstruction to a
Paley-to-prize reduction more precise. No novelty is claimed.

Let E/F be a finite-field extension, with |F|=p and |E|=q, and let the
Reed–Solomon evaluation domain D lie in F. The code consists of evaluations
of polynomials of degree less than k. The same statements apply to any
number of interleaved rows, with column agreement as the metric.
An MCA-bad challenge is defined by an agreement set for f+gamma g which
has no simultaneous degree-<k explanations for f and g on that same set.
This is the event audited in [the prize reduction note](prize-reduction-audit.md),
using [ABF's definition](https://eprint.iacr.org/2026/680).

## 1. Linear maps of the alphabet preserve the code

Every F-linear map L:E^s -> E^t preserves polynomial degree when applied
coefficientwise. Because x belongs to F at every evaluation point,

\[
 L\!\left(\sum_{i=0}^{k-1}a_i x^i\right)
 =\sum_{i=0}^{k-1}L(a_i)x^i.
\]

Thus applying L to a codeword gives a codeword in t rows. In particular,
an E-valued interleaved code on D, viewed in an F-basis of E, is exactly
an F-valued code with [E:F] times as many rows and the same column distance.
Multiplication by an extension challenge becomes its F-linear matrix on
these components; this observation does not reduce it to scalar F-folding.

## 2. Exact descent of bad challenges for base-field words

Suppose f,g take values in F^s. Then at every radius,

\[
 \boxed{\operatorname{Bad}_{E}(f,g)
       =\operatorname{Bad}_{F}(f,g)\subseteq F.}
 \tag{1}
\]

**Proof.** For gamma outside F, the vectors 1,gamma are F-linearly
independent. Choose F-linear functionals lambda,nu:E -> F satisfying
lambda(1)=1, lambda(gamma)=0 and nu(1)=0, nu(gamma)=1.
If an E-codeword h agrees with f+gamma g on S, applying lambda and nu
rowwise gives F-codewords agreeing with f and g respectively on S.
Thus gamma is not bad.

For gamma in F, choose an F-linear map pi:E -> F with pi(1)=1.
It fixes all base-field entries. An E-codeword explaining f+gamma g
on S projects to an F-codeword explaining it there. Similarly, a pair
of E-codewords explaining f,g projects to a pair of F-codewords.
The converse in both cases is inclusion of F-codewords into E-codewords.
Therefore both requirements in the MCA event agree exactly. This includes
gamma=0 and imposes no minimum size on S beyond the event's definition. QED.

Consequently, maximizing only over such received pairs gives exactly

\[
 \sup_{f,g\in(F^s)^D}\frac{|\operatorname{Bad}_E(f,g)|}{q}
 =\frac{p}{q}\,\varepsilon_{\rm mca}(RS[F,D,k]^{\bowtie s},\delta)
 \le\frac{p}{q}.
 \tag{2}
\]

The same upper bound applies if f and g belong to scalar multiples of
the base-field alphabet, or do so after subtracting separate codewords.
Adding codewords does not change the bad event; nonzero scalar factors
rescale its challenge set. This is a restriction on a class of received
words, not an upper bound over all E-valued words.

## 3. A ratio-set bound for arbitrary small alphabet subspaces

More generally suppose f(x) lies in an F-linear subspace U of E^s and
g(x) lies in an F-linear subspace V of E^s for every x in D. Define

\[
 R(U,V)=\{\gamma\in E^*:U\cap\gamma V\ne\{0\}\}.
\]

Then

\[
 \boxed{\operatorname{Bad}_E(f,g)\subseteq\{0\}\cup R(U,V).}
 \tag{3}
\]

**Proof.** If gamma is nonzero and U and gamma V intersect trivially,
define F-linear maps on U direct-sum gamma V by
L_f(u+gamma v)=u and L_g(u+gamma v)=v, and extend them to E^s.
Applying these maps to a folded explaining codeword gives simultaneous
explanations for f,g, by the degree-preservation lemma. The argument does
not apply at gamma=0, which is why zero is retained in (3). QED.

If a=dim_F U and b=dim_F V, then

\[
 |R(U,V)|\le
 \min\!\left(q-1,\frac{(p^a-1)(p^b-1)}{p-1}\right).
 \tag{4}
\]

Indeed each gamma in R has a pair (u,v) of nonzero vectors with
u=gamma v. A fixed pair determines at most one gamma, even for s>1.
Simultaneously scaling u,v by F* gives p-1 distinct pairs for that gamma.
These pair sets are disjoint for different gamma. This proves (4), also
when a or b is zero. Equations (3)-(4) give the uniform bound

\[
 \frac{|\operatorname{Bad}_E(f,g)|}{q}
 \le\min\!\left(1,
 \frac{1+(p^a-1)(p^b-1)/(p-1)}q\right).
 \tag{5}
\]

Neither this counting bound nor (3) asserts that every permitted challenge
is bad. When U=V=F in one row, R=F*, recovering the bound p/q.

## 4. What this says about the pinned prize profile

The [pinned profile](official-profile-and-trace.md) has
p=2130706433 and q=p^6. Therefore every base-field-valued received pair
has MCA error at most p^-5, and

\[
 p^5>2^{128},\qquad p^{-5}<2^{-128}.
\]

This extends the earlier bound for one monomial pair to **all**
base-field-valued pairs, at every radius. A worst-case MCA argument must
handle more general extension-valued words. Ordinary base-field Paley
cancellation by itself does not yet supply that argument.

There is no contradiction with the [large winning-set construction](list-to-winning-set.md).
That construction uses (w,0) and exploits the additional inner-product
claim. Its MCA-bad set is empty, since every folded agreement witness has
the simultaneous explanation (u,0). The protocol winning set and the MCA
bad set are different quantities.

The elementary code/projection proof is independent of the deep necklace
input and of any unproved subgroup bound. Its finite tests are in
[the verifier](../experiments/parallel_prize_subfield_2026_09_04.py) and
[the result file](../results/parallel_prize_subfield_2026_09_04.json).

The verifier exhausts all 729 base-field received pairs of length three
over F_3, and compares their exact bad sets in F_3 and F_9 at three
agreement thresholds (2187 comparisons). It separately enumerates all
25 base-field and 625 extension-field degree-at-most-one codewords on
four points in F_5 and checks 64 sampled received pairs at four thresholds
(256 comparisons). In F_81 it checks two dimension-two alphabet spaces:
their nonzero ratio set has exactly 32 elements, attaining (4), and 96
sampled received pairs pass all three threshold checks. Every challenge
and every possible agreement subset is exhausted for each checked pair.
All arithmetic is integral, and the three extension moduli are checked by
verifying that every nonzero element of each quotient is a unit.
