# A projective transfer for all three degree-two necklace labels

**Status: partial, uniform results; Paley and the full localization conjecture remain unproved.** Inversion transfers every previously proved binary singleton estimate to each of the other two choices of a two-label alphabet. In particular, all 189 length-six words using at most two of the three possible nonempty labels are now bounded. This adds 125 words to the 64 binary singleton words already covered. A separate exact identity handles the first word using all three labels and completes the degree-two necklace estimate through length three.

These are ordinary mathematical proofs. No novelty or formalization claim is made. The finite checks below supplement the proofs and do not prove the remaining uniform conjectures.

## Definitions and inputs

Let p be prime, p≡1 mod 4, and let χ be its quadratic character, with χ(0)=0. Write

\[
 S_{xy}=\chi(x-y),\qquad
 d_A(x)=\chi(x),\quad d_B(x)=\chi(x-1),\quad
 d_C(x)=\chi(x)\chi(x-1).
\]

Thus A, B, C denote the labels {0}, {1}, {0,1}, respectively. Let D_L be the diagonal matrix with entries d_L(x), set C_L=D_L S, and define

\[
 N(w)=\operatorname{tr}(C_{w_1}\cdots C_{w_k}).
\]

The monochromatic quantities are t_j=tr(C_A^j)/(p−1), j≥1. The imported analytic input is

\[
 |t_j|\le (j-1)p^{(j-1)/2},
 \tag{1}
\]

derived in [the preceding identity note](localized-necklace-identities.md) from [Lu–Zheng–Zheng, Lemma 2.1 (2.3)](https://arxiv.org/html/1305.3405v3#S2). The necklace convention is that of [Kunisky, Definition 1.13](https://arxiv.org/html/2303.16475v1). The previously proved binary formulas and their precise scope are in [the planar reduction note](planar-necklace-reductions.md).

## Uniform transfer under inversion

Let π fix A and exchange B,C. For any word w, let n_L(w) count its L positions. Then

\[
 \boxed{\bigl|N(w)-N(\pi w)\bigr|
       \le (n_B(w)+n_C(w))p^{k/2}.}
 \tag{2}
\]

This estimate holds at every length k, including lengths growing with p. No hidden constant depends on k.

First restrict every necklace variable to F_p*. Denote the result by N*(w). In the substitution x_i=1/y_i, all variables remain nonzero, and

\[
 d_A(1/y)=d_A(y),\qquad
 d_B(1/y)=d_C(y),\qquad
 d_C(1/y)=d_B(y).
\]

Moreover,

\[
 \chi(1/y_i-1/y_{i+1})
 =\chi(y_i-y_{i+1})\chi(y_i)\chi(y_{i+1}).
\]

Here χ(−1)=1 is used. Multiplying around the cycle cancels every extra χ(y_i)^2. Thus the exact identity is

\[
 \boxed{N^*(w)=N^*(\pi w).}
 \tag{3}
\]

Equivalently, on F_p* let R be the permutation matrix for inversion, let S* and D_L* be the restricted matrices, and let T=diag(χ(x)). Then R²=T²=I and

\[
 R D_L^* R=D_{\pi L}^*,\qquad
 R S^*R=T S^*T,\qquad
 R(D_L^*S^*)R=T(D_{\pi L}^*S^*)T.
\]

Taking the trace of a product gives (3), since the successive T factors cancel. No inverse or permutation is defined at zero in this argument.

It remains to restore the zero coordinate, without silently dropping its contribution. Set P=I−e_0e_0^T. Then

\[
 N^*(w)=\operatorname{tr}(P C_{w_1}\cdots P C_{w_k}).
\]

The elementary relation S²=pI−J gives ||S||=√p. Hence ||C_L||≤√p and ||P C_L||≤√p for each L. Replace C_{w_i} by P C_{w_i} successively. For the difference at position i, put A_j=P C_{w_j} for j<i. The exact rank-one term is

\[
 \operatorname{tr}(A_1\cdots A_{i-1}e_0e_0^T C_{w_i}\cdots C_{w_k})
 =e_0^T C_{w_i}\cdots C_{w_k}A_1\cdots A_{i-1}e_0.
\]

Since ||e_0||=1, its magnitude is at most the product of the k operator norms, namely p^(k/2). If w_i is A or C, the contribution is exactly zero, since their diagonal entries at zero vanish. Therefore

\[
 |N(w)-N^*(w)|\le n_B(w)p^{k/2}.
 \tag{4}
\]

For πw the corresponding count is n_C(w). Equations (3) and (4) prove (2). This argument uses a telescoping sum of at most k rank-one errors, rather than an exponentially large inclusion-exclusion bound.

Reflection x↦1−x exchanges A,B and fixes C exactly, without a restriction or error term. Together, reflection and inversion supply the two transpositions generating all permutations of the three labels. For a general permutation one can compose (2) and the exact reflection, retaining the corresponding error terms. They do not remove a label from a word that uses all three.

## Consequences for two-label words

Every word on {A,C} is the image under π of a binary word on {A,B}; reflection gives the analogous statement for {B,C}. Thus every proved bound for a binary word transfers with an additive error at most k p^(k/2).

In particular, for any two distinct labels L,M in {A,B,C}, r,s≥1, and k=r+s,

\[
 \boxed{|N(L^rM^s)|\le(k^2+k)p^{k/2}.}
 \tag{5}
\]

The same bound holds for a nonconstant word using two labels when one of them occurs at most twice. These follow from the two-block and one-/two-occurrence bounds already proved for binary singleton words. Their normalized bound is (k²+k)/p, so it tends to zero uniformly for k=o(√p).

For constant words the transfer gives instead

\[
 |N(L^k)|\le (k-1)(p-1)p^{(k-1)/2}+kp^{k/2}.
 \tag{6}
\]

The second term is unnecessary for the original singleton labels. It is retained as a convenient uniform expression for all three labels.

At length six, the previously established binary bound 36p^(7/2) and (2) give

\[
 \boxed{|N(w)|\le42p^{7/2}
   \quad\text{if }w\in\{A,B,C\}^{6},\ |\operatorname{im}(w)|\le2.}
 \tag{7}
\]

There are 3·2^6−3=189 such words. The bound is deliberately loose. The other 540 length-six words use all three labels and are not covered by this result.

## The first genuinely three-label word

There is also an exact identity

\[
 \boxed{N(ABC)=t_4.}
 \tag{8}
\]

To prove it, allow distinct anchors a,b and put A at a triangle vertex x, B at y, and C at z. Let N_{a,b} be this triangle sum. Under x=a+(b−a)u, and similarly for y,z, its seven character factors give

\[
 N_{a,b}=\chi(b-a)N(ABC).
\]

Consider the graph H on vertices a,b,x,y,z with edges

\[
 xy,yz,zx,\ ax,az,\ by,bz,\ ab.
\]

The last edge weights the anchor sum by χ(a−b). Its zero value removes a=b. Consequently its affine partition function is

\[
 Z_p(H)=\sum_{a,b}\chi(a-b)N_{a,b}
       =p(p-1)N(ABC).
 \tag{9}
\]

The same graph is a wheel: its hub is z and its four rim vertices, in cyclic order, are x,y,b,a. Fix the hub by translation and sum the rim variables. This gives

\[
 Z_p(H)=p\operatorname{tr}(C_A^4)=p(p-1)t_4.
 \tag{10}
\]

Equating (9) and (10) proves (8). Every ordering of three distinct labels is a rotation or reversal of ABC, and necklace reversal is valid because S is symmetric and all D_L are diagonal. Hence (8) covers all six orderings. By (1), their absolute values are at most 3p^(3/2).

All shorter words, and all other length-three words, use at most two labels. The preceding transfer and binary estimates therefore prove the degree-two necklace asymptotic for every word of lengths at most three. This is not a statement for all fixed lengths, nor for larger anchor sets.

## Verification and remaining gap

[The exact script](../experiments/parallel_necklace_2026_09_04.py) constructs the full Legendre matrices using Python integers. It checks the restricted inversion identity, the rank-one error estimates, all 189 length-six two-label words, longer transferred two-block examples, and all six orderings in (8). In small fields, separate literal sums over the original coordinate tuples check the full and restricted definitions. A literal wheel graph sum independently checks (9)–(10), without substituting the claimed necklace value. [The result file](../results/parallel_necklace_2026_09_04.json) records the check counts and input hashes.

The transfer controls the zero-coordinate corrections uniformly. It supplies no cancellation for a general binary word whose binary bound is still unknown, and it does not turn a word using all three labels into a two-label word. General growing-depth necklaces and the extreme-eigenvalue conclusion remain unproved. None of these statements gives the arbitrary-small-set Paley bound, the thin dyadic subgroup bound, or a bridge to the prize benchmark.
