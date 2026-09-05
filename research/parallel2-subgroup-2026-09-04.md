# Primitive-root fibers and the balanced mixed energy

**Status: the uniform subgroup bound remains unproved.** This pass gives
an exact polynomial test for the balanced mixed energy, an upper bound in
terms of its squarefree defect, and a class with the quadratic bound
B≤2k². The defining condition of that class is itself an unproved
collision restriction in the full quartic family. A new finite quartic
resonance of parent order 512 is certified below. No spectral bound,
infinite family in the growing quartic window, novelty claim, or Lean
formalization follows.

## A polynomial containing only the primitive half

Let n=2k≥4 be a power of two. Write H=μ_n⊂F_p*, where p≡1 mod n,
K=μ_k and L=H\K. Every element of L has exact order n. The balanced
energy is

\[
 B=\sum_x(1_K*1_L)(x)^2.
\]

Let ζ be a primitive complex n-th root. Define the degree-d polynomial,
where d=n/4=k/2,

\[
 R_n(Y)=\prod_{\substack{1\le j<n/2\\j\text{ odd}}}
               (Y-(1+\zeta^j)^n). \tag{1}
\]

These indices select one element from each inverse pair of primitive
roots. The expression is independent of the choice of ζ and has integer
coefficients: the Galois group permutes the roots, which are algebraic
integers. They are distinct over C, since

\[
 (1+\zeta^j)^n=-\bigl(2\cos(\pi j/n)\bigr)^n
\]

has strictly decreasing absolute value in the displayed range. The
Galois action is transitive, so R_n is irreducible over Q, though its
reduction need not be squarefree. Its constant term is 2^k.

For example,

\[
 R_4(Y)=Y+4,\qquad R_8(Y)=Y^2+136Y+16.
\]

Its exact power sums admit an elementary coefficient formula:

\[
 s_r=\sum_{R_n(\lambda)=0}\lambda^r
 =\frac k2\sum_{a=0}^{2r}(-1)^a\binom{nr}{ka},\qquad r\ge1. \tag{2}
\]

To see this, sum (1+h)^(nr) over all n-th roots and subtract the sum
over all k-th roots, leaving precisely the primitive n-th roots. The
root-of-unity filter keeps exponents divisible by n in the first sum
and by k in the second. Divide by two for inverse pairs. Newton's
identities then construct R_n using integers only.

Compared with the earlier P_n in
[kernel-discriminant.md](kernel-discriminant.md), R_n retains only the
odd-indexed roots. The even-indexed roots of P_n are the squares of
the roots of P_k. Thus the factorization is

\[
 P_n(Y)=R_n(Y)\prod_{P_k(\lambda)=0}(Y-\lambda^2)
\]

for n≥8. The new diagnostic isolates precisely the primitive half that
governs B, rather than all fourth-energy collisions.

## Exact fiber formula

Reduce R_n modulo p. All its roots lie in F_p. Write c_i≥1 for the
multiplicities of its distinct roots; then \(\sum_i c_i=k/2\). We have

\[
 \boxed{B=2k\sum_i c_i^2.} \tag{3}
\]

**Proof.** A representation x=a+b, with a∈K and b∈L, determines
h=b/a∈L and requires 1+h∈xK. Conversely any such h determines the
representation uniquely. The map h↦h⁻¹ preserves L and has no fixed
point there. For an inverse pair,

\[
 1+h^{-1}=(1+h)/h.
\]

Since h∈L, exactly one of these two elements lies in xK and the other
lies in xL whenever their common H-coset is xH. Both have n-th power
xⁿ. Each multiplicity c_i therefore equals (1_K*1_L)(x) throughout
one H-coset, containing 2k elements. There is no representation of zero.
Summing the squares proves (3). ∎

This also proves that the largest fiber multiplicity equals
\(\max_x|(K\cap(x-L))|\): the polynomial test is an exact encoding
of a shifted-coset intersection problem.

## Squarefree defect and a quadratic class

Put

\[
 a=\deg\gcd(R_n,R_n')=\sum_i(c_i-1).
\]

There is no small-characteristic exception here: p>n>deg R_n. Then

\[
 \boxed{k^2+4ka\le B\le k^2+2ka(a+1).} \tag{4}
\]

Indeed, writing t_i=c_i−1 gives
\(B-k^2=2k\sum_i t_i(t_i+1)\). For integers t_i≥0 with sum a,
the latter sum lies between 2a and a²+a. Thus
\(a=O(\sqrt{k\log k})\) is a sufficient condition for the required
near-quadratic B bound. No such uniform defect estimate is proved.

The upper inequality can be sharp for an actual subgroup: at p=2113,
n=32, the fibers comprise five singletons and one triple. Here k=16,
a=2 and B=448=k²+2ka(a+1). A small defect does not automatically
exclude a triple root.

For a more direct quadratic class, suppose

\[
 \boxed{\gcd(R_n,R_n',R_n'')=1\quad\text{in }\mathbb F_p[Y].} \tag{5}
\]

This is equivalent to having no root of multiplicity at least three.
Every c_i≤2, so (3) immediately gives

\[
 \boxed{B\le2k^2.} \tag{6}
\]

Repeated roots are allowed in (5); this class strictly extends the
collision-free case B=k². If (5) held at every level of a target dyadic
tower, the preceding
[projection induction](parallel-subgroup-2026-09-04.md) would give
\(E_2(H_s)\le2(7+4\sqrt3)s^2\) at every level.

**This implication is conditional.** Condition (5) is another expression
of the unproved pointwise collision bound \(1_K*1_L\le2\). It does
not yet give an arithmetic reason why that bound should hold uniformly.
For each fixed n, all but finitely many splitting primes even make R_n
squarefree, since its integer discriminant is nonzero; this familiar
fixed-order fact gives no growing-order quartic family. In the quartic
window, this pass proves (5) only for the finite cases checked below.

More generally, let

\[
 a_j=\deg\gcd(R_n,R_n',\ldots,R_n^{(j)}).
\]

Then \(a_j=\sum_i\max(c_i-j,0)\), and hence the exact identity

\[
 \boxed{B=k^2+4k\sum_{j\ge1}a_j.} \tag{7}
\]

Thus controlling the full derivative-defect sum at O(k log k) is
equivalent to the near-quadratic B obligation. The shorter first-defect
condition used in (4) is sufficient and can be stronger.

## Literature screen and quantitative gap

Shkredov's [Theorem 6](https://arxiv.org/html/1504.04522) gives
\(E^\times(\Gamma+x,\Pi+y)\ll|\Gamma||\Pi|\log\min(|\Gamma|,|\Pi|)
+|\Gamma|^2+|\Pi|^2\) when \(|\Gamma||\Pi|<p\) and x,y≠0.
For equal subgroups of order n this is O(n² log n), applicable in the
quartic window. This controls multiplicative energy of shifted sets,
whereas B is an additive mixed energy. In the earlier exact transfer,
\(E_2(H)\le3n^2-n+n\mathcal X(H)/3\), it only gives O(n³ log n).
It therefore does not imply (4)'s defect hypothesis or (5). Subtracting
the known trivial shifted-energy solutions does not improve that upper
bound by a power of n. No best-current estimate is asserted here.

## Finite certificates and a new order-512 resonance

The checker reconstructs R_n over the integers for n=4,8,16,32,64,128
using (2), independently reconstructs reductions from their roots,
computes polynomial gcds, and compares the fibers with literal K×L
pair sums. It verifies 46 earlier exceptional fields through order 32,
including genuine triple and quadruple fibers outside the quartic window,
and both earlier quartic witnesses at orders 64 and 128.

It additionally checks the first 2048 primes p≡1 mod n starting at
n⁴/4+1 for each n=256,512,1024. These 6144 fields are a deterministic
finite prefix, not a density theorem or a complete classification.
Primality is checked by trial division by all required primes. Every
checked field satisfies (5). One has repeated primitive fibers:

\[
 p=17189277697,\quad n=512,\quad k=256,\quad
 g=730170861.
\]

The generator has exact order 512. Its fibers are 124 singletons and
two doubletons. Therefore

\[
 a_1=2,\quad a_2=0,\qquad
 B=512(124+2\cdot4)=67584>65536=k^2.
\]

Literal pair counting independently verifies this value. It is a
counterexample to uniform collision-freeness at parent order 512, and
a nontrivial example in the quadratic class (6). The earlier order-512
example at p=17179869697 was circular; this is a different field.
The additional exact counts are E₂(K)=195840, T=0 and E₂(H)=797184.
There are no nontrivial repeated-entry orbits and exactly one
four-distinct orbit, so the parent excess is 24·512=12288.

Run `python3 experiments/parallel2_subgroup_2026_09_04.py`.
The result file records the scan endpoints, complete collision profiles,
polynomial coefficients, source hashes, and exact checks. Nothing in
these calculations bounds the centered balanced moments at logarithmic
depth required by CM or PM+.
