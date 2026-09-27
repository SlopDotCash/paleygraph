# Primary-source audit of the shifted-product excess

The four sources checked here do not provide the uniform estimate
X(H)<<n^(61/31) needed by the pass41 triangle criterion in the dyadic
quartic range. They do supply the norm characterization behind the
exceptional-prime argument. This is a review of specified statements,
not a claim that all possible literature has been exhausted.

The source lane retrieved these primary PDFs before reaching its account
usage limit. Root read the relevant primary statements and checked their
application below. The archives and hashes are recorded in
[the source scope](../results/parallel42_source_scope_2026_09_06.json).

## Circularity and primitive cyclotomic norms

[Ke and Kiechle, Overlaps in circular planar nearrings, Algebra and
Discrete Mathematics 35 (2023), 134–167](https://admjournal.luguniv.edu.ua/index.php/adm/article/download/2130/pdf),
DOI 10.12958/adm2130, gives the nontrivial shifted-product circularity
criterion in Theorem 2.1 and fixed-order finiteness in Theorem 2.2.
Proposition 3.4, pages 142–143, relates prime divisors of a primitive
cyclotomic norm to existence of a vanishing primitive root in the finite
field. The quantifier is existential: divisibility of one polynomial's
norm does not force that polynomial to vanish at every chosen primitive
root. Corollary 3.7 identifies the circularity exceptions through these
primitive norms. It supports the fixed-order finite-support argument.

Our [compatible norm budget](parallel42-excess-equations-2026-09-06.md)
uses the product over all ordered nontrivial tuples, which is Galois
invariant. At any fixed split root, each actual collision contributes
at least one valuation. Equivalently one can sum determinant nullities
over all primitive roots and divide the cyclotomic degree exactly once.
This retains multiplicities without an invalid fixed-root converse.
The sharper archimedean constant is derived locally by an exact second
moment and arithmetic-geometric mean. No novelty claim is made.

A narrow source caution matters for implementation: the sentence after
Theorem 2.2 on page 136 prints a resultant against x^k-1. The displayed
shifted-product polynomials vanish at 1, so that literal resultant is
zero. Our argument uses the primitive cyclotomic norm in Proposition 3.4
and Corollary 3.7 instead. This observation is not a rejection of those
later statements. Remark 3.6's crude bound for a broader overlap exception
set is exponential in phi(k); it does not control the quartic regime.

## Bounds on individual cyclotomic numbers

[Do Duc, Leung and Schmidt, Upper bounds for cyclotomic numbers,
Algebraic Combinatorics 3 (2020)](https://alco.centre-mersenne.org/item/10.5802/alco.86.pdf),
Theorem 1.2, uses q=ek+1 and the hypothesis

    p>(sqrt(14))^(k/ord_k(p))

to bound every cyclotomic number by 3. Here k is subgroup size, not
subgroup index. In our prime-field problem, k=n and p=1 mod n, so
ord_n(p)=1. Its sufficient lower bound on p is exponential in n and
does not cover growing n with p comparable to n^4.

Theorem 1.3 bounds those numbers by 2 under an additional prime-k
hypothesis and a different exponential condition. Prime k excludes
dyadic n>=4. The introductory asymptotic for fixed index e and growing
q also has different quantifiers from our thin-subgroup problem.

## Circularity as an assumption

[Ke and Kiechle, Circularity in finite fields and solutions of the
equations x^m+y^m-z^m=1, arXiv:2307.05586v2](https://arxiv.org/pdf/2307.05586v2),
Remark 3, recalls fixed-order finiteness and discusses observed table
proportions. That discussion is not a uniform asymptotic density theorem
as conductor grows. The associated published article is Finite Fields
and Their Applications 98 (2024), 102467, DOI 10.1016/j.ffa.2024.102467;
the publication metadata was checked on the
[authors' institutional record](https://researchoutput.ncku.edu.tw/en/publications/circularity-in-finite-fields-and-solutions-of-the-equations-xsupm/).
The archived mathematical text used here is explicitly the arXiv version.

Results that assume circularity concern the case of no nontrivial
shifted-product collisions. Assuming that condition would remove the
very excess that still needs a bound. We do not use circularity as an
unproved premise, or infer an eligible-prime denominator from the tables.

## Shifted multiplicative energy

[Macourt, Shkredov and Shparlinski, Multiplicative energy of shifted
subgroups and bounds on exponential sums with trinomials in finite
fields, arXiv:1701.06192v2](https://arxiv.org/pdf/1701.06192v2),
Theorem 1.2 and Corollary 4.1, give a piecewise shifted-energy estimate.
In the small-subgroup range n<p^(1/2)log p, Corollary 4.1 gives

    E^times(H+lambda) <= n^4/p+O(n^2 log n), lambda!=0.

Deleting zero cannot increase multiplicative energy. Thus the statement
bounds our B=E^times((H-1) minus {0}), but subtracting its exact trivial
part does not produce X<<n^(2-delta) for any delta>0 in the quartic
range. Remark 1.3 explicitly places the new collinearity improvement in
the intermediate-size range. Corollary 4.2 supplies set-size consequences,
not the missing nontrivial collision saving.

The pass42 result remains an absolute bound on the number of exceptional
primes. It supplies neither a bound for every prime nor a reduction from
this intermediate triangle estimate to the full Paley or prize targets.
