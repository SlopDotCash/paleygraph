# Independent Paley–Peisert twin review

Review date: 2026-09-05. Construction, pair-type certificates, the complete
six-subset census, and the statement in `pair_type_twins/universality.md`
are supported. No unresolved correctness error was found. Round-two
sources were read without modification.

## Independent arithmetic and census

Run `python3 tooling_lab/round3/novelty/review_twins.py`. The script uses
its own arithmetic in F7[X]/(X²+1), not the candidate field adapter.
The modulus is irreducible because -1 is not a square in F7. Direct
multiplication and the conjugate/norm inverse formula verify all 48
nonzero inverses. The element 9 encodes 2+X and has order 48; the checks
at powers 16 and 24, and enumeration of its 48 distinct powers, agree.

The independently reconstructed Paley connection set is the set of
nonzero squares, C0 union C2. The Peisert set is C0 union C1. Since
-1=g^24 belongs to C0, both connection sets are closed under negation,
so the resulting zero-diagonal matrices are symmetric. Both matrices
match the saved entries exactly. The explicit trade removes 294 edges
and adds 294 edges, as recorded.

The review reconstructs all **4,802 ordered row-pair histograms** and
inner products across the two matrices. Each graph has three histogram
types, occurring 49, 1,176 and 1,176 times. These are matching multisets;
the labeled map (i,j) to its histogram is not equal. In fact its zero-site
entries reveal the edge sign, so a complete labeled table would recover
the graph.

The candidate census's XOR method is correct: symmetry turns negative
row masks into column masks; their XOR records product parity, while
excluding the six marked rows removes the factors that vanish on the
diagonal. Thus the kernel is 43 minus twice the number of remaining
negative products. Its increasing six-index loops enumerate each
six-subset once, and all integer quantities fit the used types at q=49.

The independent C++ oracle, `twin_origin_census.cpp`, uses ordinary
integer products with their actual zeros, not XOR or popcount. It
enumerates the **1,712,304 six-subsets containing zero**. For any
translation-invariant kernel-value bin, double counting pairs (A,t)
with 0 in A+t gives `6*full_count = 49*origin_count`. This identity also
holds for the joint two-graph value bins, because both graphs are Cayley
graphs for the same additive group. Scaling every bin by 49/6 exactly
reproduces the saved **13,983,816-set marginal and joint histograms**.
No unproved orbit-size assumption is needed.

All 26 saved individual bin witnesses and all saved moments were checked.
The displayed set `[0,1,7,41,47,48]` gives values 19 and -21; the full
joint histogram confirms that its absolute difference 40 is maximal
under this fixed labeling. Evidence and source bindings are recorded in
`twins_review.json`.

## Universality statement and inference scope

The proof in `universality.md` uses precisely its listed hypotheses:
symmetric signs, zero diagonal, zero row sums, and S*Sᵀ=qI-J.
For distinct rows, the four nonzero type counts follow by solving the
four equations for their count, two sums and inner product. Multiplying
both signs by the mutual entry normalizes the table. The coefficient
with exponents (r,d-r,d-r) is unchanged because the Y/Z exponent sum
is 2(d-r). Consequently the aggregated overlap coefficients are
universal even when d is odd. This is separate from the even-degree-only
vanishing assertion repaired in round two.

The marked containment probability depends only on q,n,j,r,d.
Distance correlations are therefore universal. The known Johnson
distance eigenmatrix recovers the harmonic squared norms, since this
membership polynomial has degree at most d and d<=min(n,q-n). One may
use all distance operators; equivalently, the first d+1 rows are an
invertible polynomial-evaluation system on the distinct Johnson degree
eigenvalues. Thus the entire stated L2 energy spectrum depends only on
q,n,d. No finite-field multiplicativity assumption is needed.

The unequal kernel histograms certify nonisomorphism: simultaneous
vertex relabeling permutes the rows and the six-subsets, preserving the
value multiset. A value difference on one fixed labeled subset alone
would not prove this. The histograms also prove that the aggregated
pair-type/L2 signature does not determine the full kernel distribution
within this larger class of matrices.

Both graphs have minimum -29 and maximum 27. This example therefore
does not show unequal extrema, nor rule out a universal upper bound
using the shared data. It concerns an extension field with q=49 and
six-element inputs, not the prime-field critical regime. Neither the
Peisert matrix nor the universality argument creates a counterexample
to a Paley conjecture or a proximity-prize statement.

## Prior art and novelty classification

- **Known graphs and known nonisomorphism:** [Peisert, *All
  Self-Complementary Symmetric Graphs* (2001)](https://www.sciencedirect.com/science/article/pii/S0021869300987143)
  classifies these families. [Alexander, *Designs from Paley graphs and
  Peisert graphs* (2015), introduction](https://arxiv.org/html/1507.01289)
  explicitly gives the two constructions and cites Peisert's Lemma 6.2
  for their nonisomorphism when q differs from 9. None of this is new.
- **Known higher-structure comparisons:** Alexander constructs designs
  from subgraph counts; [Bhowmik–Barman, *Number of complete subgraphs of
  Peisert graphs and finite field hypergeometric functions* (2022)](https://arxiv.org/abs/2205.03928)
  evaluates four-clique counts and related character sums. Comparing
  higher configurations beyond strongly regular parameters is established.
- **Known executable construction:** [Mullin's author-hosted Sage
  routines](https://www.math.uwaterloo.ca/~nmullin/SAGEcode.html) already
  construct Paley and Peisert conference graphs. The present GF49
  implementation is not a new graph construction.
- **Already on the local radar:** the proximity note
  `docs/kb/deltastar-444-arxiv-charp-sweep-new-papers-2026-06-21.md`
  records [Brouwer–Goryainov–Shalaginov–Yip, *Cliques in Paley graphs of
  square order and in Peisert graphs*](https://arxiv.org/abs/2503.09914).
  Its current revision is April 2026; this review verified that revision's
  metadata but does not use its clique theorems in the certificate.

The defensible local contribution is using this standard pair as an
exact adversarial test of the newly implemented energy compiler,
including a reproducible kernel-distribution certificate. This bounded
search did not locate the same six-product histogram calculation. It
does not establish that the pairing, histogram, or universality
calculation has never appeared elsewhere. The universality argument is
an elementary consequence of established intersection-count and Johnson
scheme methods, not a claimed new association-scheme theorem.

Literal web queries on the review date:

1. `Peisert graph GF 49 nonisomorphic Paley graph primitive generator fourth powers union`
2. `Peisert all self complementary symmetric graphs 2001 Paley 49 classification pdf`
3. `Paley Peisert graphs same parameters signed character row pair counts quartic kernel subsets`
4. `"Peisert" "six" "character"`
5. `"Paley" "Peisert" "Johnson"`
6. `"Peisert" "correlation" "distribution"`

The final narrow queries yielded mostly unrelated material; only the
primary sources identified above support the historical statements.
