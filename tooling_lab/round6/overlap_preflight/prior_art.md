# Prior-art and local-abstraction audit — 2026-09-05

Classification: **known foundation, locally added support-selection and
certificate capability; historical originality not established**. The exact
intersection covering condition, its greedy optimization, and its use with
interpolation cannot be called new mathematical inventions.

| Primary work inspected | Direct relevance | Boundary |
|---|---|---|
| Gordon, Kuperberg, Patashnik, [New constructions for covering designs](https://arxiv.org/html/math/9502238v1), 1995, Sections II and VI.3 | Greedy covers and complement duality between covering designs and Turán systems. For pure k-supports our condition is T(n,s,k)=C(n,n-k,n-s). | This prototype's greedy policy and finite cover are not new foundations or a new optimal-design result. |
| Li and van Rees, [Lower Bounds on Lotto Designs](https://www.researchgate.net/profile/G-H-Van-Rees/publication/250329694_Lower_Bounds_on_Lotto_Designs/links/0046352a003de957c3000000/Lower-Bounds-on-Lotto-Designs.pdf), 1999 author paper, pp.1–4 | Their lotto-design definition requires every target set to intersect a selected block in at least a specified number of points. This is precisely the overlap condition for fixed support size, with a mixed-size extension here. | “Retaining overlap” alone is not a new combinatorial concept. The allowed blocks here come from verified polynomial tracks, rather than arbitrary blocks. |
| Coffey and Goodman, [The complexity of information set decoding](https://authors.library.caltech.edu/records/c2a7x-nzr02), IEEE TIT 36(5), 1031–1037, 1990; [author-uploaded full text](https://www.researchgate.net/publication/3077637_The_Complexity_of_Information_Set_Decoding), Section II/Theorem 1 | Information sets reconstruct codewords. The paper relates selecting them to covering error patterns and analyzes fixed-rate covering exponents. For RS codes every k-coordinate set is an information set. | Even the covering/decoding bridge and entropy cost are classical. Our certificate concerns an entire affine pencil and preserves every qualifying node, but that specialization does not establish originality. |
| Gordon, Kuperberg, Patashnik, Spencer, [Asymptotically optimal covering designs](https://arxiv.org/abs/math/9511224), 1995 | Counts blocks by the number of target subsets they cover and obtains asymptotically optimal constructions for fixed block size and strength. | Its fixed-parameter asymptotic statement must not be used as a fixed-rate theorem. Our elementary fixed-rate lower bound is derived explicitly in README; the information-set paper is the closer fixed-rate precedent. |
| Ben-Sasson, Carmon, Ishai, Kopparty, Saraf, [Proximity Gaps for Reed–Solomon Codes](https://eccc.weizmann.ac.il/report/2020/083/), revision 3, 2021 | Studies classical algebraic decoding on formal affine-space elements over rational-function fields. This is the relevant larger affine-family setting. | This tiny cover does not reproduce their proximity theorem or address all prize hypotheses. Its event is plain threshold agreement, with all scalar/codeword nodes retained. |

The most decisive local precedent is
`round2/proximity/compressed_support.py:run_census` and its README's
“Mechanism and completeness proof.” It selects a family of interpolation
bases that hits every s-set, reconstructs affine tracks, computes all scalar
buckets, and merges coincident tracks/nodes. Its README already identifies
cross-track overlap as a next optimization. `round4/support_batches/` adds
ZDD representations of basis families and removes all bases belonging to one
maximal joint support. Distinct tracks cannot have k common support points,
so those basis batches are disjoint. Neither inspected implementation
performs the new lane's greedy set cover over the actual maximal-support
catalog. This is an implementation gap within the known abstraction, not a
new completeness theorem.

The round6 partition compiler proves completeness through a disjoint
partition and a sum of root-count caps. The present certificate directly
proves all-s-set coverage and therefore bypasses that sufficient condition.
The independent barrier census proves why no partition can pass in the two
saved fixtures. Combining those artifacts creates a concrete local insight:
failure of the partition representation is distinct from failure of the
broader information-set-cover representation. The known foundations already
predict this distinction; the finite certificates make it testable on the
saved mathematical data.

Search/read log, all accessed 2026-09-05:

1. Local `rg` for `greedy|overlap|Turan|Turán|covering design` in round2
   proximity, round4 support_batches, and round6 cover_barriers; inspected
   both old README proofs and the interpolation/census implementations.
2. Web query `covering design Turan system greedy algorithm covering subsets
   Gordon Kuperberg Patashnik 1995`; opened arXiv abstract/full HTML and
   read Sections II/VI.3 plus the stated counting bounds.
3. Query `generalized covering designs intersection lottery covering m t`,
   followed by `Li van Rees Lower Bounds on Lotto Designs pdf manitoba`;
   opened the author-paper PDF, including its exact definition. Secondary
   mirrors and search snippets were discovery aids, not technical evidence.
4. Query `Reed Solomon covering information sets decoding interpolation`;
   found Freudenberger's doctoral thesis Section 2.3 as a citation lead.
   Followed its Coffey–Goodman citation with the exact-title query; opened
   the Caltech institutional record and the author-uploaded full article.
   The institutional PDF download failed via an expired redirected URL;
   the author-uploaded paper supplied readable full text instead.
5. Opened `math/9511224` to check the fixed-k,t scope of the asymptotic claim.
6. Rechecked ECCC TR20-083's current revision. Its July2021 correction fixes
   a cited lemma's proof without changing its main statements.

Unsuccessful fetches: the guessed Gordon author PDF paths returned errors;
arXiv full HTML replaced them. A Manitoba thesis PDF was inaccessible; the
original Li–van Rees paper was available. No missing fetch was treated as
positive novelty evidence. The search is bounded and does not certify the
absence of prior implementations of this exact support optimizer.
