# Independent review of unanchored polynomial-pencil tracks

2026-09-05. The complete-list compiler and its mathematical contract passed independent review. No source correction was required. The reviewed core hash is `1d2b3ca8dfaaa164155faa55d3b3a83d4e5c99dbfa012daf4911d99d77cc60f1`; the saved experiment hash is `6fadaba99d3f96ccef8d94266d48218c5d2ba9b286c6ff14d5f99a3bcf9c7216`.

[review_tracks.py](review_tracks.py) independently reconstructs coordinate equations, polynomial-collision equations, the exceptional set, the generic complement, and every finite deduplication group. It uses the previously frozen independent prime/GF9 arithmetic oracle. [Its evidence](review_tracks.json) contains1586 complete tiny codebook/scalar comparisons,3345 nodes,759 pencils without an exact code anchor,513 whole-field qualifying sets, and five rejected corruptions. The exhaustive cases include every F3 length3 pencil at k1,s2 and k2,s3, plus64 constructed F5 and64 GF9 inputs.

For both n1024 examples the review reconstructs the entire symbolic ledger with independent arithmetic; it does not rerun discovery or enumerate the full field or codebook:

| F65537,n1024,k64,s410 fixture | Tracks | Cap | Exceptions | Collisions | Generic list | Qualifying scalars | Nodes |
|---|---:|---:|---:|---:|---:|---:|---:|
| Two skew tracks |2|126|1011|0|2|65537|131074|
| Three pair-colliding tracks |3|189|3|3|0|3|3|

The medium F257 examples, saved GF9 certificate, and both successful discovery-iteration certificates are also replayed. The GF9 and F17 iteration certificates receive complete independent codebook/scalar comparisons. Large conclusions follow from the complete symbolic ledger and root-count proof, not the selected scalar spot checks in the original experiment driver.

## Completeness, collisions and no-anchor evidence

A verified disjoint track partition gives a degree<k polynomial `a_j+z*b_j` on each assigned region at every scalar z. Any polynomial differing from every such polynomial has at most `sum_j min(|R_j|,k−1)` matching coordinates. The strict cap is therefore a sufficient complete-list certificate. Overcounting regions when specialized track polynomials coincide only weakens this bound; it cannot create a false exclusion.

Every track-coordinate agreement equation is affine in z, with solution set all, empty, or one field element. The same classification applies to equality of two coefficient vectors. Globally identical tracks are merged. Because n≥k and domain elements are distinct, distinct degree<k coefficients cannot represent the same evaluated codeword. Outside the union E of coordinate roots and polynomial collisions, all supports are fixed and all track polynomials are distinct. At z∈E, deduplication by the complete coefficient vector is exact; equal polynomials necessarily have equal maximal agreement supports. Thus

```
node_count=(q−|E|)*generic_list_size + sum_(z in E) list_size(z)
```

counts every scalar/codeword pair exactly once. Including collision scalars even when no track qualifies is harmless. The independent review checks the all/empty/singleton partition at every coordinate, every collision, every finite grouping, the entire generic scalar range, maximal supports, and both total-node and qualifying-scalar counts.

The large no-anchor conclusions also have independent finite certificates. Every assigned region has at least k points. An exact codeword on the pencil would consequently equal every specialized track polynomial, since each difference has at least k zeros. For the two-track fixtures their coefficient-collision equation has no solution. For the three-track fixtures the three pairwise collision scalars are different, so their intersection is empty. At F65537 these scalars are32767,65535,65536. No scalar can make all tracks equal. This establishes absence of an exact codeword anywhere on the pencil without relying on the production `find_code_anchor` function.

## Discovery repairs and a verified failure

[review_discovery.py](review_discovery.py) and [its evidence](review_discovery.json) check the two new discovery iterations:

- A fresh F17 replay without reassignment reproduces cap3 at threshold3 and correctly refuses a list certificate. Reassignment uses only coordinates satisfying both polynomial equations, gives cap2, and matches all34 scalar/codeword nodes in a complete oracle.
- On the F257 length192 example, the default marginals give four verified joint tracks and cap28<s80, yielding two nodes at one scalar. The alternative scalar-one word has four distinct polynomial pieces with maximal supports48,49,48,48. Each exceeds3(k−1)=21. Any putative three-polynomial cover would have to contain each of those four polynomials by root counting, an impossibility. Since a four-piece cover is explicit, the minimum is exactly4. This independently explains failure of a max-three marginal discovery budget. All three exported dual inconsistency witnesses for degrees1,2,3 were also checked directly against independently constructed row equations, without rerunning the medium discovery solve.
- The nine-joint-region fixture has a verified partition but cap63 at threshold50. A fresh compiler call correctly rejects its complete-list claim. The failure supplies no evidence that the true list is empty.

The tiny compiler review took11.64s, including all saved symbolic readbacks; the separate discovery check took1.01s. No large discovery or dense matrix run was repeated.

## Scope and prior work

The source documentation accurately separates bounded discovery, verified joint covers, strict-cap completeness, ordinary agreement events, and the additional official MCA nonjointness condition. Greedy reassignment has no minimum-joint-cover guarantee. The method can fail because a marginal cover is not found or because intersections create too many regions. GF9 verifies extension arithmetic in a toy setting, with no claim of an official-scale extension-field backend.

The local [July11 fiber-Chebyshev note](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-rate-quarter-fiber-chebyshev-2026-07-11.md) explicitly records the k−1 bound for actual codeword-pair differences and an at-most-one exceptional scalar lemma. It also records failure of a root-count argument for arbitrary input-relative ratios. The present compiler respects that boundary by verifying the polynomial identities on each region before applying the bound. The [July9 coefficient-arrangement note](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-design-matrix-affine-clusters-2026-07-09.md) already develops interpolation pencils, affine clusters and scalar/codeword incidence. These are substantive local predecessors.

The primary [ECCC report for *Proximity Gaps for Reed–Solomon Codes*](https://eccc.weizmann.ac.il/report/2020/083/) describes classical algebraic decoders applied to a formal affine-space element over a rational function field. Its current record lists revision3, July2021. This source and the two local notes were checked2026-09-05; no general proximity theorem from the paper is being reproved or strengthened here.

The specific local addition is the implemented complete scalar ledger for a verified piecewise affine polynomial partition, reached by a bounded discovery pipeline that no longer requires an exact code anchor. Historical novelty, a complete decoder for arbitrary inputs, and a prize proof remain unestablished. All earlier rounds and reviews are unchanged.

Replay `python3 tooling_lab/round6/pencil_tracks_review/review_tracks.py` and `python3 tooling_lab/round6/pencil_tracks_review/review_discovery.py`.
