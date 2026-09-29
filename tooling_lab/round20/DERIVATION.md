# Exact completions, representation cost, and missing equivalence information

This round uses the centered cyclotomic digit section from rounds17–19. It studies how to represent and query the entire finite section, and how to certify lower bounds when that representation is unavailable. No assertion here proves either prize. Automata minimization, distinguishing suffixes, and observation-table refinement are established methods; the contribution is their exact arithmetic interface and the checked experiments on this project's inputs.

## 1. Fixed-order residual languages

Let L be the finite set of actual digit words in a specified coordinate order. For a prefix P of length d, define its residual language L/P = {S : PS belongs to L}. A deterministic reader at depth d can merge two prefixes exactly when these suffix sets coincide. Empty residuals can all use an implicit rejecting state, excluded from our width counts.

At depth N the only nonempty residual is the empty suffix. Inductively a nonempty residual at depth d is identified by its sorted set of pairs (next digit, child residual). Different digits partition the language, so identical signatures imply equal suffix sets; different signatures imply different suffix sets by induction. Interning these signatures gives the minimum number of nonempty states at every fixed depth for this fixed reading order. This is not a minimum over coordinate orders, arbitrary arithmetic algorithms, skipped-level representations, or unrestricted machine models.

The verifier checks distinct signatures, exact outgoing symbols and child indices, reachability of every stored node, and the suffix-count recursion. Small diagrams are also expanded to every accepted word and compared with independent dense-matrix scalar encodings. Large symbolic diagrams are checked against a separate binary-endpoint recurrence. A graph that is minimal for the wrong language would be useless; both language equality and minimality are checked.

## 2. A complete generator-two constructor

For Q = 2^N+1 and positive k dividing Q, put p = Q/k. In the canonical relation 2-X^-1, the digit section is

L(N,k) = {0} union {D: nonzero digits alternate in sign, their number is odd, and S(D) is divisible by k},

where S(D) = sum_j 2^(N-1-j) D_j. This is the exact round18 language. For composite Q and k=1 it remains a centered residue model; such cases are not asserted to be prime fields.

The producer keeps (first nonzero sign, latest nonzero sign, S-prefix mod k). On a zero digit the signs stay fixed. A nonzero digit must oppose the latest nonzero sign; it initializes the first sign if necessary. Every digit updates the weighted residue. At the endpoint, the residue is zero and the first and latest signs agree. The all-zero path is separate. There are at most 4k+1 reachable states at each depth before minimization, and O(kN) states/transitions overall. Counts use arbitrary-precision integers. This is a cofactor-dependent bound, not a uniform efficient algorithm for every short relation.

The separate verifier uses binary variables b_0,...,b_N with b_N=1-b_0 and D_j=b_(j+1)-b_j. Each nonzero accepted D has exactly one such binary path; the all-zero D has none and is added separately. It recursively counts binary completions for each endpoint/residue state, then checks every reachable diagram transition against the union of these binary possibilities. It does not import the producer's sign-state construction.

At N32,k641, the constructor covers all p=6,700,417 scalar words, including zero, without enumerating p scalars or the ternary ambient cube. It yields 36,381 nodes and 67,668 edges. The exact minimal width is 2,565=4k+1 at depths11 through21. This finite equality does not assert that every cofactor attains the bound.

## 3. Exact width five when k=1

For N>=4, the canonical language has the exact profile

`1, 3, 5, ..., 5, 4, 1`.

At depths2 through N-2, all five sign states are reachable: the zero prefix and the four pairs (first,last) in {+1,-1}^2. They are pairwise distinguishable using at most two further nonzero digits, padding any unused suffix positions with zero. The same-sign states accept the all-zero suffix; the opposite-sign states do not. The (+,-) state accepts a single + and the (-,+) state a single -. The two same-sign states are separated by the two-digit suffix (-,+) or (+,-). The zero-prefix state accepts both one-digit signs and the all-zero suffix, separating it from the others.

At depth1 only the zero prefix and the two same-sign states are reachable. At depth N-1 the two same-sign states merge: both have residual {0}. The other residuals are {-1}, {+1}, and {-1,0,+1}. At depth N only the accepting empty suffix remains. Therefore the minimal layered diagram has exactly 5N-6 nodes and 11N-17 edges. The programs check this through N128, where the language has 2^128+1 words, 634 nodes, and 1,391 edges. That modulus is not prime.

This precisely measures the retained memory that the long raw-window descriptions in round19 obscure. It is a property of an already known signed-digit mechanism, not new finite-automata theory.

## 4. Ordering as an exact control

Round18's orientation adapter sends the saved small encodings to the canonical encoding by a signed coordinate permutation. For saved g=2^e and the corresponding monomial shift h, coordinate j lands at (e*j-h) mod N. Reading the saved coordinates in the order sorted by this position restores the canonical width profile. Per-coordinate sign changes merely relabel outgoing letters and do not change residual equality.

Complete scalar censuses establish:

| Prime | Natural saved peak width | Canonical peak width | Saved encoding, reordered |
|---|---:|---:|---:|
|257|60|5|5|
|65,537|1,597|5|5|

At p65,537 the saved diagram has 6,650 nodes; reordering reduces it to74. Every scalar remains in the represented language, and original-coordinate query assignments are translated explicitly. This is a finite representation improvement, not a measured wall-clock speedup or a new ordering theorem. The general p2,013,265,921 relation lacks this known generator-two adapter.

## 5. Completions with arbitrary fixed coordinates

For fixed coordinates A and values v, define C(d,s) as the number of accepted suffixes from node s compatible with the remaining assignments. At a fixed coordinate retain only its selected edge; at an unfixed coordinate sum all child counts. Set C(N,accept)=1. Memoization visits at most every node once. Enumerating along edges with positive C returns complete words in lexicographic read order. Output size remains an unavoidable cost; the API explicitly reports a full count and whether its listed words are truncated.

The command-line API accepts coordinates in the original digit order, translates them into the stored read order, then translates output words back and recovers their scalar residues. As a concrete large query, D0=0,D7=1,D19=-1 has exactly209,394 completions at p6,700,417. The example lists only five and marks the list truncated. An independent binary recurrence checks its count, and direct centered arithmetic checks every returned scalar/digit pair.

## 6. Certified lower bounds without a complete diagram

Choose actual words W_i=P_i S_i at the same cut. The exact codec defines M_ij=1 if P_i S_j is actual and0 otherwise. All diagonal entries must be1.

The initial certificate selects a greedy clique in the graph whose edge i-j means M_ij=0 or M_ji=0. Either rejection gives a suffix accepted for one prefix and rejected for the other. Every selected pair therefore needs different reader states. A greedy clique is only a lower bound, not a maximum clique.

A stronger use of the same matrix groups prefixes by their whole row. Whenever two rows differ, some supplied suffix column explicitly distinguishes them, even if neither prefix's own suffix does. We store one representative per distinct row and an explicit differing column for every representative pair. Equal observed rows are inconclusive: an untested suffix may distinguish them.

On the saved p65,537 input at depth8, the own-suffix clique certifies144 states, all supplied suffixes certify195, and the complete language has801 states at that cut. The256 supplied scalar prefixes occupy205 of those801 true residuals. The missing ten distinctions among these supplied prefixes are recovered with nine adaptively chosen suffixes from the complete diagram. Independent arithmetic checks all2,304 new prefix/suffix words. This stops with exact equivalence classes for the supplied prefixes, while leaving596 residual classes outside the prefix sample.

At p6,700,417,256 supplied words certify245 states at depth16, against exact full width2,565. At p2,013,265,921 the natural-order certificates give203 states at cut4 and255 at cut8. At cut32 all256 selected prefixes are pairwise distinguishable under each of natural, reversed, and even-then-odd order. These are finite, order-specific lower bounds. No all-orders statement or growth rate follows, and256 is the sample ceiling.

## 7. The missing tool after membership

The general projection codec can classify any fully supplied word. It cannot certify that two prefixes have identical completion sets without considering their possible suffixes. The complete diagrams provide that equivalence oracle in the structured controls, letting the adaptive procedure produce a genuinely new test suffix whenever a false merge remains. For the general N64 input, that complete oracle is unavailable; its adaptive record explicitly says so instead of treating absence of sampled counterexamples as equivalence.

This identifies the next concrete obligation: an arithmetic completion or residual-equivalence certificate that works directly from the scalar section, without a full p-element codebook or a huge cofactor-state table. A useful next prototype must certify both empty and nonempty residuals and expose its cost on adversarial prefixes. Arbitrary norm aggregates and the full problem's quantifiers remain beyond the present tool. No proximity algorithm is changed in this round; its supplied-affine-space scope from round16 remains in force.
