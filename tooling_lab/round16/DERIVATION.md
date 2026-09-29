# Conditional query covers, actual decoding, and the cost of a missing singleton

The input is an affine polynomial space of dimension three over F_p. Coordinate i has value o_i+theta.v_i. Fix a nonzero projective class: v_i=lambda_i*w with lambda_i!=0. Let

```
b_i=(received_i-o_i)/lambda_i,
m_b=#{i in the class : b_i=b}.
```

A parameter with theta.w=b agrees with exactly m_b coordinates in this class. This histogram is received-word information that a geometry-only all-s-set query cover discards.

## Complete routing with explicit parameter lifts

For a cutoff h, put H={b:m_b>h}. If b is outside H, including values absent from the histogram, m_b<=h. Every qualifying parameter of that kind must have at least s-h agreements outside the class. One residual dimension-three cover at threshold s-h therefore contains it.

For each b in H, impose theta.w=b. Normalize w so its first nonzero entry, at pivot j, is one. The free parameters are theta_l for l!=j, and

```
theta_j=b-sum_(l!=j) w_l*theta_l.
```

The residual affine origin is o+b*v_j and its two directions are v_l-w_l*v_j. Every qualifying parameter in that fibre has at least s-m_b agreements outside the class. Its dimension-two cover supplies it. If a residual threshold exceeds the number of remaining coordinates, that branch contains no qualifying parameter and can be omitted. This prototype explicitly declines thresholds smaller than the relevant residual dimension.

The light branch excludes H, and each heavy branch has a different fixed b. Their accepted parameter sets are disjoint. Solving every query, filtering each complete branch transcript, lifting its surviving parameters and verifying total agreement therefore returns precisely the list inside the original supplied space. If puncturing shortens the domain below the original degree bound, reduction modulo the punctured vanishing polynomial preserves each parameter's evaluations; the final output still uses the original basis.

## Query count is only one cost

The initial compiler chooses h using balanced complete-block query estimates. On the three length1024 tests its query counts are71113,71113 and70892, versus136155 for the geometry-only baseline. The actual decoder confirms that the lists are unchanged. However, the third case requires1467372 explicit coordinate evaluations, versus138310 for the baseline: a query reduction accompanied by a verification regression.

The reason is visible in the transcript bounds. If a false candidate occurs in one rank-d query in a block, its initial agreement upper bound often equals the threshold. In an absent block of size b, the prior cap is d-1. With no matches, that cap first decreases after b-d+2 tested coordinates. Thus a rank-three absent block of size20 may need19 tests, while a rank-two block of size29 may need29. The first conditional family has no singleton in some branches. The independent full-word review confirms the resulting long rejection paths on these inputs.

The revised profile requires at least one unqueried coordinate in each branch. It does not reject a candidate merely for disagreeing with that coordinate. Instead, the same exact agreement-bound update uses the test outcome. The transcript's initial surplus still protects true candidates that disagree with a singleton. A tiny geometry can have no profile satisfying the extra condition, which is reported explicitly.

The resulting query counts are71520,71520 and73554. Coordinate evaluations are74419,72345 and76276, including branch output materialization and final whole-word checks. These are measured arithmetic counts, not wall-clock speedups or universal rejection bounds. The late near miss still needs661 tests within its fibre. A bounded profile choice and a favorable finite workload do not establish optimal adaptive testing or a general decoder for all ambient polynomials.

The substantive design lesson is to optimize the query certificate and its rejection information together. Minimizing the query-family cardinality alone can choose a worse computation. Existing coordinate conditioning and pruning remain prior art; this round provides exact certificates, a demonstrated cost failure, and a verified repair of that failure on the stated inputs.
