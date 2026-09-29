# Exact covering certificates: prior art and scope

This round applies established branch-and-bound and exact-certificate ideas to the arithmetic difference body derived in round22. It does not establish that these general techniques are new or that the application is a worldwide first.

Cheung, Gleixner and Steffy's [Verifying Integer Programming Results](https://arxiv.org/abs/1611.08832), submitted2016 and published at IPCO2017, proposes sequentially verifiable certificates for mixed-integer linear programming and an independent verifier. The abstract was inspected. This directly precedes the broad idea of separating fast discovery from exact verification.

The primary [VIPR repository](https://github.com/scipopt/vipr) describes exact rational checking, branch-and-cut certificates, and versions of its certificate format. Its [version1.1 specification](https://github.com/scipopt/vipr/blob/master/cert_spec_v1_1.md) supports assumptions and incomplete derivations. This prototype uses its own small JSON format and checker; it neither invokes VIPR nor claims compatibility. A production mathematical tool should consider export to a maintained certificate standard rather than indefinitely expanding a private verifier.

The norm-residual dual certificate and convex-body lattice-enumeration context are audited in [round22's prior-art note](../round22/prior_art.md), including Boyd and Vandenberghe's norm-approximation duality and Dadush, Peikert and Vempala's lattice enumeration in arbitrary convex bodies.

The project contribution being tested is narrower: compile exact scalar-visible difference constraints for this cyclotomic digit code, obtain rational cuts through the arithmetic interface, and emit an exhaustive cover of integer difference boxes without enumerating all their points. The checker verifies every rational identity, interval contraction, disjoint split, root region and unresolved leaf. A finite unique-completion theorem follows only if the entire nonzero difference domain is covered by exclusions.

The later orbit-intersection counter is closely related to the hidden-number formulation in Boneh and Venkatesan's [Hardness of computing the most significant bits of secret keys in Diffie–Hellman and related schemes](https://crypto.stanford.edu/~dabo/pubs/abstracts/dhmsb.html), CRYPTO1996. Their primary abstract was inspected; it studies recovery of a scalar from partial information about its products with known powers. Our constraints are exact centered intervals for such products. This experiment uses complete enumeration from the tightest interval plus exact filtering, not their oracle theorem or a new general hidden-number algorithm.

Queries were run on September6,2026. Search coverage is not a proof of historical absence. The model translation must still be justified: a correct integer-programming certificate proves only the encoded claim. The separate cyclotomic derivation and codebook controls supply that additional interface here.
