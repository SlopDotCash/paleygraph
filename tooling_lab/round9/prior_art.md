# Prior art and contribution boundary

Research checked September5,2026. A bounded literature search cannot prove that an approach has never been invented or privately tried. No such claim is made.

The conference identity and character convolution are standard. The iterative NTT implementation develops the lab's frozen [round5 exact NTT](../round5/trace_backend/ntt_exact.hpp); it adds transform lengths and signed CRT range rather than a new transform. The elementary insertion expansion, conditioning, quadratic contractions and Cauchy–Schwarz are known operations. [Döbler–Peccati](https://arxiv.org/abs/1802.00394) study contractions for quantitative limit theorems of symmetric U-statistics; that is a relevant precedent, not a limit theorem applied here.

The earlier [marked-moment compiler](../round4/marked_moments/conditional_moments.py) already computes general conditional second moments. This round's formula agrees with it after inclusion-exclusion. The proposed contribution is the exact fixed-deletion-pair query using selected contractions and a single character convolution, its independent large-field scalar reconstruction, and explicit actual-input witnesses. General conditional variance is not a new invention.

[Gaar–Pucher, v2](https://arxiv.org/abs/2412.12958v2), studies exact subgraph constraints and a vertex-transitive hierarchy for Paley stability bounds. This is nearby work on information retained by local constraints, with a different objective and optimization interface. We have not established that the present interface is absent from all related work.

Queries included “Paley graph local subgraph counts exchangeable pairs second order contractions” and “Reed Solomon proximity line affine pencil agreement supports overlapping information sets.” They returned these broad precedents and the folded-code paper audited in [coding_audit.md](coding_audit.md). Search results asserting broad prize resolutions were not accepted as verified theorem status.

The local [round7](../round7/README.md) and [round8](../round8/README.md) experiments are essential prior inputs. The saved non-affine twins already shared complete one-swap distributions. This round checks their actual second neighbourhoods and discovers that conditioning on the deletion pair reveals a difference even in the mean. Their global two-swap variance still agrees. The new Q term is required for the exact formula but has a rigorously tiny possible effect on the tested large-input variances. These findings argue against making Q the next dominant-obstruction hypothesis.
