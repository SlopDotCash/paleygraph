# Pass 38: energy-source discrepancy resolved

The full Paley and Proximity Prize goals remain open. This pass removes an unsupported route to a stronger imported energy bound; it does not improve the current uniform period exponent.

Two parallel reviews recovered and independently checked the [published source](parallel38-published-energy-2026-09-06.md). Its Theorem 8 uses energy exponent 32/13, whereas the older arXiv v3 states 22/9. Since 32/13 > 49/20, the published result is consistent with the later MRSS improvement already used by this project. The existing 49/20 input and repeated-six exponent 129/40 remain unchanged.

The [independent incidence audit](parallel38-independent-incidence-2026-09-06.md) identifies an added published hypothesis: the map from three subgroup cosets to their two ratios must be injective. It fails for a triple of identical sets containing multiple subgroup cosets. The older proof applies its incidence bound to such diagonal triples without this condition.

The [exact finite verifier](../experiments/parallel38_coset_multiplicity.py) checks this on the older proof's own dyadic construction at p=1153, |H|=8. One level has 24 elements; its ratio image has size 1600, while the published hypothesis requires 1728. A second example with nested subgroups in F_97 shows why retaining multiplicity matters: the correct incidence count is 48, whereas forgetting the multiplicity gives 16. These are failures of a hypothesis and of an unweighted counting substitution, not counterexamples to an asymptotic 22/9 bound.

The audit also proves the exact weighted counting identity. A level-set argument gives a conditional replacement with a factor (maximum coset multiplicity)^(1/3), provided the distinct-pair incidence estimate applies at every level. Controlling that factor is extra work; it cannot be discarded to recover the preprint's argument.

The publisher's raw PDF, HTML provenance, and two decisive page images are hashed in the [source manifest](../sources/parallel38-energy-source/manifest.json). Root and both reviewers visually checked the relevant published pages. The finite enumeration was independently reproduced. No Lean build, Prove2Me submission, public claim, or change to the main source manifest was made.

Next mathematical work must use the retained, scoped energy inputs or prove an additional structural estimate. The required positive moment bound, uniform lattice-distance deviation bound, classical two-set Paley theorem, and official prize conclusion are all unproved.
