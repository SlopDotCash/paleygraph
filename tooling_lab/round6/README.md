# Round 6: information value and unanchored codeword families

Two new extensions are implemented and independently checked. Further bounded work is testing the structural reasons they succeed or fail. This iteration remains in progress; the previous [round5 snapshot](../round5/README.md) is frozen.

| Tool | Verified result | What it does not establish |
|---|---|---|
| [Exact trace envelopes](trace_envelopes/README.md) | Rational certificates quantify the range of global third moments allowed by cheap constraints or ordinary theta moments. Three fields retain exact relaxed ambiguity; q101 has a unique relaxed histogram. | The rational histograms need not be realized by graphs. The global average does not bound individual sets. |
| [Pencils without codeword anchors](pencil_tracks/README.md) | Blind discovery and symbolic certificates give complete whole-field lists on structured n1024 inputs with no exact anchor. Independent tests cover1586 tiny pencils and both large outputs. | Discovery and the strict root-count cap can fail. Ordinary agreement lists are separate from the official MCA event. |

The first tool tests the value of an additional statistic. At q1,000,033,n31, the basic conference identities and Hasse-supported residue domain already constrain the raw third moment to relative interval width below6.50×10⁻¹⁰. Ordinary moments through6 leave a nonzero rational ambiguity, but its certified width is below2.90×10⁻¹⁹ of the actual raw moment. An exact inequality can therefore identify lost information and simultaneously show that its effect on this average is small.

The second tool tests an earlier structural restriction. A verified track is a pair of degree<k polynomials defining `h_j(z)=a_j+z*b_j`. On its coordinate block the input pencil agrees identically with that track. If the root-count cap is strict, every qualifying polynomial at every scalar must equal a track. Agreement and collision equations are linear in the scalar, so a finite exception ledger plus the generic complement certifies the entire field without an exact codeword anchor.

For n1024,k64,s410 over F65537, the two-track fixture has two qualifying codewords at every scalar, totaling131074 nodes. The three-track fixture qualifies only at three track-collision scalars, totaling three nodes. These are selected structured inputs, not a general proximity theorem.

Independent evidence is in [trace_envelopes_review](trace_envelopes_review/README.md) and [pencil_tracks_review](pencil_tracks_review/README.md). The envelope review checks5208 rational dual inequalities,104 primal equalities and56 corruption controls without a production-model import or optimizer. The pencil review independently rebuilds both large event ledgers, verifies no-anchor certificates, and checks the two discovery repairs.

The [current plan](../NEXT_ITERATION.md) continues with exact critical-size moment algebra, an obstruction to every partition-based root-count certificate on the hard toys, and a covering-design prototype retaining overlapping supports. These extensions will determine whether the remaining gap is discovery, the certificate's representation, or scale.

Moment optimization, polynomial reconstruction, root counting, affine tracks and finite-field linear equations are established. The local contribution is their precise validated interfaces, measured information limits, and failure-driven refinements. Historical uniqueness is unestablished; neither prize is proved.
