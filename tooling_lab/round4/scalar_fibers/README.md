# Exact scalar fibers after a complete codeword intersection

This iteration resolves the preceding interleaved-fixture failure by transporting every non-anchor scalar to one static direction-word decoding problem. A verified piece-polynomial partition then certifies the complete static list. The output retains the polynomial, full codeword and maximal support of every symbolic fiber, with a separate anchor node.

The mathematical foundations are elementary linear-code symmetry and polynomial root counting. Related through-code, singleton-anchor and scaling reductions already exist in the local research repository. This is an extension of that existing diagnostic into a complete executable certificate, not a newly discovered theorem or a claim of historical novelty.

Here, a scalar is counted when its word has **some codeword agreeing on at least s coordinates**. This is ordinary closeness. It does not impose the MCA event's additional exclusion of jointly explained agreements. In particular, the whole-field examples below do not claim a whole field of MCA violations.

## Transport certificate

Suppose the input pencil `w(z)=u0+z*u1` meets the code exactly at `z=z*`, with `w(z*)=f` a degree-less-than-k polynomial evaluation vector. Then

`w(z)=f+(z−z*)u1`.

For every `z≠z*`, code linearity gives a bijection

`h ↦ g(z)=f+(z−z*)h`, with inverse `h=(g−f)/(z−z*)`.

Coordinate equality is preserved exactly:

`g(z)[i]=w(z)[i]  ⇔  h[i]=u1[i]`.

Thus the complete lists at every non-anchor scalar have equal cardinality, the same support multiset and explicit codeword correspondences. At the anchor, all these tracks meet at f, whose maximal support is the whole domain. The tool retains this one anchor node separately. Because `s≥k`, any other degree-less-than-k codeword agreeing with f on at least s distinct coordinates must equal f.

The implementation rejects `s<k`; it does not silently drop extra anchor codewords. An exhaustive F5 test with n5, k2 and s1 finds21 codewords at the zero anchor, demonstrating the need for this guard.

The anchor can be found without scanning the field. Interpolate u0 and u1 on the first k points, subtract their evaluations, and solve the two residual vectors' common affine zero. A consistent solution yields the complete codeword intersection. A zero direction residual with nonzero offset residual proves there is no intersection; two zero residuals mean the entire pencil lies in the code. All domain and field encodings are validated before interpolation.

## Static-list certificate from supplied pieces

The user or another tool supplies a partition of the coordinates and a polynomial for each part. The verifier checks that these polynomials have degree less than k and agree with the direction word u1 throughout their assigned parts. Repeated polynomials are merged, and their regions are joined.

For distinct piece polynomials `p_j` on disjoint regions `R_j`, any degree-less-than-k polynomial h outside this list satisfies

`|{i:h(i)=u1(i)}| ≤ Σ_j min(|R_j|, k−1)`.

Each term follows from the root bound for the nonzero polynomial `h−p_j`. If s is **strictly larger** than the sum, every qualifying polynomial is among the supplied pieces. Evaluating each candidate on the full domain then yields the exact static list and its maximal supports, including coincidences outside the assigned region.

This verifies a supplied decomposition; it does not discover one. The scale fixtures provide their planted piece polynomials and partition explicitly. Finding such witnesses for an arbitrary hard input remains an unresolved algorithmic task.

The strict inequality matters. Over F7, the word `[0,0,3,3]` on points `[0,1,2,3]` is covered by two constant pieces. With k2 and s2, the outside polynomial `h(X)=X` also qualifies on coordinates `{0,3}`. Since the root-count cap equals2, the tool rejects this insufficient piece certificate.

## Results and comparison with the failed batching approach

All inherited aligned and interleaved fixtures retain exactly the earlier complete node outputs. Every row below uses three verified piece polynomials; anchor discovery additionally uses two k-point interpolants.

| Inherited fixture | Previous distinct tracks | Verified direction pieces | Complete output |
|---|---:|---:|---|
| Aligned n48, k6, s28 | 3 | 3 | Only anchor z=1 |
| Interleaved n48, k6, s28 | 546 | 3 | Only anchor z=1 |
| Aligned n128, k8, s80 | 3 | 3 | Only anchor z=1 |

Two further cases use the order1024 subgroup of F65537, k64 and s410, with the pieces interleaved coordinate by coordinate:

| Case | Outside-piece cap | Exact piece agreement sizes | Non-anchor list size | Complete symbolic node count |
|---|---:|---|---:|---:|
| Three pieces | 189 | 342,341,341 | 0 | 1 |
| Two pieces | 126 | 512,512 | 2 | 131,073 |

The three-piece case certifies only the anchor scalar. The two-piece case certifies ordinary closeness at every scalar: two distinct codewords at each of65,536 non-anchor scalars and one merged full-support codeword at the anchor. Its certificate represents those131,073 nodes with two polynomial fibers plus the anchor, preserving all supports. The saved construction times were0.169 and0.141 seconds; input generation and separate readback are excluded. These are finite fixture timings, not asymptotic guarantees.

The earlier n1024 selected cover contains approximately `3.57×10^47` bases. That enumeration was not attempted on the new interleaved inputs. This certificate bypasses it using the supplied decomposition and the complete static-list proof. The example does not establish that hard prize inputs have a comparably small piece witness.

The main suite checks **295 complete tiny pencil oracles** over F3, F5 and GF9, totaling1,185 exact nodes;52 cases also compare the static list directly against every coefficient vector. It includes241 whole-field ordinary-closeness cases, entire-code pencils and a two-fiber GF9 example. Six deliberate certificate corruptions are rejected. Controls also reject absent anchors and the two invalid parameter regimes above.

A separate GF(3^6) test uses n4, k1, s2, two interleaved constant pieces and non-base-field elements. It independently enumerates all **531,441 scalar/constant-codeword pairs**, comparing all1,457 nodes to the symbolic certificate. The field degree matches the profile's degree-six requirement, while its field size and code parameters remain tiny.

## Reproduction and artifacts

From the repository root:

```sh
python3 tooling_lab/round4/scalar_fibers/run_scalar_experiments.py
python3 tooling_lab/round4/scalar_fibers/degree_six_control.py
```

- [scalar_fibers.py](scalar_fibers.py): anchor discovery, verified static piece list, symbolic scalar transport, per-scalar materialization and certificate readback.
- [results.json](results.json): complete large symbolic certificates, supplied piece witnesses, inherited equality checks, tiny-oracle totals and counterexamples.
- [degree_six_results.json](degree_six_results.json): complete GF729 control and its symbolic certificate.

`materialize_scalar` evaluates one certified scalar fiber on demand. It returns all qualifying codewords with their maximal agreement supports; it does not search for candidates. The field coefficients use the existing polynomial-basis integer encoding.

The production readback reconstructs the proof data and shares field/polynomial helpers with the generator. The tiny oracle separately enumerates every polynomial and scalar using direct monomial evaluation. The GF729 constant-code oracle enumerates each field constant directly.

An [independent review](../novelty/scalar_fibers_review.json) passed1,700 complete coefficient/scalar comparisons totaling4,400 nodes,729 exhaustive anchor searches,108 static-list checks,54 strict-cap rejection controls and40 independent GF9 cases. It included904 whole-field cases,486 of them entire-code pencils, and read back six saved certificates including n1024. Its [script](../novelty/review_scalar_fibers.py) also rejected five corrupted exports. This review covers the frozen core and main results; the separate GF729 control is additional evidence from this lane.

## Prior art, complexity and remaining boundary

Local sources already include:

- [Global consistency charge, §3](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-rate-quarter-global-consistency-2026-07-10.md): the zero-defect through-code case has no non-anchor riders on sub-threshold tracks; its source contains `base_row_stack_carries_no_riders`.
- [The second-witness floor, singleton-anchor theorem](/Users/shawwalters/proximityprize/ArkLib/Data/CodingTheory/ProximityGap/Frontier/_SecondWitnessFloor.lean:326): the exact singleton witness fiber at a codeword intersection when s≥k.
- [MCAEventCosetInvariance.lean](/Users/shawwalters/proximityprize/ArkLib/Data/CodingTheory/ProximityGap/MCAEventCosetInvariance.lean) and [FarLineScalarDilation.lean](/Users/shawwalters/proximityprize/ArkLib/Data/CodingTheory/ProximityGap/FarLineScalarDilation.lean): existing code-coset and scalar-dilation invariances. Their event definitions and hypotheses should not be conflated with ordinary closeness here.

The root-count mechanism is standard polynomial reconstruction; Guruswami and Rudra describe the familiar degree-versus-number-of-roots principle in [Error Correction Up to the Information-Theoretic Limit](https://www.cs.cmu.edu/~venkatg/cacm09.html). A bounded source search did not find the same exported complete static-piece/transport certificate, but the overlapping local mathematics is already explicit. No claim of a new foundational approach is warranted.

Given r supplied pieces, the certificate uses roughly `O(k²+n*k+r*n*k)` field operations and stores `O(r*n+r*k)` field values/support indices. It does not scan scalars or interpolation bases. Field construction and arithmetic are separate costs: the inherited tiny extension adapter precomputes its elements and is unsuitable for constructing the official enormous extension field. The GF729 check validates semantics, not official-scale field performance.

This extension demonstrates the value of applying existing algebraic reductions before a generic census. Its remaining conditions are substantial: the pencil must meet the code exactly, the piece witness must be available, and the strict root-count cap must certify the static list. General far-line inputs need a different reduction or a stronger static-list certificate. The Proximity Prize remains outside this finite tooling result.
