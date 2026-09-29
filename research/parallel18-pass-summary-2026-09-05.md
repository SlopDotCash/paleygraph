# Eighteenth pass: expansion is available; the centered trace remains open

We can now apply an existing SL₂ expansion theorem to every Cartesian
family arising from polynomial-sized Paley input sets. This gives an
operator power saving, but no new bound for the original character sum.
The exact trace formula includes an overlap correction. A concrete
conjugacy average shows why a small operator norm or small higher-power
traces cannot by themselves supply the missing estimate.

The [derivation](parallel18-mobius-trace-2026-09-05.md) contains the following
results, with ordinary mathematical proofs and finite exact checks.

- For G(U,V)={[[s,st−1],[1,t]]:s∈U,t∈V}, sizes m,n, every coset of
  a proper subgroup meets G in at most max(2max(m,n),120) elements.
  Point-transporter bounds and the classical subgroup classification
  supply this estimate. For m,n≥p^ε, Lyamkin's Theorem 11 therefore
  gives an operator saving p^−κ(ε), for some unspecified κ(ε)>0.
- The multiplicative energy is exactly
  E(G)=n²E₊(U)+m²E₊(V)−m²n². No low-additive-energy assumption is used.
- In the p-dimensional Weil representation, the average T over
  G(A+2,−B) satisfies tr(T)=[S(A,B)+√p|A∩B|]/(mn).
  Its (0,0) entry is 1/√p, so the generic bound |tr(T)|≤p‖T‖
  cannot establish a trace saving. The factorized Hilbert–Schmidt
  identity and diagonal Cauchy–Schwarz bound also retain the usual
  arbitrary-set threshold.
- For g₀=[[3,−1],[1,0]], the complete conjugacy average is
  V=(I+χ(5)R)/(p+χ(5)), where Rf(x)=f(−x). It has trace exactly
  one, norm 2/(p+χ(5)), and tr(Vʲ)=((p+χ(5))/2)^(1−j).
  This is a genuine probability average in the same representation.
  It is not asserted to come from a Cartesian input family and is
  not a Paley counterexample.

The source-dependent expansion application and new deductions need
independent mathematical review. The representation and expansion
theorems are classical; no novelty claim is made for those inputs.

The [verifier](../experiments/parallel18_mobius_trace_2026_09_05.py)
uses exact integer, finite-field, rational and cyclotomic arithmetic.
The final [results](../results/parallel18_mobius_trace_2026_09_05.json)
record 744,192 right-generator kernel checks, all 2,304 group character
checks at p=5,13, and 1,028 Cartesian cases: all 961 nonempty subset
pairs at p=5 and 67 cases at p=13. These include 6,704 rational and
29,672 quadratic-extension transporter rows, 1,028 energy identities,
1,028 corrected trace identities, and 1,028 factorized norm identities.
Complete conjugacy averages at p=13,17,29 verify all 1,299 matrix entries.
Finite checks do not establish the imported asymptotic theorem or
the missing uniform Paley estimate. The [final audit](../results/parallel18_pass_audit_2026_09_05.json)
pins these artifacts and preserves the preceding pass.

The new proof inputs are primary HTML: Thomas's arXiv v3 character
formulas are archived; Lyamkin's norm definitions, subgroup classification
and expansion statement were reviewed live. Direct Lyamkin archival
requests returned 403, recorded as a metadata-only manifest entry.
No PDF was newly rendered or visually reviewed in this pass.

The research is more precise about where this approach stops, but
there is no evidence that the full proof is close. Classical two-set
cancellation, the uniform subgroup square-root target and exceptional
primes, the spectral edge, and the quantitative bridge to the official
Reed–Solomon prize remain unproved. The next argument must control the
centered trace using additional Cartesian character structure; asking
for that estimate without further reduction just restates Paley.

All three parallel workers remain terminal at usage limits in the live
check during this pass. Root completed this work locally. The full goal
remains active and unachieved.
