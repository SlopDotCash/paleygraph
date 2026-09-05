# Pass 24: representative bounds, the next spectral coupling, and finite subgroup classes

**The full goal remains unproved. This pass makes limited mathematical
progress, with no improvement to a worst-case Paley cancellation
exponent and no official Proximity Prize certificate.** Three parallel
workers and the root pursued four bounded investigations. The previous
[pass-23 assessment](parallel23-pass-summary-2026-09-05.md) was verified
and classified as progress before this work began.

We are closer in the sense that more actual arithmetic quantities are
controlled and the remaining gaps are explicit. These results do not
justify saying that a proof is near. In particular, a well-behaved
transformed representative does not bound its original input, and a
fixed-dimensional spectral estimate does not bound the full operator.

## 1. An unsigned B_h representative for every input

The [inversion lane](parallel24-inversion-orbit-upper-2026-09-05.md)
proves that every k-set C has an inverse D_z={1/(c−z):c∈C} which is
B_h and satisfies

    M_(2r)(D_z)
      ≤ [(2r−1)!! p²k^r+(2r−1)²p k^(2r)]/[p−Q_h(k)],
    Q_h(k)=k+(2h−2)binom(binom(k+h−1,h),2),

provided p>h and p>Q_h(k). The argument squares complete character
sums before averaging over poles and intersects that estimate with an
explicit count of poles preserving the B_h property. It assumes no
independence between the two conditions.

For fixed r≥3 and 2≤h≤floor((r+1)/2), at
k=floor(p^(1/(r+1))) this has Gaussian order p k^r with an explicit
fixed coefficient. At r=3,h=2, every C has a Sidon inverse with

    M_6(D_z)≤(20+o(1))p k³.

The quantifiers improve on merely knowing that almost all sets in the
global slice are good: now every input's inversion family contains a
good unsigned representative. The original moment transports instead
to the weights w_d=χ(d), for which

    M_(2r)(D_z,w)=M_(2r)(C)−|F_C(z)|^(2r)+k^(2r).

Those weights remain uncontrolled. The existence statement supplies
neither the uniform B_h input needed by SS-B* nor bounds on the
particular completions used to remove the signs in pass 20.

This is a new combination within the project of the earlier pole
count and the standard squared-complete-sum bound. No claim of
literature novelty is made.

## 2. Actual local persistence does not eliminate the exceptions

The [classical lane](parallel24-classical-exceptions-2026-09-05.md)
derives the exact sixth-moment drift under one deletion/insertion swap,
including the zero character entries. Its leading coefficient is
1−6(p−5)/[k(p−k)]. The remaining coefficients are nonnegative on
the thin slice. The localized insertion variance still contains the
actual tenth and eighth moments of the deleted set, which cannot be
replaced by their global averages for an adversarial input.

The evaluated shell-persistence argument must retain more than
k^(5/6) labels even at the largest exception allowed by the ordinary
Weil estimate. All sets retaining those labels occupy an exponentially
small part of the slice, far below the O(k²/p) exceptional allowance.
This argument therefore does not force that allowance to be zero.

For a good unsigned inverse, splitting its actual χ(d) weights into
majority and minority parts loses a factor of 8 in the normalized
sixth moment when the signs are balanced. Applying the ordinary Weil
bound to the minority set returns the old leading power. Strongly
unbalanced poles have only an O(p/k) upper count; no joint selection
of an unsigned-good and suitably unbalanced pole was proved. These
calculations assess specific attempted arguments, not all possible
ways of controlling the exceptional sets.

## 3. The actual second Krylov space is controlled, including leakage

For the two-anchor neighborhood C and compressed residue projection B,
let u be its normalized uniform vector and v the normalized component
of Bu perpendicular to u. The [spectral lane](parallel24-spectral-operator-2026-09-05.md)
evaluates the next coefficients using an exact recurrence involving
three Möbius pullbacks of the actual rank-three kernel. Their local
monodromy separates the cross correlations; the imported ingredients
are the pinned Katz irreducibility, local-form, and weight results.

In the basis (u,v), with Π the projection onto their span,

    uᵀBu→3/8, uᵀBv→1/8,
    vᵀBv→1/2, ||(I−Π)Bv||²→3/64.

All errors are O(p^(−1/2)). Consequently

    sup_(f∈span{u,Bu}, ||f||=1) fᵀBf
      =(7+√5)/16+O(p^(−1/2)),
    ||B restricted to span{u,Bu}||²
      =(15+√74)/64+O(p^(−1/2)).

The second statement includes the image outside that subspace. The
next coupling tends to √3/8, so deleting the remaining cross-block
still incurs a nonvanishing error. The new estimates hold on this
specified actual subspace, not on arbitrary vectors. The explicit
next residual h would require an additional two-variable estimate
hᵀS_C h=o(p^(7/2)); that estimate is unproved. Even it would be only
another fixed coefficient, with growing depth and other sectors still
to be controlled.

## 4. All quartic primes at subgroup orders 4, 8, and 16 are classified

The [subgroup lane](parallel24-subgroup-unbalanced-2026-09-05.md)
uses dyadic cyclotomic norms to classify actual opposite-free six-term
relations throughout three finite order classes.

| Subgroup order | Eligible quartic primes | Opposite-free ordered six-term count |
| ---: | ---: | --- |
| 4 | 16 | 0 throughout |
| 8 | 95 | 0 throughout |
| 16 | 579 | 480 at p=33713,37201,41521; 0 elsewhere |

At all three exceptional primes, the relations have multiplicity
pattern (4,1,1) and are product-unbalanced at every three-versus-three
partition. Thus this finite classification reaches an actual part
that the balanced-relations bound left open. It gives E_3≤15n³ for
these three orders, including their exceptions.

The norm criterion is exponential in the order, while the target
prime window is polynomial. It does not automatically extend to
growing orders. The lane also gives an exact intrinsic contribution
in every product-ratio fiber and verifies that the orbit-mass second
moment is the existing E_3 problem itself. No general aggregate saving
is extracted from that identity.

## 5. Review and verification

Each lane has a reviewer distinct from its author: the root reviewed
the classical lane, and workers cross-reviewed the other three.

| Lane | Main exact verification | Review |
| --- | --- | --- |
| Inversion | 573 actual sets; 3,570,540 coordinates in each signed/unsigned identity; 27 tuple expansions; checked integer-overflow limits | [Independent inversion review](parallel24-inversion-independent-review-2026-09-05.md) |
| Classical | 1,745 swap identities; 86 shell checks; 144 transported-sign checks; 7 larger combinatorial certificates | [Root review](parallel24-classical-independent-review-2026-09-05.md), with 858 independently enumerated increment moments |
| Spectral | 88 primes; 11,407 recurrence coordinates; 1,056 twisted cross-pullback bounds; 85 nondegenerate next-coefficient identities | [Independent spectral review](parallel24-spectral-independent-review-2026-09-05.md) |
| Subgroup | All 690 eligible prime fields; 1,770 norm determinants; 10,088 product-ratio correlations | [Independent subgroup review](parallel24-subgroup-independent-review-2026-09-05.md) |

These counts describe finite checks, not a numerical proof of an
asymptotic statement. The spectral asymptotics depend on the cited
primary results. Agent review is not human refereeing or Lean
verification. Existing formal elementary certificates are separate.

The inversion verifier was rerun after adding an explicit fixed-width
integer overflow guard. The spectral result identifies a metadata-only
post-run refresh; its recorded arithmetic was not changed. The
classical note explicitly handles its negative zero-overlap shell and
retains the p-dependent volume estimate. No mathematical conclusion
was upgraded by these reporting clarifications.

The [final artifact audit](../results/parallel24_pass_audit_2026_09_05.json)
pins current proof, verifier, result, and review bytes; checks their
dependencies and local links; preserves the prior pass's proof and
verification artifacts; and verifies the unchanged 32-entry source
manifest. Its [reproducible script](../experiments/parallel24_final_audit_2026_09_05.py)
records each check. No source version was advanced or official prize
state refreshed in this pass.

## 6. The full goal and the next inputs

The original goal remains active and unchanged. The missing inputs are:

1. A uniform character-sensitive bound on exceptional classical inputs,
   or a signed inversion/completion estimate sufficient for SS-B* at
   an unbounded sequence of fixed moment orders.
2. An aggregate upper bound for actual unbalanced subgroup relations at
   growing dyadic orders, including exceptional primes, followed by the
   full uniform square-root cancellation target.
3. Control of the remaining actual spectral operator with its nonzero
   coupling retained, at enough depth and across all required sectors.
4. Bounds on the arbitrary-word remainder fibers and bad-scalar sets
   in the [pinned official prize reduction](parallel23-prize-bridge-2026-09-05.md),
   meeting its certificate and separate spot-check obligations. Neither
   the classical formulation nor the thin-subgroup target has been
   shown here equivalent to the full official challenge.

No proof submission, external message, paid reset, model change, or
new user-owned task was made. One worker turn hit a transient selected
model capacity error and resumed successfully with the same model.
That event is not evidence of an account-wide usage limit.
