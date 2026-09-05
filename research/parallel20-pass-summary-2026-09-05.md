# Twentieth pass: inversion strengthens the restricted moment reduction

The restricted moment criterion now allows B_h relation order
h≤⌊(r+1)/2⌋ for moment order 2r. This replaces pass 19's sufficient
condition h/r→0. The new transfer keeps every input element and has
no partition remainder. **The uniform moment estimate remains unproved.**

The [derivation](parallel20-inversion-moments-2026-09-05.md) establishes
three ingredients.

- For a set C of size k, at most
  Q_h(k)=k+(2h−2) binom(binom(k+h−1,h),2) poles z either belong
  to C or make {1/(c−z):c∈C} fail the B_h property. Each pair
  of different h-term multisets gives a nonzero polynomial of degree
  at most 2h−2. The multiplicities remain nonzero because p>h.
- Inversion introduces signs w=χ(c−z), with the exact moment identity
  M_(2r)(D,w)=M_(2r)(C)−|F_C(z)|^(2r)+k^(2r). Its correction
  is favorable for an upper bound, but cannot be omitted.
- A B_h completion from size k to 2k expresses any signed size-k
  indicator as two averages of unsigned size-k indicators plus a
  coefficient of absolute value at most one times another size-k
  indicator. The resulting moment cost is at most 3^(2r).

Under the explicit finite conditions p>Q_h(k) and p>f_h(2k−1),
a uniform unsigned moment bound on B_h sets of size k therefore
implies the same bound on all sets of size k, multiplied by 3^(2r).
Restriction gives the reverse direction. These uniform moment
assertions are equivalent up to that constant in the stated range.

For k=⌊p^(1/(r+1))⌋, the conditions hold for sufficiently large
primes whenever 2h≤r+1, including equality. An explicit sufficient
threshold is p≥K_h^(r+1), with
K_h=max(2,4h(h−1),(h+1)2^(2h−1)+1). Each order is fixed before
p grows; the constants are not uniform estimates at growing order.

Thus an unproved sixth-moment estimate on Sidon sets at size p^(1/4)
would imply the unrestricted sixth-moment estimate and Paley above
ε=1/4. A sequence of the restricted bounds, with the allowed power
loss tending to zero, would imply the full classical conjecture.
No such character-moment bound was proved here. The Sidon fourth
moment at size p^(1/3) is outside the range established by this method.

The [exact verifier](../experiments/parallel20_inversion_moments_2026_09_05.py)
and [results](../results/parallel20_inversion_moments_2026_09_05.json)
check 56,256 nonzero relation polynomials, 646,347 rational/polynomial
equivalences, and 1,319 complete sets of poles. They also check
112,762 pointwise inversion rows, 20,636 corrected moment identities,
102 sign patterns, and 408 signed moment transfers. Thirty-two
complete examples give 128 transfers from original moments to unsigned
B_h completion moments. Twenty-two rational polynomial certificates
cover all k beyond the explicit thresholds for h=2,…,12.
These finite checks do not prove the missing uniform moment estimate.
The [final audit](../results/parallel20_pass_audit_2026_09_05.json)
pins the artifacts and preserves pass 19.

No new external theorem or source is imported; the source manifest
is unchanged. No PDF was newly rendered or visually reviewed.
Independent mathematical review is outstanding. The full classical
conjecture, subgroup square-root target, spectral edge and official
Reed–Solomon prize bridge remain unproved. Root worked locally; the
workers were stopped at account limits in the most recent live check.
