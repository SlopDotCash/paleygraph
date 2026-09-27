# Translation boxes and bounded erasure decoding through a scalar quotient

This round studies the finite centered digit section used in the existing Paley tooling. It gives exact completion algorithms and certificates for selected erasure patterns. It does not estimate the relevant norm uniformly or prove either prize. Lattice basis reduction, coordinate rounding, carry identities and erasure decoding have extensive prior art; historical novelty of this application remains unestablished.

## 1. The arithmetic section and its projection

Work in R=Z[X]/(X^N+1), with N>=4 a power of two, p an odd prime, g of order2N and u=g^-1 mod p. Let nonzero f satisfy f(u)=0 mod p. Write m=(p-1)/2 and F_a[j]=center(a*g^j mod p) in[-m,m]. The actual digit word is E(a)=f*F_a/p, with coefficients bounded by B=floor(m*||f||1/p).

The existing inverse certificate supplies A_f and k>0 such that f*A_f=kp. Consequently F_a=A_f*E(a)/k. Let w be the coefficient-zero row of multiplication by A_f:

`w = (A_f[0], -A_f[N-1], ..., -A_f[1])`.

Every actual digit word D therefore satisfies w.D divisible by k, and its scalar is

`a = (w.D/k) mod p`.

These facts hold even when round19's field calibration beta vanishes or p divides k. No good-calibration assumption is made in this decoder. The exact inverse identity, rather than an assumed LLL quality guarantee, is the arithmetic foundation.

## 2. Residual intersections become translated centering boxes

Fix a coordinate set U and two assignments P,Q on U. A shared completion supplies the same digits outside U; hence the difference of its two actual words is H, equal to Q-P on U and zero elsewhere. Define

`V = A_f*H/k`.

If V has a nonintegral coefficient, the two residual sets are disjoint. If integral, put delta=V[0] mod p. If any V[j] differs from delta*g^j mod p, they are again disjoint. Otherwise an actual word E(a) with assignment P has a shared completion with Q exactly when

`max(-m,-m-V[j]) <= F_a[j] <= min(m,m-V[j])` for every j.

Indeed these inequalities say that F_a+V is centered. Its scalar congruences are those of a+delta, so F_a+V=F_(a+delta), and applying f/p gives E(a+delta)=E(a)+H. The converse follows by subtracting the two actual centered vectors. Empty coordinate intervals certify disjointness without examining a scalar fibre.

There is also an exact carry description. Set C=(F_delta-V)/p. In the feasible-box case its coordinates belong to{-1,0,1}, and

`F_(a+delta)=F_a+F_delta-p*C`,
`E(a+delta)-E(a)=E(delta)-f*C=H`.

The complete small audit covers every pair of nonempty prefixes at every depth, plus one empty representative at each positive depth. Independent rational inverses and direct suffix sets check364029 comparisons. The audit includes squared relations and scalar p,p^2 controls, so divisibility alone is not silently substituted for the scalar-congruence check.

This reduction by itself does not enumerate a general fibre. The next step uses its necessary arithmetic congruence to construct a complete candidate set for erasures.

## 3. A congruence lattice for erased coordinates

Let U be the erased coordinates, with d=|U|; write the unknown digit vector as x and the known digits as y. Let w_U,w_K be the corresponding parts of w. Every completion satisfies

`w_U.x = -w_K.y mod k`.

Let h=gcd(k,w_U[0],...,w_U[d-1]). If the right side is not divisible by h, there is no completion. Otherwise a Bezout identity gives one integral solution x0. Every integral solution is

`x = x0 + z*L`, with z in Z^d,

where the rows of L form a basis of the congruence kernel. The kernel has index k/h in Z^d. The certificate checks that each basis row lies in the kernel and that |det L|=k/h; together these prove that L spans the entire kernel. It also checks an exact rational inverse T=L^-1.

LLL is used only to choose a useful basis. An initial SymPy1.14 LLL call failed an internal size-reduction assertion on an eight-dimensional discovery input. The prototype switched to the already installed fpylll0.6.4. Its output is accepted only after exact determinant, inverse and congruence checks. No floating-point lattice claim enters the proof.

## 4. Remove directions invisible to the recovered scalar

For each basis row L_j define its scalar step

`s_j = (w_U.L_j/k) mod p`.

This quotient is integral because L_j lies in the congruence kernel. Let a0=(w_K.y+w_U.x0)/k mod p. For any possible completion,

`a = a0 + sum_j z_j*s_j mod p`.

If s_j=0, z_j cannot change the scalar. It need not be enumerated. This is stronger than merely noticing a short basis row: the exact scalar-step identity justifies discarding it. For the main consecutive erasure patterns, many such directions are ordinary multiples of f, but the algorithm requires only zero scalar step, so it also handles non-unit and degenerate controls.

This avoids a table with k states. The large cofactor k=9985208709332560769028097 instead supplies a strong congruence. Large-integer preparation and basis inversion remain real costs; the algorithm does not claim that they are free.

## 5. First coordinate box: bounded digits

Since z=(x-x0)T and each unknown digit satisfies |x_i|<=B,

`|z_j + (x0*T)_j| <= r_digit[j] = B*sum_i |T[i,j]|`.

This gives an exact integer interval for every z_j. An empty interval rules out a completion. Enumerating only the intervals with nonzero s_j gives a finite superset of possible scalars. Re-encode each distinct scalar and retain exactly those words matching every known digit.

Completeness follows because every actual word supplies an integral z inside every certified interval, and its scalar is included in the visible-coordinate sum. Soundness follows from final actual re-encoding and equality with the known digits. The chosen x0 or an intermediate lattice vector need not itself be centered or have all coefficients in the digit box.

The number of candidate combinations is bounded uniformly in y by

`product over visible j of (floor(2*r_digit[j])+1)`.

An exceeded work limit is an incomplete search, with no completion count attached. The first version proved unique completion for15 initial erased coordinates at the largest input, but its32-coordinate candidate cap was110612791296. Four of six24-coordinate examples initially exceeded20000 candidates; increasing the limit to100000 completed all six. Both stages remain recorded.

## 6. Stronger box: retain the original centering cancellations

The first box replaces the actual encoding equations with independent digit bounds. That triangle inequality loses cancellation. Let M be the integer negacyclic multiplication matrix for f. For an actual word,

`x = M_U*F_a/p`.

For inverse-basis column T_j define the exact rational vector

`q_j = M_U^T*T_j`.

Then x.T_j=q_j.F_a/p, so the original centered bound gives

`|z_j + (x0*T)_j| <= r_center[j] = (m/p)*||q_j||1`.

The implemented radius is the smaller of r_digit[j] and r_center[j]. Every q_j coefficient and both radius calculations are checked by a separate rational-matrix implementation. This is a different certificate from summing bounds on the already compressed digits; cancellations inside q_j are retained exactly.

For the saved p=2013265921,N64, nine-term relation f, the certified uniform caps change as follows:

| Consecutive erased coordinates | Independent digit bounds | Centered pullback bounds |
|---:|---:|---:|
|15|1|1|
|16|2|1|
|20|8192|1|
|24|419904|4|
|32|110612791296|41472|

The caps are on candidate combinations, not observed list sizes or wall-clock times. They hold for every bounded assignment to the remaining coordinates for these specific certificates. All32 consecutive erasures can therefore be handled within the50000-candidate budget. This is a finite full-input guarantee at the specified prime and relation, not a statement uniform over all p,N,f or erasure patterns.

## 7. Universal uniqueness, cyclic transport, and state lower bounds

If 2*r[j]<1 for every visible direction, any two actual completions of the same known digits must have the same scalar. To see this, their unknown difference H lies in the congruence kernel, say H=t*L. Applying either difference bound gives |t_j|<=2*r[j]. Thus every visible integral t_j is zero. The scalar difference sum t_j*s_j vanishes, so the scalars and their full encodings coincide. This proves uniqueness for all known assignments; it is not inferred from successful samples.

At20 consecutive erasures the improved certificate satisfies this strict condition. Negacyclic multiplication by X commutes with f and maps E(a) to E(a/g). Hence a cyclic block can be moved to the initial block with the appropriate sign changes. The certificate gives unique recovery after **any20 consecutive cyclic erasures**, using the other44 coordinates. At32 consecutive cyclic erasures the same transport preserves the41472-candidate cap, with no uniqueness assertion. The code independently tests every cyclic start.

Equivalently, different nonempty20-digit prefixes have disjoint complete residual languages. A complete census of scalars0 through99999 finds100000 distinct such prefixes. Therefore any deterministic layered reader in the specified natural order needs at least100000 nonempty states at that cut. The erasure theorem certifies all4999950000 pairwise distinctions at once; no cross-splice matrix is built. A separate direct-power/dense-matrix census checks every prefix. This is a finite model-specific lower bound, not an all-orders or asymptotic lower bound.

## 8. Limits exposed by controls

The same16 held-out scalars and bounded edits were tested under three32-coordinate erasure geometries. Consecutive erasures completed. Alternating and seeded random erasure sets produced enormous certified boxes and remained explicitly incomplete at50000 candidates. Thus the tested algorithm's useful bound depends strongly on the coordinate geometry.

The small complete controls include many queries with multiple actual completions, including an entirely erased17-word code. They prevent replacing complete decoding with a uniqueness assumption. Degenerate controls use f squared, f=p and f=p^2, and exercise non-unit syndrome divisibility.

A second control is D=f at the main input. It lies inside the digit-height box, satisfies the lattice syndrome, and shares all known suffix digits with zero when the first15 coordinates are erased. Its inverse lift has coefficient p, outside the centered interval. The scalar projection identifies it with scalar0, and actual re-encoding retains only the zero word. Thus syndrome feasibility alone still does not establish realizability.

The broad unresolved task is to control candidate growth for general erasure geometry and to attach useful arithmetic aggregates to these exact completion mechanisms. Prefix emptiness and arbitrary residual equivalence remain unsolved in general. The proximity decoder's supplied-space scope is unchanged. No norm estimate for every scalar, spectral-edge estimate, general list-decoding theorem, or prize proof follows from these finite certificates.
