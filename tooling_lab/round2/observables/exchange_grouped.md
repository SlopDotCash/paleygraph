# Grouped exact exchange energies

`exchange_grouped.py` computes the same finite all-input L2 spectrum as
`review_exchange.py`, without iterating over the field or over input sets.
Run `python3 exchange_grouped.py`; results are in
`exchange_grouped_results.json`. This backend assumes an odd prime
`p=1 mod 4` and uses degree six in the reported experiments.

For the distinct character rows at arguments zero and one, the coordinate
types `(a,b)=(chi(-x),chi(1-x))` occur with multiplicities:

| Type | Multiplicity |
|---|---:|
| (0,1), (1,0) | 1 each |
| (1,1) | (p−5)/4 |
| (1,−1), (−1,1), (−1,−1) | (p−1)/4 each |

To derive these counts, remove the two zero sites. On the remaining `p-2`
coordinates the sums of `a`, `b`, and `ab` are all `-1`; the last equality
is the quadratic-character correlation identity. Solve the resulting four
linear equations for the four nonzero types. On the diagonal, `(1,1)` and
`(-1,-1)` each occur `(p-1)/2` times, and `(0,0)` once.

Instead of multiplying one factor per coordinate, the backend expands

```
(1+abX+aY+bZ)^m
```

for each type, using exact multinomial coefficients and retaining only
the required degree box. This yields precisely the same `K_r` overlap sums
as the row-by-row calculation.

The marked-point containment probability also has a bounded-degree form.
For degree `d`, put `u=2d-r`. With falling factorial notation `(v)_h`, it is

```
sum_(a,b=0..d-r) C(d-r,a) C(d-r,b)
  (n-j)_(u-a-b) (j)_a (j)_b / (p)_u.
```

This assigns `u` distinguished points to the common, first-only and
second-only coordinate regions of a uniform distance-`j` pair. It equals
the earlier multinomial probability but requires no factorials involving
the full field or input-set size. The Eberlein solve and rational
nonnegativity/total-energy checks remain unchanged.

For fixed degree, the polynomial coefficient work is independent of `p`
and `n` except for integer bit complexity; this prototype additionally
uses trial division to validate the supplied primes. The backend matches
**every overlap coefficient, distance correlation and exact energy** of
the earlier calculation at all six previous cases.

Three additional degree-three/five cases check the generic-degree API.
Odd-degree targets change sign under nonsquare dilation and can have
nonzero degree-two energy; the zero-level-one/two assertion is therefore
restricted to even degrees. The reported degree-six calculations are
unchanged.

| p | n | n=floor(p^(1/4))? | Exact lower-degree share, rounded |
|---:|---:|:---:|---:|
| 65,537 | 16 | Yes | 2.306125331304e−7 |
| 1,000,033 | 31 | Yes | 5.264984958370e−9 |
| 6,700,417 | 50 | Yes | 3.458075708225e−10 |
| 6,700,417 | 64 | No | 5.913835325271e−10 |

Each new exact spectrum took approximately 0.025–0.18 seconds across the
recorded runs, depending on machine load. The ratios are exact fractions
in the JSON. All degree-one and
degree-two energies are zero. The last row is explicitly an off-critical
comparison, not another critical-size test.

These results describe the uniform all-input second moment. Small global
energy does not imply a uniform bound on an individual set, and the
projection remains target-dependent. No Paley or proximity-prize proof,
new worst-case estimate, or historical novelty claim follows.
