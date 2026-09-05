# Independent review of exact spectral certificates

Review date: 2026-09-05. Source hashes are recorded in
`spectral_exact_review.json`. The reviewed C++ source is
`49b8ffd38fff67bcdf3b4873e9f72046c5efc962f01a1c247a7b7c0f41f77bee`;
the Python directional-bound checker is
`0bbd889432133b3b1228dbc2813ba3f9ec7ce0a0c7c5d21cf96ab658a7561bc9`.

No unresolved mathematical error was found within the stated input range.
One output ambiguity was reported and corrected: skipped autocorrelation
calculations now have `autocorrelation_computed=false`. In that case the
legacy values maximum=-1, shift=0, count=0 are sentinels; count=0 does not
assert that no translation has large overlap. When computed, the count
includes the identity shift zero. The updated README explains both cases.

## Independent execution

`review_spectral_exact.py` freshly compiles the inspected backend into this
review directory. Its expected convolution uses Python arbitrary-precision
integer multiplication with coefficient packing, not the candidate's NTT,
NumPy, or an FFT. Adding constant offsets makes coefficients nonnegative;
the radix exceeds the largest possible coefficient. Ordinary integer
multiplication then computes the polynomial product without carries
between coefficients. Cyclic folding and the exact constant offset
correction recover the integer convolution. This is established Kronecker
substitution, used here as an independent arithmetic oracle.

All 11 witnesses passed: both signs at p=101,401,1009,4001,65537, plus the
lower-coefficient p=65537 negative witness. Every small witness row was
also checked by direct summation (2,746 rows altogether). For each of the
three p=65537 witnesses, 32 rows were checked directly in addition to the
full packed-product quadratic-form comparison.

The 512-coefficient witness was independently checked at **all 65,537
translations**, giving:

```
norm_squared = 343365276
sum_entries = 50
signed_quadratic_form = -76160140876
maximum nonzero absolute autocorrelation = 12426901
first maximizing shift = 17303
number reaching 3/5, including zero = 1
```

The overlap ratio is strictly below 181/5000=0.03620, since
`181*343365276 - 5000*12426901 = 14609956 > 0`. The exact directional
Rayleigh lower bound 866400/1000000 was independently accepted by an
integer square-margin computation. That rational number is above
sqrt(3)/2. These facts support the stated finite witness obstruction,
without an infinite-family or universal inverse-theorem claim.

`review_million_samples.py` additionally checks the new p=1,000,033,
128-coefficient witness. A local copy of the C++ source adds diagnostic
printing only, exposing five character-convolution rows and eight
autocorrelation entries. All match independent direct integer sums,
including the reported maximizing shift 1 with numerator 2,189,980 and
norm squared 195,367,076. The full freshly compiled C++ output matches the
saved artifact, and the 867000/1000000 lower bound has positive exact
square margin. This is **sampled independent entry verification** plus
full candidate replay, not an independent full million-point convolution.
The large numerical eigensolver was not restarted or modified.

## Mathematical audit

For degree at most p-1 inputs, the linear product has degree at most
2p-2. Zero padding to length at least 2p-1 prevents NTT aliasing before
the intended cyclic fold. Modulo X^p-1, coefficient x is the sum of linear
coefficients x and x+p, with the latter included only if at most 2p-2.
The implemented wrap condition is correct, including x=p-1.

The modular prime is trial-divided, its predecessor is 2^23*7*17, and
the root is tested against each distinct prime factor. This establishes
the needed primitive roots for the supported power-of-two transform
lengths. All reviewed builds retain assertions.

The character convolution obeys |c(x)| <= L, where L is the input l1
norm. The checked condition 2L < 998244353 gives a unique centered integer
recovery. The quadratic form is then accumulated as sum z(x)c(x). For
autocorrelation, Cauchy-Schwarz gives |a(t)| <= n, where n is norm squared;
the additional condition 2n < modulus is therefore sufficient. Skipping
that second calculation does not invalidate the first calculation.

For p <= 1,000,033 and coefficients of absolute value at most 1024:

| Quantity | Absolute upper bound |
|---|---:|
| Neighborhood size | 250,007 |
| L | 256,007,168 |
| n | 262,151,340,032 |
| Quadratic-form absolute value, bounded by L squared | 65,539,670,067,380,224 |
| Product of two modular residues | 996,491,786,299,899,904 |
| Signed int64 maximum | 9,223,372,036,854,775,807 |

Thus coefficient squares, all partial quadratic-form sums, modular
products, and the factors used in correlation thresholds fit int64.
The padded length is at most 2^21, below the supported 2^23. This review
does not certify arbitrary uint32 inputs beyond the stated range.

For positive b and n, the directional comparison reduces exactly to
X*sqrt(p)>Y, where X=b*sign*A and Y=a*p*n+b*sign*sum(z)^2.
Opposite signs decide the result immediately. For two positive quantities
squaring preserves order; for two negative quantities it reverses order.
The zero and equality cases in `above` correctly preserve strictness.
The independent script checked 15,680 bounded combinations, including
positive and negative numerators, zero, perfect-square p, and equality.
Every stored witness lower bound was also checked with Python integers.
Floating values propose bounds and provide display consistency only.

## Limits

These are exact finite Rayleigh lower certificates for explicit integer
vectors, with exact transport data where the recovery precondition holds.
They do not certify completeness of a numerical Ritz span, a spectral
upper bound, a persistent edge excess along infinitely many primes, or
the nonexistence of some other structured witness. Neither the NTT nor
the independent integer-packing oracle is claimed as invented mathematics.
