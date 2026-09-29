# Independent exchange-spectrum review

**Verdict: no actionable mathematical or implementation error found.**
The original generator was not modified. `review_exchange.py` adds an
independent exact computation of the **all-input** L2 spectrum, replacing
the larger cohorts' sampled energy estimates with exact rational values.

Run `python3 review_exchange.py`. The saved result records the reviewed
generator hash. The review implementation imports no algebra from the
generator. Its approximately two-second runtime is machine dependent.

## Checked claims

1. **Degree.** For `T_6(C)=sum_(|S|=6,S subset C) E(S)`, with
   `E(S)=sum_y prod_(x in S) chi(y-x)`, membership indicators give a
   multilinear polynomial of degree six. On the uniform size-`n` slice,
   this lies in the Johnson harmonic degrees zero through six. The
   generator's seven polynomial projectors are valid on this subspace.
   They are not a complete set of orthogonal projectors for an arbitrary
   higher-degree function when `n>6`.
2. **Distance averages.** A kernel set with `r` points inside `C` and
   `6-r` outside is contained in a uniformly chosen distance-`j` neighbor
   with probability
   `C(n-j,r)/C(n,r) * C(j,6-r)/C(p-n,6-r)`.
   The sign-count histogram includes the unique zero correctly. Dividing
   the all-field polynomial `(1-t^2)^((p-1)/2)` by the inside generating
   polynomial gives the correct outside coefficients.
3. **Walk spectrum.** The normalized one-swap eigenvalue is
   `1-d*(p-d+1)/(n*(p-n))`. The radial walk probabilities for increasing
   and decreasing distance are `(n-j)*(p-n-j)/(n*(p-n))` and
   `j^2/(n*(p-n))`. The generator's finite polynomial interpolation uses
   distinct eigenvalues in its admitted range.
4. **Normalized sampling.** The even-degree target is invariant under
   every affine map of the field. Johnson projections commute with all
   coordinate permutations, so their products are affine invariant too.
   Given an arbitrary set and an ordered pair of its distinct points,
   the unique affine map sending that pair to `(0,1)` gives a set
   containing `(0,1)`. Each such normalized set has exactly `p*(p-1)`
   preimages among set/ordered-pair choices. Thus uniform normalized-set
   averages equal all-set averages for the measured invariant quantities.
   This would not justify arbitrary biased normalization procedures or
   non-invariant observables.
5. **Interpretation.** The projection uses the target itself. It diagnoses
   which harmonic degrees contain target variation; it is not a target-free
   predictor. Sampled component squares need not sum to the sampled
   centered square because cross terms need not vanish in a finite cohort.
   A sampled top-energy ratio above one is therefore not a projector error.

Known reference machinery: [Filmus, *Orthogonal basis for functions over a
slice of the Boolean hypercube*](https://arxiv.org/abs/1406.0142) and the
Johnson/Eberlein eigenvalue formula used, for example, in [Laurent's
semidefinite analysis, equation (17)](https://homepages.cwi.nl/~monique/files/laurent2.pdf).
The review does not claim new slice-harmonic theory.

## Independent exact L2 calculation

Define

```
K_r = sum_(|S|=|T|=6, |S cap T|=r) E(S) E(T).
```

For fixed character-row arguments `y,z`, put `a_x=chi(y-x)` and
`b_x=chi(z-x)`. Their contribution to `K_r` is

```
[X^r Y^(6-r) Z^(6-r)] prod_x (1 + a_x*b_x*X + a_x*Y + b_x*Z).
```

The three nonconstant choices assign a coordinate to `S cap T`, `S\T`,
or `T\S`, respectively; their signs are exactly the two kernel products.
Translation reduces diagonal row pairs to one case. An affine change of
variable reduces every distinct row pair to `(0,1)`; the total character
degree is twelve, so its scaling sign is one. Hence `K_r` is `p` times
the diagonal coefficient plus `p*(p-1)` times the distinct-row coefficient.
The truncated coefficient DP has at most 75 states for each `r`, so it
requires no enumeration of input sets.

For fixed `S,T`, let `u=12-r`. A uniform ordered pair `(C,D)` of `n`-sets
at distance `j` contains `S` in `C` and `T` in `D` with probability

```
sum_(a,b=0..6-r) C(6-r,a) C(6-r,b)
  Multinomial(p-u;
    n-j-u+a+b, j-a, j-b, p-n-j)
 / ( C(p,n) C(n,j) C(p-n,j) ).
```

Here `a,b` count points forced into `C\D` and `D\C`; negative multinomial
arguments contribute zero. Multiplying by `K_r` and summing gives the exact
distance correlation `H_j=E[T_6(C)*(A_j T_6)(C)]`.

The review then solves the independent Eberlein system

```
H_j = sum_d theta_d(j) Energy_d,

theta_d(j) = sum_t (-1)^t C(d,t) C(n-d,j-t) C(p-n-d,j-t)
             / (C(n,j) C(p-n,j)).
```

This avoids the generator's walk-polynomial projector construction. Every
energy is nonnegative; the degree-zero energy equals the squared mean;
all energies sum exactly to the second moment.

## Exact results

| p | n | Degree-six share of centered L2 | Degrees 1–5 share |
|---:|---:|---:|---:|
| 13 | 6 | 0.511380145278 | 0.488619854722 |
| 17 | 6 | 0.838498362308 | 0.161501637692 |
| 61 | 6 | 0.994772183132 | 0.005227816868 |
| 1297 | 6 | 0.999990979174 | 0.000009020826 |
| 2437 | 7 | 0.999992371929 | 0.000007628071 |
| 4129 | 8 | 0.999994697920 | 0.000005302080 |

All entries are exact rational quantities internally; the table rounds
them. Degrees one and two have exactly zero energy in every audited row.
The p=61 exact value also shows why the 512-sample ratio should not be
treated as a global energy measurement.

Verification includes every one of the **2,944,656 ordered input pairs**
at `p=13,n=6`, all 108,900 pairs of a separate `p=11,n=4,degree=2` constant
null control, independent pointwise projection of three inputs per case,
and exact agreement with the generator's two exhaustive Gram diagonals
and its saved p=13 pointwise examples. The larger-row answers are exact
computations using the proved counting formula, not exhaustive input
enumeration or Lean-checked theorems.

Small global L2 energy is not a worst-case certificate. A lower component
can still matter on particular sets, especially where the centered target
is small. No uniform pointwise bound or conjecture conclusion follows.
