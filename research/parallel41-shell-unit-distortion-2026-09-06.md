# Real units cannot make the shell recovery map a near-isometry

There is a uniform obstruction to one precise proposed use of real
units: they cannot make the full linear map from
[the real-subfield recovery note](parallel41-shell-arithmetic-recovery-2026-09-06.md)
preserve Euclidean lengths up to relative error tending to zero.

If `h` is any real norm-`p` generator, then its recovery map has
condition number at least

    sqrt(5)*p^(-2/N).

This holds after every real-unit multiplication. Roots-of-unity
factors do not change the singular values. In a quartic family the
lower bound tends to `sqrt(5)`. Already for dyadic `n=2N>=1024` and
`p<=n^4`, the condition number is strictly greater than two.

This refutes a near-isometry approach on the full coefficient space.
It does not rule out estimates on the restricted set of actual
cosets, a unit chosen separately for each vector, or a use of
constant distortion together with additional arithmetic information.
No class-number, unit-regulator, or uniform thin-annulus estimate
is assumed or proved here. The full goals remain open.

## 1. Exact singular values of the recovery map

As before, `n=2N>=4` is a power of two, `p` is an odd prime,
and

    R=Z[X]/(X^N+1),       bar(X)=X^(-1).

Let `h in R^+` be nonzero and satisfy `|Norm_+(h)|=p`.
The root condition needed for arithmetic recovery is irrelevant
to the following distortion lower bound; assuming it only narrows
the set of `h` to which the result applies.

On the full real coefficient space, consider

    Phi_h(F)=h*bar(F)/p.

Conjugation is a signed permutation of the coefficient basis.
Evaluation at the `N` primitive `n`-th roots, normalized by
`1/sqrt(N)`, is unitary. Hence the singular values of `Phi_h` are

    |h(zeta^j)|/p,       j odd modulo n.                            (1)

The determinant scale, meaning the geometric mean of these
singular values, is

    rho=p^(2/N-1),

because `Norm(h)=p^2`. Write `A_h=rho^(-1) Phi_h` and let
`kappa(h)` be the ratio of the largest to smallest singular value.
If

    S(h)=sum_(j=0)^(N-1) h_j^2,

Parseval gives the exact normalized mean-square identity

    (1/N)||A_h||_F^2 = S(h)/p^(4/N).                               (2)

Here `||.||_F` is the matrix Frobenius norm. The singular values of
`A_h` have geometric mean one. Therefore

    ||A_h||_op >= sqrt(S(h))*p^(-2/N),
    kappa(h) >= ||A_h||_op.                                        (3)

The second inequality uses that the smallest normalized singular
value is at most its geometric mean, namely one.

## 2. An integer trace obstruction valid for every real unit choice

**Claim.** If `h in R^+` and `|Norm_+(h)|` is an odd prime, then

    S(h)>=5.                                                       (4)

**Proof.** Put `D=N/2`. Write `h` in the integer real basis as

    h=a_0+sum_(j=1)^(D-1) a_j (X^j+X^(-j)).

Its coefficient squared norm is

    S(h)=a_0^2+2 sum_(j=1)^(D-1) a_j^2.                            (5)

Modulo two, `X^N+1=(X+1)^N`. Multiplication by `h` in the basis
of powers of `X+1` is triangular, with diagonal value `h(1)`
modulo two. This value equals `a_0` modulo two. Since the complex
norm of `h` is the odd integer `p^2`, `a_0` is odd. Thus (5) is
a positive odd integer.

If `S(h)=1`, then `h=+1` or `-1`, whose norm is one.
If `S(h)=3`, then `h=+/-1 +/- (X^j+X^(-j))` for one `j`.
Every such element is a unit of norm of absolute value one.
To check this directly, put `y=X^j`. Then

    1+y+y^(-1)=(y^3-1)/(y(y-1)).

The automorphism `X -> X^3` permutes the cyclotomic embeddings,
so `Norm(y^3-1)=Norm(y-1)`, while `Norm(y)=1`. The quotient has
norm one. The minus case follows by replacing `y` with `-y`;
this substitution has the same property because three is odd.
The denominators are nonzero: the relevant `j` satisfies
`1<=j<N/2`. An integral element of norm one is a unit by its
integer adjugate identity. Neither case can have real norm of
absolute value `p`. Hence `S(h)>=5`.

Combining (3) and (4) proves

    ||A_h||_op >=sqrt(5)*p^(-2/N),
    kappa(h) >=sqrt(5)*p^(-2/N).                                   (6)

For every real unit `epsilon`, the generator `epsilon*h` has
the same absolute real norm `p`. Thus the same lower bounds
hold for **every** choice of real unit, including an optimal
choice under any proposed balancing procedure. No estimate
for a unit lattice or its regulator enters this argument.

In a fixed quartic window, `log(p)/N ->0`, so the right side
of (6) tends to `sqrt(5)`. Consequently there cannot be real
unit choices making `A_h` approach an isometry in operator norm.

## 3. The obstruction persists under arbitrary scalar normalization

Suppose a positive scalar rescaling of the recovery map satisfied

    (1-delta)||F||_2 <= ||t Phi_h(F)||_2
                       <= (1+delta)||F||_2

for every real coefficient vector, with `0<=delta<1`.
Its condition number would be at most `(1+delta)/(1-delta)`.
Thus (6) implies, whenever its lower bound exceeds one,

    delta >=
      [sqrt(5)*p^(-2/N)-1]/[sqrt(5)*p^(-2/N)+1].                    (7)

Asymptotically in the quartic window, this gives

    delta >= (3-sqrt(5))/2-o(1),

an error bounded away from zero, even after arbitrary scalar
normalization.

There is also a simple exact size threshold. Let `n=2^m>=1024`
and `p<=n^4`. Then

    p^(4/N) <= 2^(16m/2^(m-1)) <= 2^(5/16) <5/4.

The middle exponent decreases for `m>=10`. The last inequality
is exactly `5^16>2^37`, with positive integer margin
`15148937153`. Therefore (6) gives `kappa(h)>2`, and (7)
forces `delta>1/3`. This is a conditional statement about every
real norm-`p` generator when such a generator exists; it is not
an existence theorem for those generators.

## 4. Roots of unity cannot improve the distortion

Multiplication by any power of `X` is a signed cyclic permutation
of the coefficient basis. Replacing `h` by `X^k h` therefore
left-multiplies `Phi_h` by an orthogonal matrix. Precomposing
the input with a power of `X` also results in a signed cyclic
shift of the output. Neither operation changes singular values
or the condition number.

These signed permutations commute with coefficientwise centering
modulo odd `p`. They preserve both coefficient squared norms
and cyclotomic norms. In particular the exact centering example
from the recovery note, whose norm defect changes from 641
to 178771841, has the same change after **every** root-of-unity
shift. Root-of-unity selection cannot repair that example.

## 5. Fixed norm data also permit arbitrarily bad real-unit choices

The preceding obstruction concerns all choices, including the
best possible one. A separate observation explains why a real
generator's mere availability supplies no automatic distortion
upper bound for an arbitrary representative.

For `n>=8`, the real element

    epsilon=1+X+X^(-1)

is a unit by the norm-one argument in Section 2. In the embedding
`X=exp(2 pi i/n)`, its value is `1+2 cos(2 pi/n)>2`.
For fixed nonzero `h`, the normalized operator norm of the maps
associated with `epsilon^t h` is therefore at least

    |h(exp(2 pi i/n))|*2^t/p^(2/N),

and tends to infinity as `t` grows. All these elements generate
the same real kernel and have the same absolute real norm.
This does not say that every unit choice is bad; it says that
the algebraic norm and ideal alone do not bound a representative's
Euclidean distortion.

At the concrete order-64 real generator from the recovery note,
three specified choices have the following exact coefficient
energies:

| Generator | Squared coefficient norm | Recovered element's squared coefficient norm |
|---|---:|---:|
| `h` | 19 | 13 |
| `epsilon*h` | 53 | 43 |
| `epsilon^2*h` | 305 | 293 |

The real generator norm remains `p`, and the recovered element
norm remains `1217*p`. The
[deterministic checker](../experiments/parallel41_shell_unit_distortion.py)
verifies these three specified cases and the integer threshold,
using the previous exact coefficient certificate and integer
norm arithmetic. It records the
[results](../results/parallel41_shell_unit_distortion_2026_09_06.json).
It performs no unit search or optimization.

## 6. What this excludes, and what remains

A proof cannot rely on choosing real units and roots of unity
to make the full recovery map a relative `1+o(1)` Euclidean
comparison. The integer trace obstruction rules that out
uniformly, even if the real norm generator is supplied.

This leaves open a materially different possibility: exploit
the specific geometry or arithmetic of the actual coset vectors
inside that coefficient space, rather than infer their behavior
from a near-isometry of the entire space. The bounds above do
not establish large distortion on an actual coset, do not
exclude a unit chosen separately for each vector, and do not
give a thin-annulus estimate. No unsupported claim about a
class number, regulator, or optimally balanced unit is needed
for the obstruction just proved.
