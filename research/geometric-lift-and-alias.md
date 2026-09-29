# A uniform moment bound for the geometric lift

The main prime-field conjecture remains open. This note proves a uniform
moment estimate for an auxiliary geometric cycle, then identifies the exact
arithmetic error that must be controlled when passing to a prime divisor.
The estimate is proved below in ordinary mathematics, not in Lean. No novelty
claim is made relative to the literature.

## The auxiliary model and theorem

Fix integers `a≥2`, `N≥2`, and put

\[
Q=a^N+1,\qquad n=2N,\qquad
H_Q=\{\pm a^j\pmod Q:0\le j<N\}.
\]

The order of `a` modulo `Q` is `2N`. Indeed, powers below `N` lie strictly
between 1 and `Q`, while `a^{N+j}≡-a^j` cannot equal 1 for `0≤j<N`.
Thus all `n` representatives are distinct units. Define the real period

\[
F_Q(b)=\sum_{h\in H_Q}e^{2\pi i bh/Q}
       =2\sum_{j=0}^{N-1}\cos(2\pi b a^j/Q).
\]

**Theorem.** For every integer `r≥1`,

\[
\boxed{\frac1{Q-1}\sum_{b=1}^{Q-1}|F_Q(b)|^{2r}
       \le (160rn)^r.}
\tag{A}
\]

The constant is independent of `a,N,r`, and `Q` need not be prime.

### Proof: changing one base-a digit has bounded effect

Let `B` be uniform on `{0,…,Q-2}={0,…,a^N-1}`. Its `N` base-`a` digits
are independent and uniform. Changing digit `j` changes `B` by `δa^j`,
where `|δ|≤a-1`. Multiplication by `a^j` permutes `H_Q`, so

\[
|F_Q(B+\delta a^j)-F_Q(B)|
\le\sum_{h\in H_Q}|e^{2\pi i\delta h/Q}-1|
\le\frac{4\pi|\delta|}{Q}\sum_{k=0}^{N-1}a^k
<4\pi.
\]

Here the first inequality is the triangle inequality after factoring out
unit complex phases, and the second uses `|e^{ix}-1|≤|x|`, choosing the
signed representatives `±a^k`. This bound holds for every digit change,
including changes of more than one in the digit's value.

Write `Z=F_Q(B)`. Reveal the independent digits one at a time and let
`D_j` be the increments of the Doob martingale for `Z`. Conditional on
the preceding digits, each `D_j` has mean zero and range length at most
`4π`: conditional expectations for two values of the next digit differ
by at most the same bounded difference, after averaging the future digits.

For completeness, a mean-zero variable of range length `c` satisfies
`E exp(tY)≤exp(t²c²/8)`. To see this, the second derivative of
`log E exp(tY)` is the variance under exponential tilting. Any variable
in an interval of length `c` has variance at most `c²/4`, including under
that tilted distribution. Integrating this second-derivative bound twice,
using value and first derivative zero at `t=0`, proves the assertion.

Applying this conditional bound successively gives

\[
\mathbb E e^{t(Z-\mathbb EZ)}\le e^{\pi^2nt^2},\qquad
\Pr(|Z-\mathbb EZ|\ge s)\le2e^{-s^2/(4\pi^2n)}.
\]

Integration of the tail yields

\[
\mathbb E|Z-\mathbb EZ|^{2r}\le2(4\pi^2n)^r r!.
\]

Character orthogonality over **all** residues modulo `Q` gives
`Σ_b F_Q(b)=0`; hence `EZ=-F_Q(-1)/(Q-1)` and `|EZ|≤n/(Q-1)≤1`.
Since `n≥4`, the elementary inequality for the power of a sum gives

\[
\mathbb E|Z|^{2r}
\le(16\pi^2n)^r r!+2^{2r-1}\left(\frac n{Q-1}\right)^{2r}
\le(16\pi^2rn)^r+n^r
\le((16\pi^2+1)rn)^r
<(160rn)^r.
\]

Finally, replace the digit interval `{0,…,Q-2}` by the nonzero residues.
This removes `b=0`, whose period is `n`, and adds `b=Q-1`, whose period
has absolute value at most `n`. The sum of even powers cannot increase.
This proves (A). The last numerical inequality follows already from `π<22/7`.

## Exact passage to the prime subgroup

Let `p` be an odd prime and let `a mod p` generate a subgroup `H` of
even order `n=2N`. Then `a^N≡-1 mod p`, so `p|Q`. Put `d=Q/p`.
Reduction maps the `n` members of `H_Q` bijectively onto `H`, and

\[
F_Q(db)=\eta_b(H).
\]

For a residue `s mod Q`, let `C_k^{(Q)}(s)` count ordered `k`-term sums
from `H_Q` equal to `s`. The prime additive energy satisfies the exact identity

\[
E_r(H)=C_{2r}^{(Q)}(0)+\sum_{t=1}^{d-1}C_{2r}^{(Q)}(tp).
\]

Separate the expected uniform contribution by defining

\[
D_r=\sum_{t=1}^{d-1}C_{2r}^{(Q)}(tp)
       -\left(\frac1p-\frac1Q\right)n^{2r}.
\tag{D}
\]

Then

\[
\frac{pE_r(H)-n^{2r}}n
=\frac pn\left[
 C_{2r}^{(Q)}(0)-\frac{n^{2r}}Q+D_r\right].
\tag{T}
\]

By orthogonality and (A), the first term in brackets is bounded above
by `(160rn)^r`. Thus a uniform upper bound

\[
D_r\le(Krn)^r
\tag{U}
\]

at `r=ceil(ln((p-1)/n))` would imply the sufficient estimate (SG) from
`subgroup-target.md`, for example with constant `2(160+K)`.
Indeed `x^r+y^r≤(x+y)^r` for nonnegative `x,y`, and `p/(p-1)≤2≤2^r`.

The converse holds up to constants as well. Write

\[
B_r=C_{2r}^{(Q)}(0)-n^{2r}/Q
    =\frac1Q\sum_{b=1}^{Q-1}|F_Q(b)|^{2r}\ge0,
\qquad Q_r(H)=\frac{pE_r(H)-n^{2r}}n.
\]

Equation (T) gives `D_r=(n/p)Q_r(H)-B_r`. If (SG) holds with
`m=(p-1)/n` and constant `K`, then

\[
D_r\le\frac np Q_r(H)
\le\frac{p-1}{p}(Krn)^r\le(Krn)^r.
\]

Consequently **(U) and (SG) are equivalent up to absolute constants** at
the chosen moment depth. The auxiliary theorem isolates a known component,
but the remaining discrepancy estimate retains the strength of the original
moment problem. It has not been shown to be an easier estimate.

**(U) is unproved.** It concerns the signed excess in a specific congruence
class, after subtracting its uniform mass. The model estimate (A) supplies
no automatic bound for that excess. Simply conditioning the uniform
frequency on being divisible by `d` loses the factor `d`, which is generally
exponentially large in `n`. The lift changes the description of the remaining
arithmetic problem; it does not solve the prime-field estimate.

The construction applies to every even-order subgroup, including every
dyadic subgroup in the working target. The value of `D_r` depends on the
chosen integer generator and lift; it is not an intrinsic subgroup invariant.

## Exact carry computation of the model and discrepancy

The new `experiments/dyadic_carries.py` counts these words without allocating
an array of size `Q`. At digit `j`, write the signed number of occurrences
of `a^j` as `e_j=u_j-v_j`. If the target has base-`a` digits `t_j`, use carries

\[
a c_{j+1}=c_j+e_j-t_j.
\]

Multiplying by `a^j` and telescoping gives

\[
\sum_{j=0}^{N-1}(e_j-t_j)a^j=a^N c_N-c_0.
\]

The boundary condition `c_N=-c_0` is therefore exactly congruence modulo
`a^N+1`. For a word of length at most `k`, all possible initial carries lie
between `-floor(k/a)-1` and `floor(k/a)+1`, so a finite enumeration is complete.
At a step using `ell` new letters, of which `u` are positive, the weight is
`binom(s+ell,ell) binom(ell,u)`, where `s` letters were already assigned.
This counts all interleavings of the ordered letters, including repetitions.

The implementation checks every target residue through order 8 against
independent dense convolutions at seven small `(a,N)` pairs. For the main
example `(a,N,p)=(2,32,6700417)`, the cofactor is 641. Its nonzero multipliers
split into ten doubling orbits of size 64; one carry calculation per orbit
therefore accounts for every additional prime relation.

The resulting counts reconstruct **all twelve** previously computed prime
energies through order 24, using a different algorithm. The saved output is
`results/dyadic_carries.json`, including a hash of the compared quotient data.

## What the calculation resolves in the current example

Here `Q=2³²+1=641p`, `n=64`. Intrinsic relations modulo `Q` already give
the entire fourth- and sixth-moment counts:

| Moment order k | Intrinsic count C_k^(Q)(0) | Additional prime relations |
|---|---:|---:|
| 2 | 64 | 0 |
| 4 | 12864 | 0 |
| 6 | 4816960 | 0 |
| 8 | 2855358016 | 12902400 |

All additional counts below order 8 vanish. An explicit shortest relation is

\[
1-2^7-2^9+2^{14}-2^{17}+2^{19}-2^{21}+2^{23}=6700417.
\]

It vanishes modulo `p` and does not vanish modulo `Q`. At order 8 there are
`12902400=64·5·8!` ordered additional words. Thus the known fourth-moment
Gaussian failure comes from the intrinsic doubling relations, before any
extra reduction modulo this prime. The new terms first enter at order 8.

For generator 2, `D_r<0` at every tested order `r=1,…,12`. That sign is
**not** a general fact about lifts: generator 8 produces the same subgroup
but the modulus `Q'=8³²+1`. Its intrinsic fourth count is 12096, and

\[
D_2'=12864-12096-\frac{64^4}{6700417}+\frac{64^4}{8^{32}+1}>765>0.
\]

Both signs are exact rational comparisons. Consequently a proposed proof
cannot remove (D) by simply declaring it nonpositive for every generator.

## Relation to the lacunary literature

[Aistleitner, Frühwirth, Hauke, and Manskova (2025)](https://arxiv.org/html/2502.20930v2)
prove a quadratic leading term for the moment-generating function of lacunary
cosine sums, with a cubic error, under Lebesgue sampling. Their doubling
example has nonzero cubic and positive quartic corrections. The theorem's
sampling measure differs from the prime grid. The proof of (A) above instead
uses independent digits on a specific finite geometric cycle and is given
in full here. Neither result alone controls the prime-reduction error (D).
