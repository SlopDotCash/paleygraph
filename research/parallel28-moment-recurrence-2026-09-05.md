# Combining the third moment with an invariant-set energy estimate

**The Paley conjecture and the Proximity Prize remain unproved.**
This note derives a stronger unconditional bound for the project's thin
multiplicative-subgroup problem from existing literature inputs. This is
an ordinary mathematical argument, not a Lean proof or a claim of a new
result in the literature. Separate-author review of this combination is
still outstanding.

Let H be a multiplicative subgroup of F_p*, n=|H|, and

    r_s(x) = #{(h_1,...,h_s) in H^s : sum h_i=x},
    E_s = sum_x r_s(x)^2,
    M = max_{a != 0} |sum_{h in H} exp(2 pi i a h/p)|.

For n>=4 and n^4/4<=p<=n^4, the argument below gives

    M << n^(71/72) (log n)^(1/18).                         (A)

The constant is absolute and uniform in this window. The argument does
not need n to be dyadic. This improves the *baseline recorded in this
project*, n^(2849/2880+epsilon), by 1/320 in the power of n, after allowing
an arbitrarily small positive epsilon. It is far from the requested
square-root bound, n^(1/2+o(1)). No general two-set Paley estimate or
transfer to the prize parameters follows.

## Literature inputs and a direct check

The three inputs use the same sum-equality energy E_s:

1. [MRSS, Corollary 7](https://arxiv.org/html/1712.00410v1):
   E_3 << n^4 log n for n=O(sqrt(p)).
2. [Shkredov, Lemma 9, equation (20)](https://arxiv.org/html/1705.09703v1):
   for nonempty H-invariant Q contained in F_p*,

       E_2(Q) <= C_* (|Q|^4/p + |Q|^3/sqrt(n)), C_*>=1.    (B)

   The [published paper](https://msp.org/moscow/2019/8-1/moscow-v8-n1-p03-s.pdf)
   has the same statement as equation (23), on printed page 20.

3. [Konyagin's inequality, Shkredov Lemma 7 (15)](https://arxiv.org/abs/1311.5726):

       M^(2kl) <= p E_k E_l n^(2kl-2k-2l).                 (C)

The [earlier project derivation](sigma-subgroup-2026-09-05.md) of (C)
also applies to odd k or l. Its Fourier inverse of |eta|^k can be signed;
only its squared norm is used. Complex Parseval uses absolute squares.

As a check that the power improvement does not depend on the refinement
below, Shkredov's equation (25), for s=3, already supplies

    E_6 << (log n)^4 (n^12/p + n^(11/2) E_3).

Combining this with (C) at (k,l)=(3,6) gives
M << n^(71/72)(log n)^(1/6) in the quartic window. The recurrence is
stated for every integer s>=2 inside the proof of Theorem 12. It is
equation (28) on printed pages 22-24 of the published version. It is not
restricted to dyadic s, and its proof precedes the size condition used
for iteration. Neither a centered-function theorem nor an omission of
the origin is needed here.

The following elementary lifting of (B) removes the logarithmic factor
from that recurrence. It improves the logarithm in (A), not its power.
The PDF checks and exact source hashes are recorded in the
[source ledger](../results/parallel28_source_scope_2026_09_05.json).

## Lifting an invariant-set estimate to nonnegative weights

For a nonnegative function g supported outside zero, define

    E(g) = sum_z (sum_x g(x) g(z-x))^2.

With the unnormalized additive Fourier transform,
E(g)=p^(-1) sum_a |g_hat(a)|^4. Thus E(g)^(1/4) is a norm and satisfies
the triangle inequality. Suppose g is H-invariant, and put
A=sum g and B=sum g^2. The zero function is immediate; otherwise A,B>0.
Its superlevel sets P_t={x:g(x)>t} are H-invariant subsets of F_p*.
The identity g=integral_0^infinity 1_(P_t) dt is a finite step-function
decomposition, so the norm triangle inequality and (B) give

    E(g)^(1/4)
      <= C_*^(1/4) [p^(-1/4) A + n^(-1/8) I],
    I = integral_0^infinity |P_t|^(3/4) dt.

Here |P_t|<=min(A/t,B/t^2). Split the integral at u=B/A:

    I <= integral_0^u (A/t)^(3/4) dt
           + integral_u^infinity (B/t^2)^(3/4) dt
      = 4 A^(3/4) u^(1/4) + 2 B^(3/4) u^(-1/2)
      = 6 A^(1/2) B^(1/4).

Consequently (x+y)^4<=8(x^4+y^4) yields

    E(g) <= 8 C_* A^4/p + 10368 C_* A^2 B/sqrt(n).         (D)

This statement explicitly excludes a mass at zero. It makes no assertion
about the unqualified signed-function theorem discussed in the
[earlier source audit](analytic-bounds-and-amplification.md).

## The origin and a recurrence uniform in the moment order

Apply (D) to g=r_s-r_s(0)1_{0}, for any integer s>=1. The invariance
r_s(hx)=r_s(x) follows by multiplying each summand by h. Write z=r_s(0).
We have

    sum g <= n^s,  sum g^2 <= E_s,
    z <= n^(s-1),  z^2 <= E_s,
    E(r_s)=E_(2s), E(z1_0)=z^4.

A second norm triangle inequality gives

    E_(2s) <= 8 E(g) + 8 z^4
      <= 64 C_* n^(4s)/p
         + 82944 C_* n^(2s-1/2) E_s
         + 8 n^(2s-2) E_s
      <= 82952 C_* [n^(4s)/p+n^(2s-1/2) E_s].            (R)

The last step uses n>=1 and C_*>=1. In particular the contribution at
zero is retained and absorbed with an explicit inequality. The constant
in this ordinary proof is independent of s; no assertion that it is
optimal is made. This is a consequence of (B), not a new incidence input.

## Substitution and scope of the improvement

Take L=log n>=1 and E_3<=A_0 n^4 L, increasing A_0 to be at least one.
With C_0=82952 C_*, (R) at s=3 and (C) at (3,6) imply

    M^36 <= p n^18 E_3 E_6
      <= C_0 A_0 n^34 L + C_0 A_0^2 p n^(63/2) L^2.

Since p<=n^4, n>=1 and L>=1, each term is bounded by its constant
times n^(71/2)L^2. Taking the 36th root proves (A).
The lower end p>=n^4/4 is used to guarantee n<sqrt(p), so that the MRSS
input applies. A prime cannot equal the endpoints that are perfect
powers. In fact n>=4 is more than sufficient for n^4/4>n^2.

The exact comparison is

    2849/2880 - 71/72 = 1/320,
    1 - 71/72 = 1/72.

The proof uses no conjectural D6 bound. Conversely it does not provide
one: the six-variable energy remains bounded by n^4 log n, and (R)
controls the twelve-variable energy E_6. Confusing E_6 with the count
of six-term relations would invalidate the claim.

The preceding literature survey did not combine the arbitrary-s
recurrence with E_3 and the mixed (3,6) Konyagin inequality. Its earlier
"no new exponent" assessments remain historical reports; this note
supersedes that conclusion for this specific subgroup baseline.
It does not establish literature novelty, current-best status, a useful
numerical constant at the finite examples, or square-root cancellation.
No source author or independent reviewer has been contacted.

## What iteration of these inputs supplies

Ignoring logarithms for this exponent comparison only, set
d_s=2s-e_s when E_s is bounded by n^e_s. The seed inputs give
d_1=1, d_2=31/20 and d_3=2. Recurrence (R) increases the deficit
by 1/2 under doubling, capped at 4 by the principal term in this window.
Log-convexity of the Fourier moments allows linear interpolation of the
deficits. The resulting concave, piecewise linear envelope has knots

    (1,1), (2,31/20), (3,2), (6,5/2),
    (12,3), (24,7/2), (48,4),

and is constant at 4 thereafter. Doubling preserves this envelope: on
[3,48] it increases by 1/2 until capped, and direct comparison on
[1,3] gives an increment of at least 1/2. The power-of-two chain from
the second-moment seed lies below this envelope. Thus further doubling
and interpolation of these seeds do not improve it.

At two orders k,l the saving supplied by (C) is
(d_k+d_l-4)/(2kl). On each pair of linear pieces this is bilinear in
1/k and 1/l, so its maximum occurs at endpoints. For k or l above 48,
a positive saving can only decrease as that variable grows. Checking
the 49 pairs of knots therefore exhausts this specific envelope. The
maximum is 1/72, at (3,6), (3,12), (6,3), (6,6), and (12,3).

This explains why merely running the recurrence to greater depths does
not approach exponent 1/2. It is a limit of the estimates supplied by
these inputs and this use of (C), not a lower bound on M or an obstruction
to other arguments. The [rational ledger](../results/parallel28_verification_2026_09_05.json)
checks the endpoint arithmetic. Its finite checks do not constitute a
formal verification of the analytic proof.
