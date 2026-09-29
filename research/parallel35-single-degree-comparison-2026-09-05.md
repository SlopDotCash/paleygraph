# A single-degree criterion for the centered opposite-free count

**The Paley conjecture and the Proximity Prize remain unproved.**
The positive upper estimate is still missing. This pass proves an
unconditional Gaussian-scale lower estimate and replaces the pass34
hierarchy criterion with a comparison at one degree. It uses an
established theorem on derivatives of real-rooted polynomials, applied
to each Fourier row before averaging. The resulting argument is an
ordinary proof with exact finite checks; separate-author review and
Lean verification remain outstanding. No novelty is claimed.

## 1. Statement

Let p be an odd prime, A=-A subset F_p^*, n=|A|>0, and N=n/2.
For every integer s>=1, define

    eta(a) = sum_(x in A) exp(2 pi i a x/p),
    T_s = (1/p) sum_(a != 0) |eta(a)|^(2s),
    O_(2s) = number of ordered distinct zero-sum 2s-words in A
             with no pair of opposite entries,
    B_s = O_(2s) - 2^(2s)(N)_(2s)/p.

Falling factorials are zero when 2s>N. Then

    B_s >= -(256sn)^s,                                    (L)
    |B_s| <= 4^s [T_s+(64sn)^s],                          (A)
    T_s <= 16^s B_s+2(4096sn)^s.                          (R)

These statements require neither multiplicative closure nor a relation
between p, n, and s. In particular they do not need the small-epsilon
condition from pass33 or estimates at smaller degrees from pass34.
The constants are deliberately conservative.

If an upper bound B_s<=(Ksn)^s is proved at one degree, then (R) gives

    T_s <= [(16K+8192)sn]^s.                              (C)

Conversely, T_s<=(Ksn)^s implies |B_s|<=[4(K+64)sn]^s by (A).
Thus the two upper-bound problems are equivalent up to absolute
constant changes at the same degree. The negative side is controlled
unconditionally by (L); the positive side remains the obstruction.

## 2. The external input and its exact hypotheses

[Ravichandran, arXiv:1609.04187v2, Theorem 4.4](https://arxiv.org/html/1609.04187v2#S4.Thmtheorem4)
states the following. If a real-rooted polynomial of degree d has all
roots in [-1,1], whose sum is zero, and k/d>=1/2, then the roots of
its k-th derivative lie in

    [-2 sqrt((k/d)(1-k/d)), 2 sqrt((k/d)(1-k/d))].         (D)

We use an integer k with k<d. The statement and its hypotheses were
checked in the primary HTML and visually in the PDF on printed page12.
The [source record](../results/parallel35_source_scope_2026_09_05.json)
pins both archives. No theorem concerning averaged polynomials is
being imported.

## 3. A pointwise comparison for bounded real rows

Let y_1,...,y_N be arbitrary real numbers in [-2,2]. Put

    F = sum_i y_i,       n=2N,       r=2s,
    H = r! e_r(y_1,...,y_N),

where e_r is the elementary symmetric polynomial, zero for r>N.
We prove the pointwise forms of (L), (A), and (R), replacing B_s by H
and T_s by F^(2s).

First suppose r<=N/2. Write mu=F/N and z_i=y_i-mu, so
sum z_i=0 and |z_i|<=4. Consider

    P(x)=product_i (x+z_i),
    e_r(x+z_1,...,x+z_N)=P^(N-r)(x)/(N-r)!.

The second identity follows either by differentiation or by choosing
which factors contribute z_i. Apply (D) to the polynomial with roots
-z_i/4, taking d=N and k=N-r. It follows that the roots lambda_j
of this degree-r derivative satisfy

    |lambda_j| <= 8 sqrt((r/N)(1-r/N)).

After evaluating at x=mu and multiplying by r!, we obtain

    H = a product_(j=1)^r (F-theta_j),
    a = (N)_r/N^r,       theta_j=N lambda_j,
    2^(-r) <= a <= 1,
    |theta_j| <= 8 sqrt(r(N-r)) <= R,
    R^2=64rN=64sn.                                      (P)

The bound on a uses j/N<1/2 in each factor (1-j/N).

Since r is even, H>=0 whenever |F|>=R. When |F|<=R each factor
has magnitude at most 2R. Consequently

    H >= -(2R)^r = -(256sn)^s.                            (P1)

For all F, convexity gives

    |H| <= (|F|+R)^r
         <= 2^(r-1)(F^r+R^r)
         <= 4^s [F^(2s)+(64sn)^s].                       (P2)

If |F|>=2R, every factor has magnitude at least |F|/2 and their
product is nonnegative. Thus H>=a F^r/2^r>=F^r/4^r, or
F^r<=16^s H. If |F|<2R, let L=(2R)^r=(256sn)^s. Then F^r<=L
and H>=-L, so

    F^r <= 16^s H+(1+16^s)L
         <= 16^s H+2(4096sn)^s.                          (P3)

Now suppose r>N/2. Then n<8s, and the elementary estimate
|H|<=r! binomial(N,r)2^r<=n^r applies when r<=N; for r>N one has
H=0. Also |F|<=n. Hence both |H| and F^r are at most
(8sn)^s<=(256sn)^s. The lower and absolute bounds above follow
immediately, and the same argument using H>=-L and F^r<=L proves
(P3). This covers all N>=1 and s>=1, including r>N.

This extends the project's [pointwise sign-row approach](parallel21-squarefree-moments-2026-09-05.md)
to bounded real entries. Here (D) supplies the root control; an exact
three-term Krawtchouk recurrence is not assumed for the cosine entries.

## 4. Fourier averaging preserves the principal subtraction

Choose one representative h_i of every opposite pair in A and set

    y_i(a)=exp(2 pi i a h_i/p)+exp(-2 pi i a h_i/p).

These values are real and belong to [-2,2], and sum_i y_i(a)=eta(a).
The coefficient e_r(y(a)) selects r different opposite classes and a
sign in each. Additive character orthogonality therefore gives

    (r!/p) sum_a e_r(y(a)) = O_r.

At a=0, every entry is 2, so r!e_r(y(0))=2^r(N)_r. It follows
exactly that

    B_s = (1/p) sum_(a != 0) H(a).                       (F)

Average (P1), (P2), and (P3) over the nonzero frequencies, with mass
1/p at each. The total mass is (p-1)/p<=1, giving (L), (A), and (R).
No coefficientwise absolute value is taken before removing a=0.

For K>=0, the elementary inequalities

    (16K)^s+2*4096^s <= (16K+8192)^s,
    4^s(K^s+64^s) <= [4(K+64)]^s

prove the single-degree implications in Section1.

For a multiplicative subgroup, nonzero Fourier values repeat on
H-cosets. Thus M^(2s)<=(p/n)T_s. If the still-unproved hypothesis
B_s<=(Ksn)^s holds, then

    M <= (p/n)^(1/(2s)) sqrt((16K+8192)sn).

Taking s>=log(p/n) would give M<=sqrt(e(16K+8192)sn). In the
quartic window this uses s=O(log n). This is a conditional subgroup
criterion, not a proof of the required upper hypothesis or a complete
reduction to the classical arbitrary-set Paley conjecture or prize.

## 5. A stronger lower estimate at degree ten

Equation (L) gives, for all eligible fields and symmetric sets,

    O_10 >= 2^10(N)_10/p - (1280n)^5.                    (T10)

For n>=90,

    2^10(N)_10 = product_(j=0)^9 (n-2j)
               >= n^10(1-90/n).

The product bound uses product(1-u_j)>=1-sum u_j for u_j in [0,1].
When p<=n^4, (T10) implies

    O_10 >= n^6-(90+1280^5)n^5,
    90+1280^5 = 3435973836800090.                         (T11)

Thus the relative possible shortfall below the principal term is
O(1/n) at degree ten, improving the O(n^(-1/5)) lower estimate from
pass34 asymptotically. The explicit constant is much larger, so the
earlier n>=2^35 sufficient threshold for O_10>=n^6/2 remains better.
This is a lower estimate only. O_10 can exceed its principal term by
an amount this argument does not control. No existence of a quartic
prime/subgroup pair at each large order is asserted.

## 6. Why averaging real-rooted polynomials is not a shortcut

Each row polynomial product_i(1+w y_i(a)) is real-rooted. Its
average over a!=0 need not be. An actual quartic subgroup example is

    p=1153, n=8,
    H={1,75,123,140,1013,1030,1078,1152}.

Its four opposite classes have representatives 1,75,123,140.
Enumerating their 3^4 choices (absent, plus, minus) shows that the
empty selection is the only opposite-free zero-sum selection. Hence

    sum_(a != 0) product_i(1+w y_i(a)) = 1153-(1+2w)^4.

The substitution t=1+2w gives 1153-t^4, with two real and two
nonreal roots. This refutes an averaged-real-root assertion even for
an eligible subgroup. It does not refute a Paley or moment conjecture.
The proof above applies the root theorem to each shifted row separately.

The direct coefficient recurrence also retains an unestimated term.
Writing E_r(a)=e_r(y(a)) and P_j(a)=sum_i y_i(a)^j, Newton's identity is

    r E_r(a)=sum_(j=1)^r (-1)^(j-1) P_j(a) E_(r-j)(a).

For the paired Fourier row,

    P_j(a) = 1_(j even) (n/2) binomial(j,j/2)
             + sum_(0<=l<j/2) binomial(j,l) eta((j-2l)a).

Frequencies are interpreted modulo p. After averaging, the j=1 term
is the correlation of eta(a) with E_(r-1)(a), still of total degree r.
The scalar means of lower coefficients do not estimate that correlation
in this calculation. The new lower estimate and one-degree comparison
do not turn this recurrence into the missing positive upper bound.

## 7. Verification and scope

The [verifier](../experiments/parallel35_verify_2026_09_05.py) passes
1,053 exact checks. It tests 114 bounded rational rows, including both
size regimes, asymmetric rows, and large-sum cases. For 78 rows it
uses exact Sturm arithmetic at quadratic-radical endpoints to count
all distinct roots inside the claimed radius. Repeated roots are
handled through the terminal gcd in the Sturm chain. It also checks
the normalization and shifted polynomial identity, the inequalities
on 43 previously counted finite-field cases, and all81 selections
in the averaged-polynomial counterexample.
The paired power-sum formula and Newton recurrence are also compared
as full group-algebra coefficient maps in four prime-field cases.

The [record](../results/parallel35_verification_2026_09_05.json)
uses only rational/integer arithmetic and takes about1.71seconds in
the recorded run. These checks validate finite implementations, not
the universal derivative-root theorem or the full conjecture. No Lean
build or hosted proof submission was launched. The positive upper
estimate on B_s remains unproved, so the period exponent71/72 and
the full requested targets are unchanged.
