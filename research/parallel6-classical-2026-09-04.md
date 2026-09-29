# Exact fourth-moment variance on the single-size slice

**Status: an upper bound on a proportion tending to one of the exact
single-size slice is proved. No uniform single-size estimate or Paley
bound is proved.** An exact variance formula also controls the signed
quartic aggregate for a typical set. The literal Gaussian coefficient
3 fails on the same slice in a fixed field. That counterexample does
not refute (SS), whose constant is unspecified.

This pass tried to turn the global variance into an upper bound for
every set. That step did not close: the variance allows exceptional
sets at the scale of the existing termwise Weil bound. The results
below separate the proved upper estimate from that unresolved step.
No novelty or formalization claim is made.

Let p>=11 be prime, chi(0)=0, and choose B uniformly from all n-element
subsets of F_p. Write F_B(x)=sum_(b in B) chi(x-b) and
M_4(B)=sum_x F_B(x)^4. Let (a)_j=a(a-1)...(a-j+1), and put

\[
 v=n(p-n),\qquad f=(n)_3(p-n)_3.
\]

## 1. Exact mean and variance

The mean is

\[
 \boxed{\mu_{p,n}=\mathbb E M_4(B)
       =\frac{v(3v-2p+1)}{p-2}.}                         \tag{1}
\]

The variance depends only on p modulo 4 and n. For p=1 modulo 4 it is

\[
 \boxed{\sigma_{p,n}^2=
 \frac{24f(p-1)
  [(p-5)(p+3)v-3(p^3-6p^2+4p+3)]}
 {(p-2)^2(p-3)(p-4)(p-6)(p-7)}.}                         \tag{2+}
\]

For p=3 modulo 4 it is

\[
 \boxed{\sigma_{p,n}^2=
 \frac{24f(p+1)[(p+1)v-3(p^2-3p+3)]}
 {(p-2)^2(p-4)(p-5)(p-6)}.}                              \tag{2-}
\]

These formulas include n=0,1,2,p-2,p-1,p, where the variance vanishes.
For example, every two-element set has M_4=8p-22. Both formulas are
invariant under n↔p-n, as required by F_(B^c)=-F_B.

Here is an exact derivation, including a reproducible polynomial
certificate for the algebra in (2). For distinct rows, affine change
of variables reduces their joint even moments to x=1,y=0: both row
sums acquire the same factor chi(x-y), which disappears from even
powers. Set epsilon=chi(-1) and define the population power sums

\[
 P_{a,b}=\sum_{z\in\mathbb F_p}\chi(1-z)^a\chi(-z)^b.
\]

When just one exponent is positive, this is zero for an odd exponent
and p-1 for an even exponent. When both are positive, its value is

\[
 \begin{array}{c|cc}
    &b\text{ even}&b\text{ odd}\\ \hline
 a\text{ even}&p-2&-\varepsilon\\
 a\text{ odd}&-1&-1
 \end{array}                                             \tag{3}
\]

The last entry is the two-point character correlation -1. The other
entries retain the excluded zero columns explicitly. Thus (3) uses
only the elementary two-point identity, not an estimate for quartic
curves or any assumption of independent rows.

For a list of positive-degree exponent pairs u_1,...,u_k, let
D(u_1,...,u_k) sum the corresponding products over distinct columns.
With D(empty)=1, exclusion of repeated columns gives the exact recursion

\[
 D(u_1,\ldots,u_k)=P_{u_1}D(u_2,\ldots,u_k)
 -\sum_{j=2}^kD(u_2,\ldots,u_j+u_1,\ldots,u_k).           \tag{4}
\]

For each set partition pi of a labelled list of a x-slots and b
y-slots, record in each block its pair u=(number of x-slots, number
of y-slots). Inclusion probability for its k distinct columns is
(n)_k/(p)_k. Therefore

\[
 J_{a,b}:=\mathbb E[F_B(1)^aF_B(0)^b]
 =\sum_\pi\frac{(n)_{|\pi|}}{(p)_{|\pi|}}
                    D((u_C)_{C\in\pi}).                  \tag{5}
\]

For b=0 the same formula is the one-row moment. Equations (3)-(5)
give

\[
 \mathbb E M_4=pJ_{4,0},\qquad
 \mathbb E M_4^2=pJ_{8,0}+p(p-1)J_{4,4}.                 \tag{6}
\]

Expanding (4)-(6), and subtracting (1)^2, gives (2+)-(2-).
The companion script reconstructs all 4,140 partitions of eight
labelled slots. It verifies the resulting formulas as **identities of
integer polynomials in p and n after clearing denominators**, separately
for epsilon=1 and -1. The residual polynomial is exactly zero. This
algebraic certificate is independent of its subsequent finite-field
tests and is not a fit or interpolation from sampled primes.

## 2. A proved upper estimate on the precise SS slice

Uniformly along n→infinity and n=o(p), equations (1)-(2) imply

\[
 \mu_{p,n}=3pn^2(1+O(1/n+n/p)),\qquad
 \sigma_{p,n}^2=24pn^4(1+O(1/n+n/p)).                     \tag{7}
\]

For the more specific range n=o(sqrt(p)), direct expansion of (1)
gives

\[
 3pn^2-\mu_{p,n}=2pn(1+O(1/p+n^2/p)).                    \tag{8}
\]

Consequently Chebyshev's inequality proves

\[
 \boxed{\frac{\#\{B:|B|=n,\ M_4(B)>3pn^2\}}{\binom pn}
     \le(6+o(1))\frac{n^2}{p}}
 \quad(n\to\infty,\ n=o(\sqrt p)).                      \tag{9}
\]

In particular, for **exactly** n=floor(p^(1/3)), a proportion at least
1-(6+o(1))p^(-1/3) of all n-element sets satisfy M_4(B)<=3pn^2.
This is (SS) at r=2, beta=0, C=3 for that class of sets. The proportion
statement holds for every sufficiently large prime, with the stated
asymptotic as p grows; no averaging over primes is used.

For every such B and every nonempty A, Holder gives simultaneously

\[
 \frac{|S(A,B)|}{|A|n}
       \le\left(\frac{3p}{|A|n^2}\right)^{1/4}.           \tag{10}
\]

At n=floor(p^(1/3)), this gives a power saving when
|A|>=p^(1/3+eta), for fixed eta>0. This conclusion concerns the
specified typical B, with arbitrary A. It is not the uniform
two-set assertion.

## 3. The signed quartic pairing also has a typical upper bound

Retain the exact decomposition M_4(B)=R_B+D_B from
[the earlier cross-ratio proof](parallel2-classical-2026-09-04.md),
where D_B is the centered signed elliptic pairing and

\[
 R_B=p(3n^2-2n)-6n^3+14n^2-9n-6A_B
          +\frac{3(n)_4}{p-2},\qquad
 A_B=\sum_{a\in B}F_B(a)^2.
\]

Conditioning on a∈B leaves a uniform (n-1)-subset of the other p-1
columns, with equally many positive and negative characters. The
elementary second-moment calculation therefore gives

\[
 \mathbb E A_B=\frac{n(n-1)(p-n)}{p-2},\qquad
 \mathbb E R_B=\mu_{p,n}.
\]

It follows exactly that

\[
 D_B=M_4(B)-\mu_{p,n}+6(A_B-\mathbb E A_B),\qquad
 |D_B|\le|M_4(B)-\mu_{p,n}|+6n^3.                        \tag{11}
\]

Thus, for any t>0,

\[
 \boxed{\mathbb P(|D_B|>t+6n^3)\le\sigma_{p,n}^2/t^2.}   \tag{12}
\]

For n=floor(p^(1/3)) and any fixed 0<eta<1/2, taking
t=p^(1/2+eta)n² proves

\[
 |D_B|\le p^{1/2+\eta}n^2+6n^3
\]

outside a proportion at most (24+o(1))p^(-2eta) of sets. For example,
eta=1/6 gives |D_B|<=7p^(4/3) outside O(p^(-1/3)) of this slice.
This is an actual upper bound for the signed pairing on that class;
it is smaller than the uniform sufficient scale pn²≈p^(5/3).
It is still not a bound on its exceptional members.

## 4. A real exception, and the failed variance-to-maximum step

The literal constant 3 in (9) cannot be extended to every member of
the same slice. In the prime field p=1009, take

\[
 B=\{23,44,288,319,336,393,408,485,511,599\}.
\]

Here |B|=10=floor(1009^(1/3)), and direct integer character summation
gives

\[
 \boxed{M_4(B)=316554>302700=3p|B|^2.}                    \tag{13}
\]

The eight rows 1008,997,982,736,916,535,826,312 each equal 10.
The full moment in (13), not just those eight summands, certifies the
strict inequality. This is a fixed counterexample to the coefficient
3, not an asymptotic counterexample to any possible SS constant.

More fundamentally, the exact mean and variance do not supply a
uniform upper bound of the needed scale. For positive mu and variance
sigma², and any T>=sigma²/mu, the two-point random variable taking

\[
 \mu+T\quad\hbox{with probability }\frac{\sigma^2}{\sigma^2+T^2},
 \qquad
 \mu-\frac{\sigma^2}{T}\quad\hbox{otherwise}
\]

is nonnegative and has exactly that mean and variance. On the SS
slice, take mu and sigma² from (1)-(2) and T=sqrt(p)n^4. The upper
value exceeds pn² by a factor of order p^(1/6), while its probability
is of order n^(-4)=p^(-4/3). Its lower value remains positive for
sufficiently large p. Thus positivity and these two statistics are
compatible with an outlier at the old termwise Weil scale.

This synthetic distribution is not asserted to be realized by Paley
sets, or to have the uniform finite sample space of n-subsets. It
pinpoints the missing logical step: mean and variance estimates alone
cannot make the exceptional probability zero. Nor does (9) imply
that every larger adversarial set contains enough good n-subsets to
use the sampling reduction from pass 5.

The typical-set conclusion also should not be mistaken for a new
worst-case cancellation result: elementary random-set concentration
already gives strong cancellation for a random B against every A.
The additional information here is the exact fourth-moment variance,
the Gaussian coefficient on the prescribed slice, and a quantitative
signed-pairing bound. Uniform control of the exceptional sets in
(SS), (CS), and (LM) remains open in this workspace.

## Verification

The [verifier](../experiments/parallel6_classical_2026_09_04.py) checks
the formulas symbolically as described above, then independently
enumerates small fields to verify fixed-size means and variances.
It also checks the centered-remainder identity symbolically and the
explicit exception (13) by literal character sums. No prime scan,
numerical fit, or floating-point arithmetic is used. Counts and
note/script hashes are stored in
[the result](../results/parallel6_classical_2026_09_04.json).
The completed run verifies two symbolic variance identities, the
symbolic mean and centered-remainder identities, and 37 fixed-size
mean/variance classes comprising 24,678 exhaustively enumerated subsets.
