# Triple fibers: exact certificates and an aggregate energy criterion

The full Paley graph conjecture and Proximity Prize remain unproved.
This pass rules out a proposed maximum-fiber-two condition inside the
quartic prime range. It completely classifies that condition at fixed
dyadic orders through 128, and isolates an arithmetic quantity sufficient
for the intermediate energy and triangle bounds. No uniform estimate
for that quantity is proved.

## 1. Reuse the existing polynomial

The integer-polynomial and discriminant reductions already appear in
[the kernel note](kernel-discriminant.md) and
[the primitive-fiber note](parallel2-subgroup-2026-09-04.md).
The power sums in pass 47 recover the same underlying root data; they
are not a new integer-polynomial construction.

For n=2k>=4 dyadic, d=n/4, and a primitive complex n-th root zeta, set

    R_n(Y)=product_(1<=j<n/2, j odd) [Y-(1+zeta^j)^n].

This is a monic irreducible polynomial in Z[Y] of degree d, with distinct
negative real roots. Galois transitivity and their distinct absolute
values prove irreducibility, as in the earlier note. For a prime
p=1 mod n, it splits over F_p. If c_b are its positive root
multiplicities, and K=mu_k, L=mu_n\K, the existing exact identity is

    B_n=||1_K*1_L||_2^2=n sum_b c_b^2,   sum_b c_b=d.       (1)

Each c_b is an actual additive representation count on one mu_n coset.
In particular c_b<=2 for every b would imply B_n<=2k^2. This sufficient
pointwise condition was explicitly unproved in the earlier note.

## 2. A finite list containing every triple-fiber prime

For n>=8, define the positive integer

    G_n=gcd(|Res(R_n,R_n')|, |Res(R_n,R_n'')|).             (2)

Both resultants are nonzero. Indeed R_n is irreducible over Q, and its
first and second derivatives are nonzero polynomials of smaller degree.
Set G_4=1, since R_4 has degree one and cannot have a triple root.

For p=1 mod n one has p>n>d. If R_n modulo p has a root of multiplicity
at least three, that root is common to R_n, R_n', and R_n''. Therefore

    some c_b>=3  implies  p divides G_n.                  (3)

The converse is not asserted in general: separate resultants can vanish
at different roots. A complete factorization of G_n, followed by a
modular gcd or direct root check at each eligible prime factor, is a
complete fixed-order classification. It does not depend on a search
cutoff for p.

There is also a general way to remove this possible overcount without
factoring. Put Z_n(t)=Res(R_n,R_n'+t R_n''), a polynomial in t of degree
at most d. Modulo p, it is the product over the d roots lambda_i of
the linear expressions R_n'(lambda_i)+t R_n''(lambda_i). Thus it is
identically zero exactly when some root is common to all three
polynomials. Since p>d, evaluations at t=0,...,d detect whether Z_n is
zero. Consequently the primes p=1 mod n dividing

    gcd_(t=0,...,d) |Res(R_n,R_n'+t R_n'')|

are exactly the triple-fiber primes. This observation is an exact
criterion, not a height bound. The finite computation below uses the
simpler (2), with all of its candidate primes checked individually.

## 3. Count only the fibers beyond two

Define

    Y_n(p)=sum_(c_b>=3) c_b(c_b-2),
    s_1=#{b:c_b=1}.

Equation (1) gives the exact decomposition

    B_n=2k^2+n[Y_n(p)-s_1] <= 2k^2+nY_n(p).              (4)

Thus double fibers cost nothing above the baseline 2k^2. The arithmetic
certificate controls this positive contribution:

    Y_n(p)<=v_p(G_n).                                    (5)

To prove (5), lift the primitive n-th roots uniquely to Z_p; the roots
lambda_i of R_n then lie in Z_p and remain distinct in characteristic
zero. A reduction class of size c contributes at least c(c-1) to
the valuation of Res(R_n,R_n'), which is the discriminant up to sign.
At any root lambda_i,

    R_n''(lambda_i)=2 sum_(j!=i)
                       product_(ell!=i,j)(lambda_i-lambda_ell).

If its class has size c>=3, every term contains at least c-2 factors
divisible by p. The c members therefore contribute at least c(c-2)
to v_p(Res(R_n,R_n'')). Other roots have nonnegative contributions.
Both resultant valuations are at least Y_n(p), which proves (5).
Cancellation in the displayed sum can only increase the valuation.
Equality in (5) is not assumed.

## 4. A sufficient condition for the energy target

Fix a dyadic N>=4 and p=1 mod N. All smaller dyadic groups are in the
same prime field. Write E_s=E(mu_s). The existing exact tower identity
and Cauchy give

    E_s=2E_(s/2)+6B_s+8T_s,
    T_s^2<=E_(s/2) B_s,
    E_s<=3E_(s/2)+22B_s.                                 (6)

The last inequality follows from 8 sqrt(E B)<=E+16B. Combining (4)
with (6), and dividing by s^2, gives

    E_s/s^2 <= (3/4) E_(s/2)/(s/2)^2 + 11 + 22Y_s(p)/s.

Since E_2=6, iteration proves

    E_N <= 44N^2 + 22N^2 V_N(p),                         (7)
    V_N(p)=sum_(s=4,8,...,N)
             (3/4)^(log_2(N/s)) Y_s(p)/s.

The same upper bound holds with Y_s(p) replaced by v_p(G_s).
In the quartic prime range, a uniform V_N(p)<<N^(1/3) would imply
E_N<<N^(7/3), and then the already scoped single-coset bound in
[pass 43](parallel43-energy-prime-exceptions-2026-09-06.md) gives the
intermediate triangle power17/3. For example, uniformly
v_p(G_s)<<s^(4/3) at every smaller level would suffice by geometric
summation. Neither this valuation condition nor the weaker condition
on V_N is proved. The triangle target itself is not the full Paley or
prize goal.

The basic height estimate remains inadequate. The roots of R_n lie in
(-2^n,0), so |Res(R_n,R_n')|<=2^(n d(d-1)). It gives only
v_p(G_n)<=n d(d-1)log(2)/log(p), of order n^3/log(n) at the largest
quartic level. This is worse than even the trivial O(n^2) bound on Y_n.
The new decomposition does not turn that height estimate into a saving.

## 5. Complete fixed-order classification

The checker recomputes both integer resultants, verifies their gcd's
complete prime factorization, and certifies every factor by trial
division through its square root. All odd factors in these five cases
are splitting primes and have a positive triple gcd; there are no
false candidates among these particular factorizations.

| n | Number of triple-fiber primes | Largest such prime | Triple-fiber primes in [n^4/4,n^4] |
|---:|---:|---:|:---|
| 8 | 0 | none | none |
| 16 | 1 | 17 | none |
| 32 | 5 | 2,113 | none |
| 64 | 17 | 697,601 | none |
| 128 | 47 | 973,021,549,697 | 77,796,353; 118,593,281; 181,312,129 |

The order-128 interval is [67,108,864,268,435,456]. At each of its three
exceptional primes, the 32 primitive roots have multiplicity histogram

    27 singletons, one doubleton, one tripleton.

Thus Y_128=3, s_1=27, and

    B_128=128(27+4+9)=5120 <8192=2*64^2.                 (8)

At all other eligible primes in that interval, every fiber is at most
two. Hence B_128<=8192 throughout the entire interval, despite the
failure of the maximum-fiber-two condition at exactly three primes.
Every smaller dyadic level also has fibers at most two there, since its
largest possible triple prime is at most697601. This is a finite tower
conclusion, not a statement uniform in growing N.

One explicit triple is at p=77796353 with generator g=7906121 of order
128. For K=<g^2>, L=gK, the number x=725440 has the three representations

    (65398459,13123334), (49300312,29221481),
    (75076016,3445777)

in K x L, with sums interpreted modulo p. Direct enumeration verifies
there are exactly three. This is an actual quartic-field obstruction to
the proposed pointwise condition, not a counterexample to B<=2k^2,
the needed energy bound, or the Paley conjecture.

## 6. Reproduction and review scope

The [checker](../experiments/parallel48_collision_eliminants.py) uses
the existing integer Newton construction and python-flint0.9.0 for
exact polynomial resultants. A second algorithm computes determinants
of the integer multiplication matrices and agrees in all five orders.
The standard-library Bareiss implementation also agrees through order64.
Each candidate factor receives a separate trial-division certificate;
factorization-library output alone is not treated as a primality proof.
All70 splitting-prime cases have independent literal pair counts and
modular polynomial checks. The detailed
[certificate](../results/parallel48_collision_eliminants_2026_09_06.json)
includes the complete factorizations, integers, root hashes and explicit
quartic witnesses.

Dependencies are pinned in
[requirements](../experiments/requirements-exact-arithmetic.txt).
The observed local invocation is:

    /opt/miniconda3/bin/python3 experiments/parallel48_collision_eliminants.py --flint-path tmp/arithmetic/python-flint-0.9.0-cp313

That isolated backend was installed with uv using the existing Conda
Python3.13. Its initial system-temporary directory disappeared between
continuations, so the pinned backend was restored under the project's
already ignored tmp/arithmetic directory. The system Homebrew Python's installer had a dynamic-library
error; its installation was left unchanged. The earlier slow Python
resultant calculation finished, and a separate exploratory order256
factorization was stopped after preserving the complete order128 data.
No order256 certificate is claimed. No Lean job was launched or changed.
The [source and tool scope](../results/parallel48_source_scope_2026_09_06.json)
separates reused mathematical arguments from arithmetic-library APIs.

Root completed the ordinary derivations and exact checks. No separate
agent, Lean, external review or novelty claim is made. The next task is
to control the aggregate V_N or another sufficient uniform input, rather
than rely on the now refuted maximum-fiber-two hypothesis. All existing
uniform exponents and the full goals remain unchanged.
