# Independent review of the global third-moment compiler

2026-09-05. The accelerated computation agrees exactly with every independent complete-subset oracle reviewed here. No numerical or mathematical defect was found in the third-moment algorithm. This review does not claim historical originality or a bound for the prize problems.

The reviewed source hashes are:

- `third_moment.py`: `f21114388e4a5b4e2b929b3db14ddd2987ecf22e70fc96051944422044bf207c`
- `union_coefficients.cpp`: `9f95bef6f2f4866de7014ae9d16bfeadfacc7db1a12e13ffb6b602783c4d0de6`

## Independent evidence

[direct_third_oracle.py](direct_third_oracle.py) constructs the prime-field matrices and both GF49 matrices independently, then directly evaluates row elementary coefficients on each subset. It does not import the accelerated compiler. Its [small output](direct_third_oracle.json) includes all Paley13/17 subsets of sizes6,7,8 and all55 combinations of subset size and degree on Paley9. It also checks all eight bulk sign-type formulas against literal matrices for37,934 normalized or unordered triples.

The GF49 census uses [origin_third_census.cpp](origin_third_census.cpp), with two49-bit row masks and a small exact elementary-coefficient table. Each origin-containing subset is visited once. Translation invariance gives

```
n * sum_(all n-sets C) f(C) = 49 * sum_(n-sets C containing0) f(C)
```

for each power of T6. This double count does not require translation orbits to be free, which matters at n7. Each translation orbit meets the origin-containing slice, so minima and maxima are also global. The code bounds the full worst-case sum of absolute cubes by `1164709865022979968 < 2^63` for its supported n≤8. It rejects larger sizes instead of assuming the integer budget extends.

| Matrix | n | Exact E[T6³] | Minimum | Maximum |
|---|---:|---:|---:|---:|
| Paley49 |6|44555/35673|−29|27|
| Peisert49 |6|6145/3243|−29|27|
| Paley49 |7|−99128449/1533939|−59|77|
| Peisert49 |7|−1672321/1533939|−59|61|
| Paley49 |8|−4737196096/1533939|−132|116|
| Peisert49 |8|−4199178304/1533939|−124|132|

The saved outputs are [n6](twins_n6_third_oracle.json), [n7](twins_n7_third_oracle.json), and [n8](twins_n8_third_oracle.json). Their origin-slice counts are1,712,304;12,271,512;73,629,072. These certify the full13,983,816;85,900,584;450,978,066 subsets. Observed C++ runtimes were0.59s,3.23s,24.83s respectively; compilation and concurrent work add overhead recorded separately.

[review_union_backend.py](review_union_backend.py) expands one literal column at a time, using the seven possible nonempty row selections. It shares neither Newton identities nor the grouped power-sum/cover-number calculation with the backend. All156 arbitrary ternary-column cases, degrees0–6, agree in681 exact union coefficients. [Evidence](review_union_backend.json).

[review_global_third.py](review_global_third.py) compares both direct-inventory and interpolated compiler modes against40 supported exhaustive cases, yielding80 exact third-moment agreements. Its25 odd-degree cases have zero direct third moment and are correctly rejected by the even-degree normalized interface. It additionally checks every ordered distinct row triple against the normalized joint inventory, a total of227,388 triples. This checks aggregate normalization independently of the previously verified explicit Peisert semilinear maps. All64 interpolated classes agree,54 extra-node checks pass, and24 small cases exercise the direct-inventory fallback. The synthetic raw-coefficient check finds zero seventh finite difference for all36 tested union coefficients in the four sign classes. [Evidence](review_global_third.json).

## Mathematical audit

The cover coefficient in the C++ backend is

```
sum_(h=0)^u (-1)^(u-h) binom(u,h) binom(h,a) binom(h,b) binom(h,c).
```

It counts three subsets of a u-element set whose respective sizes are a,b,c and whose union is the whole set. This is exactly the coefficient of `t1^a t2^b t3^c` in `((1+t1)(1+t2)(1+t3)-1)^u`. Multiplying by the appropriate sign power sum and applying the ordinary Newton recurrence therefore computes the elementary polynomials in the column factors. The final replacement `z^k -> (n)_k/(q)_k` is the probability that a fixed set of k distinct columns lies in a uniform n-set; there is no extra binomial normalization.

The degree bound is structural: in the formal logarithm, a tau-dependent term must have odd positive exponent in each of the three row variables. At most d such factors can reach target degree(d,d,d). Hence every extracted union coefficient, and its inclusion-weighted sum, has tau degree at most d at fixed q and mutual signs. The code interpolates from d+1 admissible integral nonnegative sign inventories and checks another when available. These are synthetic inventories; no assertion of graph realizability is needed. The extra-node and finite-difference tests are implementation checks, not substitutes for the degree proof. Small fields use their actual inventory when insufficient nodes exist.

The repeated-row terms retain complete type inventories, as required when n>d. Their two mutual-sign branches agree for even degree by common sign complementation. The normalized distinct-row count is q(q−1) per parameter, with no additional factor6. Odd degree requires signed normalization accounting and is correctly outside the current normalized API.

The new inexpensive inventory checks also follow directly from the conference identities. Let `v_x=S_0x S_1x`, with `S_01=1`, and `tau=S v`. Then `sum v=-1`, `||v||²=q−2`, and `tau_0=tau_1=−1`. The identities `S1=0` and `S²=qI−J` give, after removing those two singular parameters,

```
sum tau=2;
sum tau²=(q−3)(q+1);
sum S_0t*tau_t = sum S_1t*tau_t =2.
```

The normalized pair has(q−5)/4 bulk++ entries and(q−1)/4 entries in each other sign class. These conditions are necessary diagnostics; they do not certify a supplied abstract inventory's realization or its normalization symmetries. The API states that contract explicitly. Optional `edge_power_sums` metadata now requires exactly four unique sign classes with matching values. The final replay independently rejects missing, duplicate, empty and incorrect metadata; all80 numerical comparisons still pass at the final source hash.

All polynomial coefficients use arbitrary-precision integers. The small fixed cover-number calculations have u≤18 and target degrees≤6. Histogram power sums have magnitude at most the validated population≤10^12; all remaining large products enter `cpp_int`. No fixed-width overflow issue was found within those declared bounds.

## What the extra information does and does not establish

The compiler now measures a distributional statistic that separates actual twins whose full Johnson L2 spectra coincide. This is a demonstrated additional capability. It does not determine the distribution or the extrema. The new complete n7/n8 data gives a particularly concrete limitation: Peisert's third moment is larger at both sizes, while its maximum is smaller at n7 and larger at n8. Third-moment ordering alone therefore does not order even these two finite upper extrema.

The arithmetic count backend, exact type compiler, proved input budget, and independent acceptance suite form a new local instrument. Classical finite-population moments, Krawtchouk coefficients, Newton identities, interpolation, Legendre traces and their Hecke moment formulas are established ingredients. The current evidence supports this specific implementation and specialization, not a claim that no one has considered the combination.

## Focused prior-art follow-up

[Feinsilver–Schott, *Krawtchouk transforms and convolutions*](https://link.springer.com/article/10.1007/s13373-018-0132-2), published22 October2018, §2 explicitly represents Krawtchouk coefficients as elementary symmetric functions of signs; §3 derives product linearization and convolution identities through generating functions. This is direct prior art for the sign-coefficient/product machinery, although it does not establish the present conference inventory compiler or the fixed-tau input budget. [Brouwer–Martin, *Triple intersection numbers for the Paley graphs*](https://aeb.win.tue.nl/preprints/p3g.pdf),2021 preprint/2022 publication, already identifies Paley triple intersection data with cubic character sums and elliptic point counts. Those identifications are not new here. The broader Legendre/Hecke prior art and normalization qualifications remain in the frozen [preflight ledger](../preflight/README.md).

A focused local scan again found the same close predecessors listed in preflight. In addition, [parallel28-moment-recurrence](../../../research/parallel28-moment-recurrence-2026-09-05.md) already combines an additive-subgroup third energy with invariant-set estimates, but its random object is additive representation multiplicity, not T6 over uniform column subsets. This distinction avoids both overlooking local moment work and falsely identifying two different moments. No exact predecessor of this generic compiler was found in the bounded scan; absence is not established.

Queries and reads on2026-09-05: `Paley graph triple intersection numbers Legendre elliptic curves moments Brouwer Martin 2021`; `Krawtchouk polynomial finite population U statistics third moment elementary symmetric Newton identities`; full primary Feinsilver–Schott HTML, Brouwer–Martin author PDF; local `rg` for third moment, tau-degree and Newton in research, experiments and rounds2–4, followed by the parallel28 note. A Griffiths2016 primary-paper search result also described the elementary-symmetric representation, but the direct page returned HTTP429 and is not used as an additional verified source.

Replay the three scripts with `python3 tooling_lab/round5/review/direct_third_oracle.py`, its `--twins-n 6|7|8` options, `review_union_backend.py`, and `review_global_third.py`. All writes stay in this review directory except that the compiler's normal build helper may rebuild its own binary when stale; the audited runs used the already-built current binary.
