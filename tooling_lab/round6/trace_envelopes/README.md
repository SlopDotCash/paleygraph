# How much arithmetic information does the third moment add?

This tool certifies the range of the global third moment allowed by selected trace constraints. It finds an exact ambiguity after ordinary moments at three tested primes, and also measures how small that ambiguity is. At q1,000,033,n31, even the conference identities and the Hasse-supported residue domain constrain the raw third moment to a range less than6.50×10⁻¹⁰ of its actual magnitude.

These are bounds for a **nonnegative histogram relaxation** of a global average. The feasible rational histograms need not be integer-valued or realized by graphs. The result is neither a bound on individual column sets nor an assertion that richer arithmetic information is useless for other targets.

## Exact model

The frozen round5 compiler supplies two rational polynomials P_m(theta), P_u(theta), of degree at most6, for the monochromatic and mixed triangle classes. Their sixth-degree coefficients coincide and their fifth-degree coefficients vanish. For a prime q≡1 modulo4, retain every theta with

```
theta² <= 4q,
theta mod8 = (15−q) mod8 for monochromatic,
theta mod8 = (q−7) mod8 for mixed,
```

and with integral nonnegative triple-sign counts. Hasse applies to these prime-Paley correlations. It is not used for Peisert inputs; this prototype rejects nonprimes.

Let N_m(theta), N_u(theta) be real nonnegative histogram counts. We use scaled variables `w_m=N_m/q` and `w_u=N_u/(3q)`. The factor3 is essential: the mixed normalized class consists of three equally weighted edge patterns. With R=1+floor(sqrt(q)), the four equality rows are

```
sum w_m = (q−5)/(4q),
sum w_u = (q−1)/(4q),
sum theta*(w_m−w_u)/R = 2/(qR),
sum theta²*(w_m+3w_u)/R² = (q−3)(q+1)/(qR²).
```

The distinct-row objective is

```
q²(q−1) * sum_theta [P_m(theta)*w_m(theta)
                     +3*P_u(theta)*w_u(theta)].
```

Add the exact repeated-row constant to obtain the full third moment. The actual certified inventory is checked to be feasible and to reproduce the compiler's exact value.

The second model also fixes ordinary theta powers1,3,4,5,6 to their actual values. Powers0 and2 are already fixed by the four rows. No class-weighted higher moments are supplied. It therefore tests precisely what can be lost when the two residue classes are merged into ordinary moments through6.

## Numerical discovery, exact proof

HiGHS proposes an optimizer and dual multipliers. Each exported dual is rationalized, then repaired by increasing its class-constant rows until its inequality holds **exactly at every admissible theta**. Dotting those verified inequalities with any feasible nonnegative weights proves the bound. Primal weights are reconstructed by exact rational linear algebra on the proposed support, and all equalities and nonnegativity conditions are checked. No floating-point solver status establishes a mathematical inequality here.

Before numerical optimization the tool subtracts an exactly known linear combination of the constraint rows from the objective. On large inputs the remaining variation is much smaller than the original objective; optimizing the unmodified numbers initially obscured it. This subtraction preserves every feasible objective difference. The exact constraint contribution is added back to the final dual. Saved certificates expose this operation and all repairs.

The output gives both dual outer bounds and feasible primal witnesses, with their exact gaps. It does not label a bound exactly optimal merely because that gap is small. Standardized interval widths are stored through their **exact square**, avoiding the loss caused by subtracting two rounded, nearly equal displayed values.

## Results

The following widths are rounded from exact rational bounds in [results.json](results.json). The last column is the full interval width divided by the absolute actual raw third moment, not a tail probability.

| q | n | Supplied constraints | Raw third-moment interval width | Relative width |
|---:|---:|---|---:|---:|
|101|8|Conference and support|202.539|6.9904×10⁻²|
|1297|8|Conference and support|730.052|2.9968×10⁻¹|
|65537|16|Conference and support|235668.931|1.8746×10⁻⁷|
|1000033|31|Conference and support|21977945.727|6.4909×10⁻¹⁰|
|1297|8|Also ordinary powers through6|0.016596750|6.8129×10⁻⁶|
|65537|16|Also ordinary powers through6|0.009861340|7.8441×10⁻¹⁵|
|1000033|31|Also ordinary powers through6|0.009810834|2.8976×10⁻¹⁹|

The two exact rational feasible histograms at q1297 have identical ordinary powers0–6 and different global third-moment values. Their difference is

```
579595511175839744 / 34922229933188280115.
```

The [independent compact witness](../trace_envelopes_review/rational_ambiguity_witness_1297.json) exposes the actual rational counts. Similar exact relaxed ambiguities occur at65537 and1,000,033. This proves that ordinary moments through6 do not determine the objective on those relaxed feasible sets. It does not establish two actual graphs with those measurements.

At q101 the optimizer's two primal values coincide. Independent review proves a stronger explanation: the ordinary-moment constraints have rank9 on10 support points, and the one remaining affine direction is forced to zero by nonnegativity at two zero-count positions. Thus that particular relaxed histogram is unique. A tiny nonzero dual repair is not evidence of ambiguity.

The finite critical-size bounds justify testing a universal leading-term explanation before investing in more exact global-moment scans. They do not by themselves prove an asymptotic law; that algebraic question is pursued separately.

## Review and reproduction

[verify_results.py](verify_results.py) performs a source-bound readback using the production model helpers, without invoking an optimizer. It checks8 cases,5208 rational dual inequalities,104 primal equalities and rejects16 corrupted certificates.

The [independent reviewer](../trace_envelopes_review/README.md) rebuilds the support, constraints, objective and actual weights without importing either production file or calling an optimizer. It verifies all eight cases, quotient/repair metadata, exact central and standardized values, and56 corruption controls. It also supplies the separate q101 uniqueness and q1297 ambiguity certificates. All source and input hashes are retained.

```
/opt/miniconda3/bin/python3 tooling_lab/round6/trace_envelopes/trace_envelopes.py
/opt/miniconda3/bin/python3 tooling_lab/round6/trace_envelopes/verify_results.py
```

The numerical proposer uses SciPy/HiGHS; exact reconstruction uses SymPy and Python fractions. The saved review can be checked without trusting the proposer. Earlier mathematical outputs are read-only inputs.

## Existing mathematics

Finite-support moment bounds by linear programming are classical: [Prékopa, The discrete moment problem and linear programming (1990)](https://www.sciencedirect.com/science/article/pii/0166218X9090068N) treats fixed-support moment constraints and dual bounds. [Bertsimas–Popescu, Optimal inequalities in probability theory (2005)](https://www.mit.edu/~dbertsim/papers/MomentProblems/Optimal-inequalities-in-probability-theory-A-convex-optimization-approach-SIAM15.pdf) develops the broader optimization approach. Subtracting equality-row components and verifying rational duals are established linear-algebra and optimization techniques.

Earlier local rounds already used realization fibers and norm envelopes. The specific addition is the residue-sensitive conference model, its exact sufficient-statistic interface, the measured critical-scale sensitivity, and portable ambiguity/uniqueness certificates. No historical-first claim is supported.
