# Arithmetic lift diagnostic — 2026-09-05

This experiment resolves a specific failure of round one's syzygy diagnostic. The known F41 sextic configuration and the structural binomial/coset examples all lose their kernel after every permitted one-root replacement. Their first-order ambient deformation signatures can also agree. Those tests do not distinguish the origin of their rank defect.

The new diagnostic changes characteristic precision while **keeping every root on the same roots-of-unity domain**. It asks which mod-p syzygies can lift to mod-p², mod-p³ and mod-p⁴. This uses established Hensel lifting and local Smith theory, not a new general mathematical theory. Its proposed contribution is a domain-preserving matrix-kernel test, joined to the existing obstruction examples, with explicit lift-convention controls and portable exact certificates.

## Observed distinction

| Configuration | Kernel survival at p,p²,p³,p⁴ using genuine domain lifts | Survival at p² using naive residue representatives |
|---|---|---:|
| F41, μ20, sextics with degree-one cofactors | 1,0,0,0 | 0 |
| Three μ4 cosets inside μ16, p=17 | 1,1,1,1 | 0 |
| Same coset construction, p=1009 | 1,1,1,1 | 0 |
| Three μ6 cosets inside μ24, p=97 | 1,1,1,1 | 0 |
| Three μ12 cosets inside μ48, p=193 | 1,1,1,1 | 0 |

The F41 finite-field syzygy has an exact nonzero obstruction already modulo 41². Its local diagonal valuations at that precision are `(0,0,0,0,0,1)`. Thus no nonzero initial syzygy lifts. Each structural control retains its kernel through all tested precisions. Their persistence is expected from the binomial identity, independently of this finite test: lifted cosets still have polynomials `X^d-c_i`, and three such polynomials have the constant syzygy with coefficients `(c_2-c_3,c_3-c_1,c_1-c_2)`.

The convention control is decisive: simply retaining each integer residue from 0 to p-1 kills **every structural control** at p². Such integer representatives usually fail `x^n=1 mod p²`. A lift-fragility score computed on those representatives would misclassify a field-uniform construction. The tool therefore requires the domain equation to be independently checked at every precision.

This is a finite separation of examples, not a theorem that every accidental kernel fails at p² or every persistent kernel has coset form. Deeper accidental p-adic coincidences can exist, and finite survival does not by itself prove infinite lifting. No link to the magnitude of complex Gauss sums or any prize bound is inferred.

## Interface and elementary justification

For p∤n and a root a of `X^n-1 mod p`, the simple-root Hensel step lifts it uniquely to a root at each p^e. Given three lifted root sets, form the same cofactor coefficient matrix A_e as in round one. We compute invertible U_e,V_e satisfying

```
U_e A_e V_e = diag(p^v1,...,p^vr,0,...,0) mod p^e,
where every v_i<e.
```

The image in Fp of `ker(A_e mod p^e)` has dimension `number_of_columns-r`. Indeed, a nonzero diagonal entry forces its coordinate to be divisible by p, whereas every zero diagonal coordinate has arbitrary reduction modulo p; V_e is invertible modulo p. These are the **initial** syzygies that survive, not the total number of modular solutions. The latter also includes solutions divisible by p that reduce to zero.

At the first step a shorter obstruction is available. For an initial kernel vector c, choose its canonical coefficient lift and set `b=A_2 c/p mod p`. An extension `c+p*d` exists precisely when `A_1 d=-b`. Any left-kernel vector l with `l*b !=0` certifies failure. The saved F41 record contains c,l,b and the nonzero pairing. A different coefficient representative changes b by an image vector, so this obstruction class is well-defined for the chosen lifted matrix.

The diagonalization uses minimum-valuation pivots and inverts units only. All operations preserve invertibility; their transformations are saved and checked directly. The method is ordinary local Smith elimination, with Fitting ideals providing the presentation-independent interpretation of the minor ideals.

## Local prior art and a recovered convention mismatch

The local file `/Users/shawwalters/proximityprize/scripts/probes/sw1_transv_multiplicity_profile.py` already Hensel-lifts roots to evaluate the valuations of **scalar** signed root-sums across degree-one cyclotomic primes. Hensel lifting itself and valuation profiling are therefore already implemented locally. The added object here is the whole cofactor-kernel image, with matrix-level obstruction certificates and the domain-preservation audit.

`docs/kb/deltastar-466-ten-by-ten-paley-assault-2026-07-10.md` already proposes Hensel stratification of collision varieties. `docs/kb/NEW_MATHEMATICS-2026-06-15.md` mentions Smith-form torsion primes. The new executable experiment must not be presented as the invention of these approaches.

The historical note `docs/kb/deltastar-444-CHARP-TRANSFER-NEW-GROUND-2026-06-22.md` uses the example `1+g^(n/2)=p` alongside a discussion of Teichmüller lifts. This value applies to canonical integer representatives `1+(p-1)`. Under the actual domain-preserving lift for even n, `g^(n/2)=-1 mod p^e`, hence the sum is zero at every precision. The prototype checks this distinction at p=17,41,97 through p³. It invalidates using that particular representative computation as a Teichmüller identity. This does not validate or refute every broader assertion in the old note; those were not independently audited here. The source files were not edited.

## Checks and files

`arithmetic_lift.py` generates five profiles, with twenty saved U*A*V certificates. Its projected-kernel implementation was exhaustively compared with all 256 two-by-two matrices over Z/4, enumerating every possible kernel vector and its reduction modulo two.

`verify_arithmetic_lift.py` imports no discovery code. It reconstructs all coefficient matrices from the roots, verifies 384 root/domain/congruence checks, independently checks invertibility modulo p, verifies every modular matrix identity, and replays the nonlifting witness. All checks pass; results are in `arithmetic_lift_results.json` and `arithmetic_lift_validation.json`.

Run from any directory:

```sh
python3 /Users/shawwalters/Desktop/paleygraph/tooling_lab/round2/novelty/arithmetic_lift.py
python3 /Users/shawwalters/Desktop/paleygraph/tooling_lab/round2/novelty/verify_arithmetic_lift.py
```

The generator reuses round-one finite-field coefficient routines; the separate verifier is standard-library only. The current matrices have at most 13 rows and six columns. Nothing here scales the full prize instance.

## What to test next

Freeze the diagnostic before selecting a broader family. Collect independently generated special-field syzygies and classify lift death depth, surviving rank, cofactor degree and exact stack incidence. Include non-coset constructions that lift, accidental examples surviving p² but not p³, and different characteristic embeddings of the same symbolic pattern. A falsifiable hypothesis would predict a specific relation between those independently measured features; no such relation is assumed yet.

A useful next algorithm would combine sparse local Smith computation with an exact support-family generator, while retaining a certificate of which initial dependencies survive. A new theorem would need to control the distribution or classification of these defects, not just compute them. The present result justifies **using lawful arithmetic lifts as a discriminator before interpreting deformation stability**.

## Primary-source and query ledger

All sources were searched/read on 2026-09-05. This is a bounded comparison; global novelty remains unestablished.

- [Stacks Project, Henselian local rings, tag04GE](https://stacks.math.columbia.edu/tag/04GE): simple-root lifting and uniqueness. Read the definition and uniqueness argument. This is the standard base of the root lift.
- [Stacks Project, Fitting ideals, tag07Z6](https://stacks.math.columbia.edu/tag/07Z6): ideals of presentation minors and independence from presentation. This is the standard invariant-theoretic interpretation; no Fitting-ideal invention is claimed.
- Elsheikh–Giesbrecht–Novocin–Saunders, [Fast Computation of Smith Forms of Sparse Matrices Over Local Rings](https://arxiv.org/abs/1201.5365), initial2012-01-25: efficient local-ring Smith computation, including Z/p^e. Abstract and metadata inspected. It is a concrete scaling source, not a theorem already implemented by this small dense prototype.
- [Relating p-adic eigenvalues and the local Smith normal form](https://arxiv.org/abs/1401.1773): adjacent prior art located; no result from it is used. Eigenvalue valuations and Smith valuations require hypotheses and must not be identified automatically.

Literal queries:

1. `site.stacks.math.columbia.edu Hensel simple root etale lifting unique modulo ideal`
2. `site.arxiv.org Smith normal form local rings p-adic matrix lifting kernel modulo prime power`
3. `site.arxiv.org Reed Solomon syzygy p-adic lifts roots unity Smith normal form`
4. `site.stacks.math.columbia.edu fitting ideals minors base change presentation module`
5. `"syzygy" "Teichmuller" "Smith"`
6. `"Reed-Solomon" "Hensel" "Smith"`

Local searches included `Teichm`, `Hensel`, `Witt vector`, `Smith normal`, `Fitting ideal`, p-adic/lift/syzygy/rank combinations, and mod-p² lifting phrases across both projects' probes, research notes and formalization frontier. No exact prior matrix-survival plus lift-convention experiment was located in that limited search. That is a search result, not proof of nonexistence.
