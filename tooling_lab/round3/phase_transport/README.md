# Phase transport microscope: a controlled negative result

This prototype asks whether ordinary translation correlations missed useful modulation information in the saved Paley witnesses. It computes the complete finite Weyl ambiguity plane at primes 101, 401, and 1009. **The tested generic phase-plane summaries do not provide a new separator.** The strongest apparent signals are amplitude statistics or the already measured translation axis; an initial fourth-moment excess largely disappears after the null controls preserve the known S3 action.

The reusable output is a verified phase-resolved interface, complete complex data, and controls that distinguish amplitude concentration, known symmetry, and actual phase information. This is a local experimental tool built from established mathematics, not an invention claim, an inverse theorem, an upper spectral bound, or a proof of either prize problem.

## What is tested

Inputs are the six saved round2 integer witnesses `witness_{p}_{positive|negative}.bin`, with their SHA256 hashes recorded in `results.json`. The vectors are supported on

\[
C=\{x\in\mathbb F_p:\chi(x)=\chi(x-1)=1\},\qquad
H_{xy}=\chi(x-y)/\sqrt p-1/p\quad(x,y\in C).
\]

The original vectors came from separate positive and negative one-sided spectral searches. They are **selected near-edge witnesses**, not six examples above the putative asymptotic edge. Only the negative p401 witness among these six exceeds \(\sqrt3/2\) in magnitude. Its integer Rayleigh certificate is in round2; this experiment does not replace it with a numerical certificate.

Write \(f=z/\sqrt n\), extended by zero outside C, where \(n=\sum z_x^2\). The full phase plane is

\[
D(a,b)f(x)=e_p\bigl(b(x+a/2)\bigr)f(x+a),\qquad
A_f(a,b)=\langle f,D(a,b)f\rangle.
\]

All arithmetic inside \(e_p(t)=\exp(2\pi i t/p)\) is modulo p, including \(1/2\). Arrays use `A[a,b]`. The symmetric phase fixes the convention, so

\[
D(a,b)D(c,d)=e_p((ad-bc)/2)D(a+c,b+d),\qquad
A_f(-a,-b)=\overline{A_f(a,b)}.
\]

The b=0 slice is exactly the previous ordinary cyclic correlation. Its integer numerator is recomputed at **every shift** by Python integer sums. The nonzero modulation slices are floating-point complex FFT computations, not exact cyclotomic certificates.

## Identities that prevent false discoveries

For every unit vector, independently of any Paley structure,

\[
A_f(0,0)=1,\quad
\sum_{a,b}|A_f(a,b)|^2=p,\quad
\sum_b|A_f(a,b)|^2=p\sum_x|f(x)|^2|f(x+a)|^2.
\]

The last formula shows that a row's total phase energy is determined entirely by amplitudes. In particular,

\[
\sum_{b\ne0}|A_f(0,b)|^2=p\sum_x|f(x)|^4-1.
\]

Thus the vertical line is a coordinate-space participation statistic. A coefficient permutation preserves that line energy exactly. A change of phases with fixed amplitudes preserves the entire a=0 slice exactly.

The p+1 projective lines through the origin partition the nonzero phase points. Their energies excluding the origin always average \((p-1)/(p+1)\). For finite slope s, the energy including the origin is

\[
\sum_a|A_f(a,sa)|^2
=p\sum_k\left|\mathcal F\left[f(x)e_p(sx^2/2)\right](k)\right|^4,
\]

where \(\mathcal F\) is the unitary Fourier transform with negative exponent. This is an ordinary Parseval identity in the known chirp/mutually unbiased basis. `largest_line_slope=p` denotes the vertical line a=0; slopes 0 through p−1 denote b=sa.

For a chirp \(f_c(x)=e_p(cx^2)f(x)\),

\[
A_{f_c}(a,b)=A_f(a,b+2ca).
\]

Chirping can therefore change the ordinary translation slice while preserving the global maximum, total fourth moment, and maximum projective-line energy. Under the unitary Fourier transform,
\(A_{\mathcal Ff}(a,b)=A_f(b,-a)\).
These are standard symplectic covariance properties, not new tools discovered here.

## Initial experiment and measured revision

The pilot compared each original vector with eight seeded coefficient shuffles within C, one random-phase vector with the same amplitudes, and one constant-amplitude chirp masked to C. Full-field chirps and delta functions supply closed-form controls: their ambiguity planes are respectively 1 on b=−2a and a=0, and 0 elsewhere. The fourth-moment statistic is explicitly \(p^{-1}\sum_{a,b}|A_f(a,b)|^4\); it is related to established stabilizer entropy quantities and is not presented as a new invariant.

| p | side | H quotient | max ordinary overlap | max full phase overlap | location of full maximum |
|---:|:---|---:|---:|---:|:---|
| 101 | positive | 0.604850 | 0.276585 | 0.539799 | (0,45) |
| 101 | negative | −0.777099 | 0.391187 | 0.650267 | (0,36) |
| 401 | positive | 0.799576 | 0.246808 | 0.379812 | (0,86) |
| 401 | negative | −0.877493 | 0.330214 | 0.475968 | (0,121) |
| 1009 | positive | 0.819481 | 0.222847 | 0.332480 | (0,345) |
| 1009 | negative | −0.858021 | 0.215940 | 0.386184 | (0,69) |

Every larger full-plane maximum occurs on the amplitude-only a=0 axis. Every largest projective-line energy is also vertical. The pilot's maximum-line statistic was unchanged by all coefficient shuffles, so it cannot separate those witnesses from those controls.

The revision removes a=0 and b=0 from the phase statistic and strengthens the null model. C carries the already known S3 action generated by x↦1−x and x↦1/x. Arbitrary coordinate permutations destroy this symmetry. `refine_phase.py` instead permutes regular six-element orbits equivariantly, fixing exceptional smaller orbits. Each permutation is checked for bijectivity and exact commutation with both generators. This preserves the coefficient histogram and the vector's transformation behavior under the known action. It does not preserve H or the Rayleigh quotient, which is explicitly recomputed.

| p | side | off-axis fourth-power sum | eight S3-null range | S3-null median |
|---:|:---|---:|:---|---:|
| 101 | positive | 2.923125 | 2.9386–3.6541 | 3.1697 |
| 101 | negative | 3.072897 | 2.5846–5.1724 | 3.1640 |
| 401 | positive | 3.152932 | 3.0075–3.2241 | 3.1355 |
| 401 | negative | 3.177875 | 3.1520–3.8259 | 3.3014 |
| 1009 | positive | 3.107443 | 3.0613–3.3227 | 3.2243 |
| 1009 | negative | 3.175493 | 3.0172–3.2322 | 3.1093 |

There is no consistent fourth-moment excess against the stronger controls. The largest nonvertical line in five of six cases is simply b=0, the ordinary translation aggregate. This does not establish that all phase information is unhelpful; it rejects the interpretation of these particular generic summaries as evidence of a new phase signature. The controls are eight exploratory seeded samples, not a calibrated statistical test or a complete conditional ensemble.

A second diagnostic keeps H fixed and chirps f. At p401 on the negative witness, c=1 changes the quotient from −0.8774933056 to −0.0130228238 while leaving all three global shear-invariant summaries unchanged. Thus those summaries cannot *suffice* to identify a large quotient for this fixed operator. This does not rule out a necessary condition, operator-aligned phase geometry, or a diagnostic that also uses information these summaries discard. Chirped vectors are allowed complex vectors on the same support; this example must not be mistaken for a real-vector-only separation claim.

![Phase planes and matched controls](phase_controls.png)

The two top panels use the same absolute-overlap color scale; the identity entry is omitted. The cyan a=0 line contains amplitude-only information. The bottom panels show the symmetry-preserving nulls and the change in the fixed H quotient under a chirp. The global summaries preserved in the latter comparison are the maximum overlap, total fourth moment, and maximum projective-line energy; individual named coordinates and the off-axis fourth moment need not be preserved.

## Verification and complexity

`validation.json` records an additional implementation that checks 32 independently computed direct complex sums and 16 Weyl group products per witness. Maximum direct-sum error is below 1.2e−16; group errors are below 5.1e−16. It also checks the Wigner transform

\[
W_f(q,k)=\frac1p\sum_a f(q+a/2)\overline{f(q-a/2)}e_p(-ka)
\]

for realness, position and Fourier marginals, and \(\sum W_f^2=1/p\). All reported errors are below 7e−17. Fourier covariance is checked over the entire p101 plane, while chirp covariance is checked over every plane in the pilot. Wigner \(\ell^1\) values are exported only as numerical diagnostics; no Wigner positivity theorem is asserted about these vectors.

The root agent added a separate scalar audit, `review_phase_root.py` / `root_review.json`. It checks every ambiguity and Wigner entry of six generic complex vectors at p=7 and p=13 (654 phase points), independently checks Fourier covariance, and recomputes the two p401 fixed-operator chirp quotients by scalar sums. All checks pass, with maximum error below 5.6e−16. This extends the real-witness audit to generic complex phase conventions; it remains a numerical validation rather than an exact nonzero-phase certificate.

One full plane uses p FFTs of length p: O(p² log p) work and O(p²) complex storage. At p1009 one complex plane is 16,289,296 bytes; the routines hold several such arrays. A hard guard refuses p>1009. This is a small-field microscope, not a million-prime algorithm. Computing k chosen displacement rows alone would require O(k p log p) work and O(p) streaming workspace, but that extension is not implemented here. Exact b=0 verification uses O(p²) integer products; the prior round's exact NTT implementation already supplies the scalable ordinary slice.

Reproduce from this repository with the runtime recorded in `results.json`:

```sh
OPENBLAS_NUM_THREADS=1 /opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round3/phase_transport/ambiguity_lab.py
OPENBLAS_NUM_THREADS=1 /opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round3/phase_transport/refine_phase.py
OPENBLAS_NUM_THREADS=1 /opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round3/phase_transport/verify_phase.py
OPENBLAS_NUM_THREADS=1 /opt/miniconda3/bin/python tooling_lab/round3/phase_transport/review_phase_root.py
OPENBLAS_NUM_THREADS=1 /opt/miniconda3/bin/python tooling_lab/round3/phase_transport/export_results.py
```

The pilot took 18.89 seconds. Solver runtime: Python 3.14.6, NumPy 2.4.4, SciPy 1.17.1. The plotting runtime is separately recorded in `summary.json`: Python 3.13.2, NumPy 2.4.5, Matplotlib 3.10.3. The inputs and sources are hashed. Saved `phase_plane_*.npz` files contain `original` and its c=1 `chirped` complex planes; JSON files retain all p+1 line energies, exact b=0 numerators, null summaries, and validation errors. `summary.json` and `summary.csv` provide compact reusable comparisons. The final SVG and PNG were visually inspected.

## Prior art and the full-problem gap

The primary-source comparison used:

- Gurevich, Hadani, and Sochen, *The finite harmonic oscillator and its associated sequences*: finite Heisenberg/Weil representations, ambiguity functions, chirps, and low-correlation sequences. [Author-hosted paper](https://people.math.wisc.edu/~sgurevich/PNAS-Version.full.pdf), [construction supplement](https://arxiv.org/html/0808.1417v1).
- Gross, *Hudson's Theorem for finite-dimensional quantum systems*: odd-dimensional discrete Wigner functions, stabilizer states, and Clifford covariance. [Primary paper](https://arxiv.org/abs/quant-ph/0602001).
- Leone, Oliviero, and Hamma, *Stabilizer Rényi entropy*: moments of Pauli expectation values as established nonstabilizerness diagnostics. We define our own displayed normalization and do not silently transfer qubit theorems to odd-prime dimensions. [Primary paper](https://arxiv.org/abs/2106.12587).
- Chaturvedi, *Aspects of mutually unbiased bases in odd prime power dimensions*: the established finite-field chirp/MUB setting. [Primary paper](https://arxiv.org/abs/quant-ph/0109003).

A bounded local text search found no previous ambiguity/Weyl/chirp/Wigner implementation in `research`, `experiments`, or `tooling_lab` before this directory was created. That is a local search result, not proof that this application has never been tried anywhere. No mathematical foundation used here is claimed as newly invented.

The next substantive missing step is an **operator-aligned** phase constraint: a feature must connect the specific localized character operator with a signed Rayleigh excess, after accounting for amplitude statistics and known S3 symmetry. Merely finding a large ambiguity coefficient, a concentrated projective line, or a low phase entropy does not supply that connection. This experiment establishes no uniform implication, no upper bound on all vectors, and no asymptotic improvement. Its value is avoiding a misleading phase-only route while providing the measured interface needed to test a more specific one.
