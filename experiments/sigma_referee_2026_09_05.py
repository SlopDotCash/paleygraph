#!/usr/bin/env python3
"""sigma pass, direction `referee` (2026-09-05): exact numerical audit of

 (a) research/parallel9-all-degrees-2026-09-05.md, claim (T):
       |tr(D_{Z_1} S ... D_{Z_k} S)| <= 3a(2a+2)^{k-2} p^{(k+1)/2}   (k>=3, p = 1 mod 4,
       a anchors, nonempty labels Z_i inside the anchor set),
 (b) research/parallel10-word-aggregate-2026-09-05.md, claims (4)/(10) (simultaneous Gram
       bound for raw word kernels) and (13) (dimension obstruction).

Parts (all arithmetic exact; floats appear only inside FFT digit planes whose rounding
residual is asserted < 0.1, and in reported ratios/slopes):
 A. dense exact necklace traces, a in {1,2,3}, k in 1..8, p = 1 mod 4 up to ~3400;
    exact comparison with (T), with the trivial bound p^{k/2+1}, log-log slopes,
    non-vacuity thresholds p_0(a,k) = (3a(2a+2)^{k-2})^2.
 B. exact fixed-variable inner sums I_k(y) (the (k-1)-variable sums that the note bounds
    by (2a-1)(2a+2)^{k-2} p^{(k-1)/2}) at p ~ 1e4, 1e5, 1e6 by integer FFT convolution.
 C. Task 3: the exact indicator expansion E_Z = 2^{-|Z|} sum_{W subset Z} D_W - B_Z, the
    resulting trace expansion, the clique case (parallel2 (8)), and S^2 = pI - J.
 D. independent implementation of Katz RLS 3.3.6/3.3.7 (tame local monodromy and rank of
    middle convolution): checks of the rank-growth lemma (R)/(R'), the pseudoreflection
    at y, the closed formulas for block counts, the involution, word recovery (b)(1)-(2),
    the raw-rank recurrence (6) and the word sums (7)-(8).
 E. raw-kernel Gram matrices on U (parallel10 (10)): diagonal -> 1 (irreducibility and
    self-duality test), off-diagonal O(1/sqrt p) (non-isomorphism test), non-asymptotic
    bound in non-vacuous cases; also k=1 (zero) and k=2 (a^2 p) checks.
 F. summary of the previous worker's larger-p data (scratchpad, optional).

Writes results/sigma_referee_2026_09_05.json. Target runtime < 10 minutes.
"""
import json, math, os, sys, time, random, itertools
from collections import Counter
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'results', 'sigma_referee_2026_09_05.json')
SCRATCH = ('/private/tmp/claude-501/-Users-shawwalters-Desktop-paleygraph/'
           '1b733267-6e28-4860-bf0c-a4c3dcfdbfb5/scratchpad/referee')
FAST = os.environ.get('SIGMA_FAST', '0') == '1'
T0 = time.time()
WITNESSES = []          # every violated inequality lands here
COUNTS = Counter()      # counts of checks performed


def log(*a):
    print(f'[{time.time()-T0:7.1f}s]', *a, flush=True)


def check(name, ok, witness=None):
    COUNTS[name] += 1
    if not ok:
        WITNESSES.append(dict(check=name, witness=witness))
        log('VIOLATION', name, witness)
    return ok


# ---------------------------------------------------------------- basic finite field tools
def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f, r = 3, int(n ** 0.5)
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


def primes_1mod4(lo, hi):
    return [n for n in range(max(lo, 5), hi + 1) if n % 4 == 1 and is_prime(n)]


def next_prime_1mod4(n):
    while not (n % 4 == 1 and is_prime(n)):
        n += 1
    return n


def chi_table(p):
    chi = -np.ones(p, dtype=np.int64)
    chi[0] = 0
    chi[(np.arange(1, p, dtype=np.int64) ** 2) % p] = 1
    return chi


def S_matrix(p, chi):
    idx = (np.arange(p)[:, None] - np.arange(p)[None, :]) % p
    return chi[idx]


def mask_vec(p, chi, Z):
    """d_Z(x) = prod_{z in Z} chi(x - z), Z a tuple of field elements (may be empty)."""
    x = np.arange(p)
    d = np.ones(p, dtype=np.int64)
    for z in Z:
        d *= chi[(x - z) % p]
    return d


def bound_T_sq(a, k):
    """(3a(2a+2)^{k-2})^2, so that (T) reads N^2 <= bound_T_sq * p^{k+1} exactly in integers."""
    return (3 * a * (2 * a + 2) ** (k - 2)) ** 2


# ---------------------------------------------------------------- Part A: dense exact traces
def prefix_traces_float(p, chi, S, Zs):
    """Exact integer prefix traces tr(D_{Z_1} S ... D_{Z_j} S), j = 1..len(Zs), by float64 matmul.

    Exactness: every entry of M_j = prod_{i<=j} D_{Z_i} S is an integer of absolute value at most
    ||M_j||_op <= p^{j/2}; every partial sum in forming M_j = (M_{j-1} D_{Z_j}) S has at most p terms
    of size <= p^{(j-1)/2}, so all intermediate values are integers of size < p^{(j+1)/2}. We assert
    p^{(j+1)/2} < 2^53, under which float64 addition and multiplication of integers is exact."""
    Sf = S.astype(np.float64)
    M, out = None, []
    for j, Z in enumerate(Zs, 1):
        assert p ** ((j + 1) / 2) < 2 ** 53, (p, j)
        d = mask_vec(p, chi, Z).astype(np.float64)
        M = d[:, None] * Sf if M is None else (M * d[None, :]) @ Sf
        diag = np.diag(M)
        r = np.rint(diag)
        assert np.max(np.abs(diag - r)) < 1e-6
        out.append(int(r.astype(np.int64).sum()))
    return out


def prefix_traces_int(p, chi, S, Zs, dtype=np.int64):
    M, out = None, []
    for Z in Zs:
        d = mask_vec(p, chi, Z).astype(dtype)
        DS = d[:, None] * S.astype(dtype)
        M = DS if M is None else M @ DS
        out.append(int(np.trace(M)))
    return out


PATTERNS = [  # (name, a, word as anchor-index tuples of length >= 8, or None for random)
    ('a1_const', 1, [(0,)] * 8),
    ('a2_full', 2, [(0, 1)] * 8),
    ('a2_alt', 2, [(0,), (1,)] * 4),
    ('a2_0full', 2, [(0,), (0, 1)] * 4),
    ('a2_00_1', 2, ([(0,), (0,), (1,)] * 3)[:8]),
    ('a2_rand', 2, None),
    ('a3_full', 3, [(0, 1, 2)] * 8),
    ('a3_cyc', 3, ([(0,), (1,), (2,)] * 3)[:8]),
    ('a3_pairs', 3, ([(0, 1), (1, 2), (0, 2)] * 3)[:8]),
    ('a3_rand', 3, None),
    ('a3_rand2', 3, None),
]


def random_word(a, seed, length=8):
    rng = random.Random(seed)
    return [tuple(sorted(rng.sample(range(a), rng.randint(1, a)))) for _ in range(length)]


def pattern_words():
    words = {}
    for name, a, w in PATTERNS:
        words[name] = (a, w if w is not None else random_word(a, hash(name) % 1000 + 7))
    return words


def slope(xs, ys):
    if len(xs) < 3:
        return None
    return float(np.polyfit(xs, ys, 1)[0])


def part_A():
    log('Part A: dense exact necklace traces')
    # exactness cross-checks: float path vs int64 path (small p) vs object ints (tiny p)
    rng = random.Random(1)
    for p in (13, 17, 29, 53, 101, 197):
        chi = chi_table(p)
        S = S_matrix(p, chi)
        for _ in range(4):
            a = rng.randint(1, 3)
            anchors = rng.sample(range(p), a)
            Zs = [tuple(anchors[i] for i in sorted(rng.sample(range(a), rng.randint(1, a)))) for _ in range(6)]
            t_float = prefix_traces_float(p, chi, S, Zs)
            t_int = prefix_traces_int(p, chi, S, Zs)
            check('A.float_vs_int64_prefix_traces', t_float == t_int, dict(p=p, Zs=Zs))
            if p <= 29:
                t_obj = prefix_traces_int(p, chi, S, Zs, dtype=object)
                check('A.int64_vs_object_prefix_traces', t_int == t_obj, dict(p=p, Zs=Zs))
    # literal k-fold sum at p=13 for two words (definition check of Kunisky's Definition 1.13)
    p = 13
    chi = chi_table(p)
    S = S_matrix(p, chi)
    for Zs in ([(0,), (0,), (0,)], [(0, 5), (5,), (0,)], [(0,), (5,), (0, 5), (5,)]):
        lit = 0
        k = len(Zs)
        for xs in itertools.product(range(p), repeat=k):
            t = 1
            for i in range(k):
                t *= int(chi[(xs[(i + 1) % k] - xs[i]) % p])
                for z in Zs[i]:
                    t *= int(chi[(xs[i] - z) % p])
            lit += t
        check('A.literal_definition_1.13_equals_trace', lit == prefix_traces_int(p, chi, S, Zs)[-1], dict(p=p, Zs=Zs, lit=lit))

    words = pattern_words()
    cand = primes_1mod4(100, 3450)
    # ~18 primes, roughly geometric spacing
    targets = [int(101 * (3400 / 101) ** (i / 17)) for i in range(18)]
    plist = sorted({min(cand, key=lambda q: abs(q - t)) for t in targets})
    if FAST:
        plist = [q for q in plist if q < 1200]
    rows = []
    for p in plist:
        chi = chi_table(p)
        S = S_matrix(p, chi)
        prng = random.Random(p)
        anchors_by_a = {1: [0], 2: [0, prng.randrange(1, p)], 3: [0] + prng.sample(range(1, p), 2)}
        kmax = 8 if p ** 4.5 < 2 ** 53 else 7
        for name, (a, w) in words.items():
            anchors = anchors_by_a[a]
            Zs = [tuple(anchors[i] for i in T) for T in w[:kmax]]
            tr = prefix_traces_float(p, chi, S, Zs)
            row = dict(p=p, name=name, a=a, anchors=anchors, traces=tr, nonvacuous=[], viol=[])
            check('A.k1_trace_zero', tr[0] == 0, dict(p=p, name=name))
            check('A.k2_bound_a2p', tr[1] ** 2 <= a ** 4 * p ** 2, dict(p=p, name=name, N=tr[1]))
            for k in range(3, kmax + 1):
                N = tr[k - 1]
                ok = N * N <= bound_T_sq(a, k) * p ** (k + 1)
                check('A.bound_T', ok, dict(p=p, name=name, k=k, N=N))
                check('A.trivial_bound', N * N <= p ** (k + 2), dict(p=p, name=name, k=k, N=N))
                if bound_T_sq(a, k) < p:
                    row['nonvacuous'].append(k)
                if not ok:
                    row['viol'].append(k)
            rows.append(row)
        log(f'  p={p} done; kmax={kmax}')
    # regression summaries
    summary = {}
    for name, (a, w) in words.items():
        for k in range(3, 9):
            pts = [(r['p'], abs(r['traces'][k - 1])) for r in rows if r['name'] == name and len(r['traces']) >= k and r['traces'][k - 1] != 0]
            xs = [math.log(q) for q, _ in pts]
            ys = [math.log(v) for _, v in pts]
            big = [(q, v) for q, v in pts if q >= 1000]
            summary[f'{name}:k={k}'] = dict(
                a=a, k=k, n_nonzero=len(pts),
                loglog_slope=slope(xs, ys), claimed_exponent=(k + 1) / 2, generic_exponent=k / 2,
                median_ratio_to_p_half_k_plus_1=float(np.median([v / q ** ((k + 1) / 2) for q, v in big])) if big else None,
                max_ratio_to_p_half_k_plus_1=max([v / q ** ((k + 1) / 2) for q, v in big], default=None),
                median_ratio_to_p_half_k=float(np.median([v / q ** (k / 2) for q, v in big])) if big else None,
                max_ratio_to_bound_T=max([v / (math.sqrt(bound_T_sq(a, k)) * q ** ((k + 1) / 2)) for q, v in pts], default=None),
                nonvacuous_primes_tested=[r['p'] for r in rows if r['name'] == name and k in r['nonvacuous']])
    thresholds = {f'a={a},k={k}': bound_T_sq(a, k) for a in (1, 2, 3) for k in range(3, 9)}
    return dict(primes=plist, rows=rows, summary=summary, nonvacuity_threshold_p0=thresholds)


# ---------------------------------------------------------------- Part B: exact FFT inner sums
BITS = 15
BASE = 1 << BITS


def to_planes(v):
    """Signed base-2^15 digit planes of an integer vector (int64 or object)."""
    if v.dtype == object:
        rem = [int(x) for x in v]
        planes = []
        while any(rem):
            lo = [((x + BASE // 2) % BASE) - BASE // 2 for x in rem]
            planes.append(np.array(lo, dtype=np.int64))
            rem = [(x - l) >> BITS for x, l in zip(rem, lo)]
        return planes or [np.zeros(len(v), dtype=np.int64)]
    rem = v.astype(np.int64).copy()
    planes = []
    while np.any(rem):
        lo = ((rem + BASE // 2) % BASE) - BASE // 2
        planes.append(lo)
        rem = (rem - lo) >> BITS
    return planes or [np.zeros(len(v), dtype=np.int64)]


class Convolver:
    """Exact circular convolution with chi via float64 FFT on digit planes (each plane's
    convolution is an integer of size <= p * 2^14 ~ 2^35 for p <= 2^21, far inside float64
    exactness once the FFT rounding residual is checked)."""

    def __init__(self, p, chi):
        self.p = p
        self.fchi = np.fft.rfft(chi.astype(np.float64))

    def conv(self, v):
        p = self.p
        outs = []
        for pl in to_planes(v):
            r = np.fft.irfft(np.fft.rfft(pl.astype(np.float64)) * self.fchi, n=p)
            ri = np.rint(r)
            assert np.max(np.abs(r - ri)) < 0.1, 'FFT rounding residual too large'
            outs.append(ri.astype(np.int64))
        if len(outs) * BITS + 40 < 63:
            total = np.zeros(p, dtype=np.int64)
            for i, pl in enumerate(outs):
                total += pl << (BITS * i)
            return total
        total = np.array([0] * p, dtype=object)
        for i, pl in enumerate(outs):
            total = total + np.array([int(x) << (BITS * i) for x in pl], dtype=object)
        return total


def chain_at_y(conv, p, chi, Zs, y, kmax):
    """u_1 = chi(x - y); u_j = chi (*) (d_{Z_j} u_{j-1}).  Returns [(u_k)_y for k = 1..kmax] and u_kmax.
    (u_k)_y = [S D_{Z_2} S ... D_{Z_k} S]_{yy} is the (k-1)-variable sum with x_1 = y fixed, and
    tr(D_{Z_1} S ... D_{Z_k} S) = sum_y d_{Z_1}(y) (u_k)_y."""
    x = np.arange(p)
    u = chi[(x - y) % p].copy()
    vals = [int(u[y])]
    for j in range(2, kmax + 1):
        d = mask_vec(p, chi, Zs[j - 1])
        masked = d * u if u.dtype != object else np.array([int(a) * int(b) for a, b in zip(u, d)], dtype=object)
        u = conv.conv(masked)
        vals.append(int(u[y]))
    return vals, u


def part_B():
    log('Part B: exact inner sums at large p')
    words = pattern_words()
    # cross-check chain sums against dense traces at small p
    for p in (101, 401):
        chi = chi_table(p)
        S = S_matrix(p, chi)
        conv = Convolver(p, chi)
        for name, (a, w) in words.items():
            anchors = [0, 3, 7][:a]
            Zs = [tuple(anchors[i] for i in T) for T in w[:8]]
            dense = prefix_traces_float(p, chi, S, Zs)
            d1 = mask_vec(p, chi, Zs[0])
            total = [0] * 8
            for y in range(p):
                if d1[y] == 0:
                    continue
                vals, _ = chain_at_y(conv, p, chi, Zs, y, 8)
                for k in range(8):
                    total[k] += int(d1[y]) * vals[k]
            check('B.fft_chain_sum_equals_dense_trace', total == dense, dict(p=p, name=name))
    plist = [10009, next_prime_1mod4(100000)] + ([] if FAST else [next_prime_1mod4(1000000)])
    rows = []
    for p in plist:
        chi = chi_table(p)
        conv = Convolver(p, chi)
        prng = random.Random(p + 1)
        anchors_by_a = {1: [0], 2: [0, prng.randrange(1, p)], 3: [0] + prng.sample(range(1, p), 2)}
        n_generic = 2 if p > 500000 else 3
        for name, (a, w) in words.items():
            anchors = anchors_by_a[a]
            Zs = [tuple(anchors[i] for i in T) for T in w[:8]]
            ys = [('generic', prng.randrange(3, p)) for _ in range(n_generic)]
            exc = [v for v in anchors if v not in Zs[0]]
            if exc:
                ys.append(('anchor_not_in_Z1', exc[0]))
            ys.append(('anchor_in_Z1', Zs[0][0]))
            for kind, y in ys:
                kmax = 8 if (p < 500000 or kind == 'generic') else 5
                vals, _ = chain_at_y(conv, p, chi, Zs, y, kmax)
                for k in range(3, kmax + 1):
                    I = vals[k - 1]
                    c = (2 * a - 1) * (2 * a + 2) ** (k - 2)
                    if kind == 'generic':
                        ok = I * I <= c * c * p ** (k - 1)
                        check('B.inner_sum_note_bound_generic_y', ok, dict(p=p, name=name, k=k, y=y, I=I))
                        nonvac = c * c < p
                    else:
                        ok = I * I <= p ** k
                        check('B.inner_sum_trivial_bound_anchor_y', ok, dict(p=p, name=name, k=k, y=y, I=I))
                        nonvac = False
                    rows.append(dict(p=p, name=name, a=a, k=k, kind=kind, y=y, I=I,
                                     ratio_to_p_half_k_minus_1=abs(I) / p ** ((k - 1) / 2),
                                     ratio_to_p_half_k=abs(I) / p ** (k / 2),
                                     note_bound_nonvacuous=nonvac, ok=ok))
            log(f'  p={p} {name} done')
    summary = {}
    for name, (a, w) in words.items():
        for k in range(3, 9):
            pts = [(r['p'], abs(r['I'])) for r in rows if r['name'] == name and r['k'] == k and r['kind'] == 'generic' and r['I'] != 0]
            byp = {}
            for q, v in pts:
                byp[q] = max(byp.get(q, 0), v)
            xs = [math.log(q) for q in sorted(byp)]
            ys = [math.log(byp[q]) for q in sorted(byp)]
            summary[f'{name}:k={k}'] = dict(
                a=a, k=k, loglog_slope_of_max_over_y=slope(xs, ys), claimed_exponent=(k - 1) / 2,
                ratios_to_p_half_k_minus_1={str(q): byp[q] / q ** ((k - 1) / 2) for q in sorted(byp)},
                anchor_y_ratios={f"{r['p']}:{r['kind']}": r['ratio_to_p_half_k_minus_1'] for r in rows if r['name'] == name and r['k'] == k and r['kind'] != 'generic'})
    return dict(primes=plist, rows=rows, summary=summary)


# ---------------------------------------------------------------- Part C: indicator expansion
def part_C():
    log('Part C: indicator expansion (Task 3)')
    rng = random.Random(3)
    rows = []
    for p in (13, 17, 29, 37, 53):
        chi = chi_table(p)
        S = S_matrix(p, chi)
        x = np.arange(p)
        check('C.S_squared_equals_pI_minus_J', np.array_equal(S @ S, p * np.eye(p, dtype=np.int64) - np.ones((p, p), dtype=np.int64)), dict(p=p))
        for trial in range(6):
            a = rng.randint(1, 3)
            anchors = rng.sample(range(p), a)
            k = rng.randint(2, 4)
            Zs = [tuple(sorted(rng.sample(anchors, rng.randint(1, a)))) for _ in range(k)]
            # 2^{|Z|} E_Z and 2^{|Z|} B_Z as integer vectors
            def E_scaled(Z):
                e = np.ones(p, dtype=np.int64)
                for z in Z:
                    e *= 1 + chi[(x - z) % p]
                for z in Z:
                    e[z] = 0
                return e
            def B_scaled(Z):
                b = np.zeros(p, dtype=np.int64)
                for z in Z:
                    if all(chi[(z - z2) % p] == 1 for z2 in Z if z2 != z):
                        b[z] = 1 << (len(Z) - 1)
                return b
            for Z in Zs:
                rhs = sum(mask_vec(p, chi, W) for r_ in range(len(Z) + 1) for W in itertools.combinations(Z, r_)) - B_scaled(Z)
                check('C.pointwise_E_Z_expansion', np.array_equal(E_scaled(Z), rhs), dict(p=p, Z=Z))
                clique = all(chi[(z1 - z2) % p] == 1 for z1, z2 in itertools.combinations(Z, 2))
                if clique:
                    # parallel2 (8): B_Z = (1/2) 1_Z  when Z is a clique
                    check('C.clique_boundary_is_half_indicator', all(B_scaled(Z)[z] == (1 << (len(Z) - 1)) for z in Z), dict(p=p, Z=Z))
            # trace expansion: 2^{sum|Z_i|} tr(prod E_{Z_i} S) = sum_{W_i subset Z_i} Sigma(W) + boundary
            M = np.eye(p, dtype=np.int64)
            for Z in Zs:
                M = M @ (E_scaled(Z)[:, None] * S)
            lhs = int(np.trace(M))
            pure = 0
            for Ws in itertools.product(*[[W for r_ in range(len(Z) + 1) for W in itertools.combinations(Z, r_)] for Z in Zs]):
                pure += prefix_traces_int(p, chi, S, list(Ws))[-1]
            # pure part must equal 2^{sum|Z|} tr(prod (E+B) S) by multilinearity
            M2 = np.eye(p, dtype=np.int64)
            for Z in Zs:
                M2 = M2 @ ((E_scaled(Z) + B_scaled(Z))[:, None] * S)
            check('C.pure_part_equals_E_plus_B_chain', pure == int(np.trace(M2)), dict(p=p, Zs=Zs))
            boundary = lhs - pure
            rows.append(dict(p=p, anchors=anchors, Zs=Zs, scaled_trace=lhs, pure=pure, boundary=boundary,
                             all_clique=all(chi[(z1 - z2) % p] == 1 for z1, z2 in itertools.combinations(anchors, 2))))
            COUNTS['C.trace_expansions'] += 1
    return dict(rows=rows)


# ---------------------------------------------------------------- Part D: Katz local monodromy bookkeeping
# A tame local representation with eigenvalues +-1 is a Counter {(eigen, blocksize): multiplicity}.
def rdim(rep):
    return sum(l * n for (e, l), n in rep.items())


def nblocks(rep, e):
    return sum(n for (s, l), n in rep.items() if s == e)


def n1blocks(rep, e):
    return rep.get((e, 1), 0)


def twist(rep):
    return Counter({(-e, l): n for (e, l), n in rep.items()})


def quot_inv(rep):
    """V / V^I: unipotent blocks shrink by one (size-one ones vanish), quadratic blocks unchanged."""
    out = Counter()
    for (e, l), n in rep.items():
        if e == 1:
            if l > 1:
                out[(1, l - 1)] += n
        else:
            out[(e, l)] += n
    return out


def reconstruct(q, dim):
    """The unique V with V/V^I = q and dim V = dim (unipotent blocks grow by one; pad with trivial 1-blocks)."""
    out = Counter()
    for (e, l), n in q.items():
        out[(e, l + 1 if e == 1 else l)] += n
    fill = dim - rdim(out)
    assert fill >= 0, 'inconsistent local data'
    if fill:
        out[(1, 1)] += fill
    return out


def katz_MC(local, finite):
    """Katz RLS 3.3.6/3.3.7 for quadratic chi: input local data of a tame irreducible middle
    extension F (finite singular points `finite`, plus 'inf'); output local data and rank of MC_chi(F)."""
    m = sum(rdim(local[v]) - nblocks(local[v], 1) for v in finite)   # sum_s rank(F(s)/F(s)^I) = rank M(inf)
    d = m - nblocks(local['inf'], -1)                                   # minus rank((F(inf) tensor chi)^I)
    new = {v: reconstruct(twist(quot_inv(local[v])), d) for v in finite}
    M = reconstruct(local['inf'], m)
    new['inf'] = quot_inv(twist(M))
    assert rdim(new['inf']) == d
    return new, d


def initial_local(a):
    loc = {v: Counter({(1, 1): 1}) for v in range(a)}
    loc['y'] = Counter({(-1, 1): 1})
    loc['inf'] = Counter({(-1, 1): 1})
    return loc


def twist_by_label(local, T, a):
    g = dict(local)
    for v in T:
        g[v] = twist(local[v])
    if len(T) % 2:
        g['inf'] = twist(local['inf'])
    return g


def freeze(local, a):
    return tuple(tuple(sorted(local[v].items())) for v in list(range(a)) + ['y', 'inf'])


def t_fin(rep):
    return nblocks(rep, 1) - nblocks(rep, -1)


def t_inf(rep):
    return nblocks(rep, -1) - nblocks(rep, 1)


def recover_word(local, a):
    """parallel10 (1): peel labels by inverse MC and the signs of t_v."""
    finite = list(range(a)) + ['y']
    word = []
    while rdim(local['y']) > 1:
        g, e = katz_MC(local, finite)
        diffs = [t_fin(g[v]) for v in range(a)]
        assert all(diffs), 'zero difference: label not recoverable'
        T = tuple(v for v in range(a) if diffs[v] < 0)
        assert T, 'empty recovered label'
        local = twist_by_label(g, T, a)
        word.append(T)
    assert freeze(local, a) == freeze(initial_local(a), a)
    return tuple(reversed(word))


def word_ranks(a, word):
    """principal rank d_w and raw rank r_w (parallel10 (6)) of a word of anchor-index tuples."""
    finite = list(range(a)) + ['y']
    loc = initial_local(a)
    d = 1
    c = [0] * a
    r = 1
    for T in word:
        loc, d = katz_MC(twist_by_label(loc, T, a), finite)
        c = [r if v in T else c[v] for v in range(a)]
        r = 1 + sum(c)
    return d, r


def part_D():
    log('Part D: Katz local-monodromy bookkeeping, rank growth, word recovery, raw ranks')
    depth_by_a = {1: 14, 2: 7, 3: 5, 4: 4} if not FAST else {1: 10, 2: 5, 3: 4, 4: 3}
    results = {}
    for a, depth in depth_by_a.items():
        finite = list(range(a)) + ['y']
        labels = [tuple(v for v in range(a) if m >> v & 1) for m in range(1, 1 << a)]
        signatures = {}
        level_stats = {j: dict(words=0, min_rank=None, max_rank=0, min_jump=None, sum_r=0, sum_r2=0, sum_c2=0, sum_d2=0, sum_rd2=0) for j in range(0, depth + 1)}
        level_stats[0].update(words=1, min_rank=1, max_rank=1, sum_r=1, sum_r2=1, sum_d2=1)
        signatures[freeze(initial_local(a), a)] = ()

        def visit(local, d, delta, word, cvec, r):
            if len(word) == depth:
                return
            i = len(word) + 1
            for T in labels:
                g = twist_by_label(local, T, a)
                new, dn = katz_MC(g, finite)
                jump = dn - d
                # (R'): d_i >= i+1, d_i <= a d_{i-1} + 1, jumps nondecreasing and >= 1
                check('D.rank_at_least_i_plus_1', dn >= i + 1, dict(a=a, word=word + (T,), d=dn))
                check('D.rank_at_most_a_d_plus_1', dn <= a * d + 1, dict(a=a, word=word + (T,), d=dn))
                check('D.jump_nondecreasing_and_positive', jump >= max(delta, 1), dict(a=a, word=word + (T,), jump=jump, prev=delta))
                if i == 1:
                    u = len(T)
                    check('D.first_rank_formula', dn == u + 1 - (1 if u % 2 == 0 else 0), dict(a=a, T=T, d=dn))
                else:
                    # (R): d_{i+1} - d_i = -delta_i + sum_{v in F(T)} t_v(F_i)
                    pred = -delta + sum(t_fin(local[v]) for v in T) + (t_inf(local['inf']) if len(T) % 2 else 0)
                    check('D.formula_R', jump == pred, dict(a=a, word=word + (T,), jump=jump, pred=pred))
                # block-count formulas of the rank lemma (F = MC(G), delta = rank F - rank G)
                dl = dn - rdim(g['y'])
                for v in range(a):
                    check('D.finite_block_formula', nblocks(new[v], 1) == dl + nblocks(g[v], 1) and nblocks(new[v], -1) == nblocks(g[v], 1) - n1blocks(g[v], 1), dict(a=a, word=word + (T,), v=v))
                    check('D.t_v_at_least_delta', t_fin(new[v]) >= jump, dict(a=a, word=word + (T,), v=v))
                check('D.infinity_block_formula', nblocks(new['inf'], -1) == dl + nblocks(g['inf'], -1) and nblocks(new['inf'], 1) == nblocks(g['inf'], -1) - n1blocks(g['inf'], -1), dict(a=a, word=word + (T,)))
                check('D.t_inf_at_least_delta', t_inf(new['inf']) >= jump, dict(a=a, word=word + (T,)))
                # pseudoreflection at y
                expect_y = Counter({(-1, 1): 1, (1, 1): dn - 1}) if i % 2 == 0 else Counter({(1, 2): 1, (1, 1): dn - 2})
                expect_y = +expect_y
                check('D.pseudoreflection_at_y', +new['y'] == expect_y, dict(a=a, word=word + (T,), got=sorted(new['y'].items())))
                # involution: MC(MC(G)) = G with rank d
                back, dback = katz_MC(new, finite)
                check('D.involution', dback == rdim(g['y']) and all(+back[v] == +g[v] for v in finite + ['inf']), dict(a=a, word=word + (T,)))
                # word recovery and distinct signatures
                key = freeze(new, a)
                check('D.distinct_signature', key not in signatures, dict(a=a, word=word + (T,), clash=signatures.get(key)))
                signatures[key] = word + (T,)
                check('D.word_recovery', recover_word(new, a) == word + (T,), dict(a=a, word=word + (T,)))
                # raw ranks (6)
                cn = [r if v in T else cvec[v] for v in range(a)]
                rn = 1 + sum(cn)
                check('D.raw_rank_at_least_principal', rn >= dn, dict(a=a, word=word + (T,), r=rn, d=dn))
                st = level_stats[i]
                st['words'] += 1
                st['min_rank'] = dn if st['min_rank'] is None else min(st['min_rank'], dn)
                st['max_rank'] = max(st['max_rank'], dn)
                st['min_jump'] = jump if st['min_jump'] is None else min(st['min_jump'], jump)
                st['sum_r'] += rn
                st['sum_r2'] += rn * rn
                st['sum_c2'] += sum(x * x for x in cn)
                st['sum_d2'] += dn * dn
                st['sum_rd2'] += (rn - dn) ** 2
                visit(new, dn, jump, word + (T,), cn, rn)

        visit(initial_local(a), 1, 0, (), [0] * a, 1)
        # (7): closed recurrence for W, R, Q, Z
        rec_ok = True
        if a >= 2:
            N, h = 2 ** a - 1, 2 ** (a - 1)
            A2 = 2 ** (a - 2) * (a * a + 3 * a - 1) - 1
            W, R, Q, Z = 1, 1, 1, 0
            for m in range(depth + 1):
                st = level_stats[m]
                ok = (st['words'], st['sum_r'], st['sum_r2'], st['sum_c2']) == (W, R, Q, Z)
                rec_ok &= ok
                check('D.recurrence_7', ok, dict(a=a, m=m, got=(st['words'], st['sum_r'], st['sum_r2'], st['sum_c2']), rec=(W, R, Q, Z)))
                W, R, Q, Z = N * W, (h * (a + 1) - 1) * R + h * W, A2 * Q + h * (a + 2) * R + (h * Z) // 2 + (h * W) // 2, a * h * Q + (h - 1) * Z
            lam = 2 ** (a - 3) * (a * a + 3 * a + 1 + math.sqrt((a * a + 3 * a - 3) ** 2 + 8 * a)) - 1
            perron = max(np.linalg.eigvals(np.array([[A2, h / 2], [a * h, h - 1]], dtype=float)).real)
            check('D.Lambda_a_is_Perron_eigenvalue', abs(lam - perron) < 1e-9, dict(a=a, lam=lam, perron=perron))
            ratios = [level_stats[m + 1]['sum_r2'] / level_stats[m]['sum_r2'] for m in range(depth)]
        else:
            lam, ratios = None, None
            for m in range(depth + 1):
                check('D.a1_ranks', level_stats[m]['sum_r'] == m + 1 and level_stats[m]['sum_d2'] == (m + 1) ** 2, dict(m=m))
        results[f'a={a}'] = dict(depth=depth, levels=level_stats, distinct_signatures=len(signatures), Lambda_a=lam, Q_ratios=ratios)
        log(f'  a={a} depth={depth}: {len(signatures)} signatures')
    return results


# ---------------------------------------------------------------- Part E: Gram matrices of raw kernels
def all_words(a, m):
    labels = [tuple(v for v in range(a) if mask >> v & 1) for mask in range(1, 1 << a)]
    return list(itertools.product(labels, repeat=m))


def raw_paths(conv, p, chi, anchors, words_by_len, y):
    """R_w(., y) for all words in words_by_len (dict length -> list of words), sharing prefixes."""
    x = np.arange(p)
    paths = {(): chi[(x - y) % p].copy()}
    out = {}
    maxlen = max(words_by_len)
    frontier = {(): paths[()]}
    for step in range(1, maxlen + 1):
        nxt = {}
        for w, v in frontier.items():
            for T in set(u[step - 1] for L in words_by_len if L >= step for u in words_by_len[L] if u[:step - 1] == w):
                d = mask_vec(p, chi, tuple(anchors[i] for i in T))
                masked = d * v if v.dtype != object else np.array([int(a_) * int(b_) for a_, b_ in zip(v, d)], dtype=object)
                nxt[w + (T,)] = conv.conv(masked)
        frontier = nxt
        for w, v in nxt.items():
            if len(w) in words_by_len and w in words_by_len[len(w)]:
                out[w] = v
    return out


def gram_case(p, chi, conv, a, anchors, words_by_len, y, label):
    paths = raw_paths(conv, p, chi, anchors, words_by_len, y)
    U = np.array([t for t in range(p) if t not in set(anchors) | {y}])
    words = [w for L in sorted(words_by_len) for w in words_by_len[L]]
    Q = np.empty((len(words), len(U)))
    D_W = L_W = 0
    ranks = {}
    for i, w in enumerate(words):
        v = paths[w]
        if v.dtype == object:
            v = np.array([float(int(t)) for t in v])
        Q[i] = ((-1) ** len(w)) * v[U].astype(np.float64) / p ** (len(w) / 2)
        d, r = word_ranks(a, w)
        ranks[str(w)] = (d, r)
        D_W += d * d
        L_W += (r - d) ** 2
    G = Q @ Q.T / p
    dev = float(np.linalg.norm(G - np.eye(len(words)), 2))
    eps = (a * D_W + 1) / math.sqrt(p)
    gam = math.sqrt(L_W / p)
    bound = eps + 2 * gam * math.sqrt(1 + eps) + gam * gam
    off = np.abs(G - np.diag(np.diag(G)))
    diag = np.diag(G)
    # per-entry claimed constants: |G_uv - delta| <= (a d_u d_v + [u=v])/sqrt p + (r_u-d_u)(r_v-d_v)/p + cross terms;
    # we record the crude entrywise ratio |G_uv - delta_uv| sqrt(p) / (a d_u d_v + [u=v]) as a diagnostic only.
    ent = []
    for i, u in enumerate(words):
        for j, v in enumerate(words):
            du, dv = ranks[str(u)][0], ranks[str(v)][0]
            ent.append(abs(G[i, j] - (1 if i == j else 0)) * math.sqrt(p) / (a * du * dv + (1 if i == j else 0)))
    nonvac = bound < 1
    ok = dev <= bound * (1 + 1e-9)
    check('E.gram_bound_10' if nonvac else 'E.gram_bound_10_vacuous', ok, dict(label=label, p=p, dev=dev, bound=bound))
    log(f'  Gram {label}: p={p} words={len(words)} D_W={D_W} L_W={L_W} ||G-I||={dev:.4f} bound={bound:.4f} '
        f'nonvacuous={nonvac} dev*sqrt(p)={dev*math.sqrt(p):.2f} diag=[{diag.min():.4f},{diag.max():.4f}] maxoff={off.max():.5f} maxoff*sqrt(p)={off.max()*math.sqrt(p):.2f}')
    return dict(label=label, p=p, a=a, anchors=anchors, y=y, n_words=len(words), D_W=D_W, L_W=L_W,
                op_norm_dev=dev, claimed_bound=bound, nonvacuous=nonvac, satisfied=ok, dev_times_sqrt_p=dev * math.sqrt(p),
                diag_min=float(diag.min()), diag_max=float(diag.max()), max_offdiag=float(off.max()),
                max_offdiag_times_sqrt_p=float(off.max() * math.sqrt(p)),
                max_entrywise_ratio_to_a_dudv=max(ent), ranks=ranks)


def part_E():
    log('Part E: raw-kernel Gram matrices')
    out = []
    configs = [
        (10009, 1, [0], {L: all_words(1, L) for L in range(1, 7)}, 'a1_len1to6'),
        (next_prime_1mod4(100000), 1, [0], {L: all_words(1, L) for L in range(1, 7)}, 'a1_len1to6'),
        (next_prime_1mod4(100000), 2, [0, 1], {2: all_words(2, 2)}, 'a2_m2'),
        (next_prime_1mod4(100000), 2, [0, 1], {3: all_words(2, 3)}, 'a2_m3'),
        (next_prime_1mod4(100000), 2, [0, 1], {4: all_words(2, 4)}, 'a2_m4'),
        (next_prime_1mod4(100000), 2, [0, 1], {1: all_words(2, 1), 2: all_words(2, 2), 3: all_words(2, 3)}, 'a2_len1to3'),
        (next_prime_1mod4(100000), 3, [0, 1, 2], {2: all_words(3, 2)}, 'a3_m2'),
        (next_prime_1mod4(100000), 3, [0, 1, 2], {3: all_words(3, 3)}, 'a3_m3'),
    ]
    if not FAST:
        big = next_prime_1mod4(1000000)
        configs += [
            (big, 1, [0], {L: all_words(1, L) for L in range(1, 7)}, 'a1_len1to6'),
            (big, 2, [0, 1], {2: all_words(2, 2)}, 'a2_m2'),
            (big, 2, [0, 1], {3: all_words(2, 3)}, 'a2_m3'),
            (big, 3, [0, 1, 2], {2: all_words(3, 2)}, 'a3_m2'),
        ]
    cache = {}
    for p, a, anchors, wbl, label in configs:
        if p not in cache:
            chi = chi_table(p)
            cache[p] = (chi, Convolver(p, chi))
        chi, conv = cache[p]
        y = 5
        out.append(gram_case(p, chi, conv, a, anchors, wbl, y, label))
    return out


# ---------------------------------------------------------------- Part F: previous worker's data
def part_F():
    out = {}
    try:
        d2 = json.load(open(os.path.join(SCRATCH, 'explore2_a1.json')))['a1']
        per_k = {}
        for k in range(3, 9):
            pts = [(r['p'], abs(r['N'][str(k)])) for r in d2 if str(k) in r['N'] and r['N'][str(k)] != 0]
            xs = [math.log(q) for q, _ in pts]
            ys = [math.log(v) for _, v in pts]
            per_k[k] = dict(n=len(pts), max_p=max(q for q, _ in pts), slope=slope(xs, ys),
                            max_ratio_p_half_k_plus_1=max(v / q ** ((k + 1) / 2) for q, v in pts),
                            max_ratio_to_bound_T=max(v / (3 * 4 ** (k - 2) * q ** ((k + 1) / 2)) for q, v in pts),
                            max_ratio_to_kunisky_k_p=max(v / (k * q ** ((k + 1) / 2)) for q, v in pts))
        out['explore2_a1_homogeneous_up_to_2e6'] = per_k
    except Exception as e:
        out['explore2_a1_homogeneous_up_to_2e6'] = f'unavailable: {e}'
    try:
        d3 = json.load(open(os.path.join(SCRATCH, 'explore3.json')))
        summ = {}
        for name in sorted(set(r['name'] for r in d3)):
            for k in range(3, 9):
                pts = [(r['p'], max(abs(v) for v in r['I'])) for r in d3 if r['name'] == name and r['k'] == k]
                pts = [(q, v) for q, v in pts if v]
                xs = [math.log(q) for q, _ in pts]
                ys = [math.log(v) for _, v in pts]
                summ[f'{name}:k={k}'] = dict(n=len(pts), max_p=max(q for q, _ in pts), slope=slope(xs, ys), claimed=(k - 1) / 2,
                                            max_ratio=max(v / q ** ((k - 1) / 2) for q, v in pts))
        out['explore3_inner_sums_a2_a3_up_to_5e5'] = summ
    except Exception as e:
        out['explore3_inner_sums_a2_a3_up_to_5e5'] = f'unavailable: {e}'
    try:
        d4 = json.load(open(os.path.join(SCRATCH, 'explore4.json')))
        out['explore4_dense_1601_to_3413'] = dict(entries=len(d4), primes=sorted(set(r['p'] for r in d4)),
                                                  violations=[r for r in d4 if r['violations']])
    except Exception as e:
        out['explore4_dense_1601_to_3413'] = f'unavailable: {e}'
    return out


def main():
    res = dict(status='running')
    res['A_dense_traces'] = part_A()
    res['D_katz_bookkeeping'] = part_D()
    res['C_indicator_expansion'] = part_C()
    res['B_inner_sums'] = part_B()
    res['E_gram'] = part_E()
    res['F_previous_worker_data'] = part_F()
    res['counts'] = dict(COUNTS)
    res['total_checks'] = sum(COUNTS.values())
    res['witnesses'] = WITNESSES
    res['runtime_seconds'] = time.time() - T0
    res['status'] = ('no violation of any tested inequality' if not WITNESSES else f'{len(WITNESSES)} violations recorded')
    res['fast_mode'] = FAST
    json.dump(res, open(OUT, 'w'), indent=1, default=str)
    log('written', OUT, 'checks', res['total_checks'], 'witnesses', len(WITNESSES), 'status:', res['status'])


if __name__ == '__main__':
    main()
