"""Primitive-coset polynomial and balanced-energy checks with integer arithmetic."""

from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def primes_through(limit):
    sieve = bytearray(b"\1") * (limit + 1)
    sieve[:2] = b"\0\0"
    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i :: i] = b"\0" * ((limit - i * i) // i + 1)
    return [i for i in range(limit + 1) if sieve[i]]


PRIMES = primes_through(600000)


def isprime(p):
    assert math.isqrt(p) <= PRIMES[-1]
    if p < 2:
        return False
    for q in PRIMES:
        if q * q > p:
            return True
        if p % q == 0:
            return p == q
    raise AssertionError("Sieve too small")


def generator(p, n):
    assert (p - 1) % n == 0
    for a in range(2, p):
        g = pow(a, (p - 1) // n, p)
        if pow(g, n // 2, p) == p - 1:
            assert pow(g, n, p) == 1
            return g
    raise AssertionError("No generator")


def primitive_labels(p, n, g):
    h, step = g, g * g % p
    out = []
    for _ in range(n // 4):
        out.append(pow(1 + h, n, p))
        h = h * step % p
    return out


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def remainder(a, b, p):
    a = trim([x % p for x in a])
    b = trim([x % p for x in b])
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        shift = len(a) - len(b)
        c = a[-1] * inv % p
        for i, x in enumerate(b):
            a[i + shift] = (a[i + shift] - c * x) % p
        trim(a)
    return a


def gcd(a, b, p):
    a, b = trim([x % p for x in a]), trim([x % p for x in b])
    while b:
        a, b = b, remainder(a, b, p)
    if a:
        inv = pow(a[-1], -1, p)
        a = [x * inv % p for x in a]
    return a


def derivative(a):
    return [i * a[i] for i in range(1, len(a))]


def from_roots(roots, p):
    out = [1]
    for x in roots:
        new = [0] * (len(out) + 1)
        for i, c in enumerate(out):
            new[i] = (new[i] - x * c) % p
            new[i + 1] = (new[i + 1] + c) % p
        out = new
    return out


def primitive_polynomial(n):
    k, degree = n // 2, n // 4
    powers = [0]
    for r in range(1, degree + 1):
        powers.append(k // 2 * sum((-1) ** a * math.comb(n * r, k * a) for a in range(2 * r + 1)))
    elementary = [1]
    for j in range(1, degree + 1):
        numerator = sum((-1) ** (i - 1) * elementary[j - i] * powers[i] for i in range(1, j + 1))
        assert numerator % j == 0
        elementary.append(numerator // j)
    out = [(-1) ** (degree - i) * elementary[degree - i] for i in range(degree + 1)]
    assert out[0] == 2**k and out[-1] == 1
    return out


def profile(p, n, literal=False, integer_poly=None):
    assert isprime(p)
    g = generator(p, n)
    labels = primitive_labels(p, n, g)
    mult = Counter(labels)
    k = n // 2
    defects = [sum(max(c - j, 0) for c in mult.values()) for j in range(1, max(mult.values()) + 1)]
    energy = 2 * k * sum(c * c for c in mult.values())
    assert energy == k * k + 4 * k * sum(defects)
    d = defects[0]
    assert k * k + 4 * k * d <= energy <= k * k + 2 * k * d * (d + 1)
    poly = from_roots(labels, p)
    if integer_poly is not None:
        assert poly == [a % p for a in integer_poly]
    running = poly
    deriv = poly
    gcd_degrees = []
    for j in range(1, min(max(mult.values()) + 1, len(poly))):
        deriv = derivative(deriv)
        running = gcd(running, deriv, p)
        gcd_degrees.append(len(running) - 1)
        assert gcd_degrees[-1] == defects[j - 1]
    if max(mult.values()) <= 2:
        assert energy <= 2 * k * k
    out = {
        "p": p, "parent_order": n, "generator": g,
        "quartic_window": n**4 <= 4 * p <= 4 * n**4,
        "balanced_energy": energy, "baseline": k*k,
        "fiber_multiplicity_histogram": dict(sorted(Counter(mult.values()).items())),
        "squarefree_defect": d, "successive_gcd_degrees": gcd_degrees,
        "maximum_fiber": max(mult.values()),
    }
    if literal:
        K = [pow(g, 2 * j, p) for j in range(k)]
        L = [g * a % p for a in K]
        w = Counter((a + b) % p for a in K for b in L)
        assert sum(c * c for c in w.values()) == energy
        assert max(w.values()) == max(mult.values())
        assert Counter(w.values()) == Counter({c: n * count for c, count in Counter(mult.values()).items()})
        out["literal_pair_check"] = True
        aa = Counter((a + b) % p for a in K for b in K)
        hh = Counter((a + b) % p for a in K + L for b in K + L)
        child_energy = sum(c*c for c in aa.values())
        parent_energy = sum(c*c for c in hh.values())
        t = sum(c*w.get(x, 0) for x, c in aa.items())
        assert parent_energy == 2*child_energy + 6*energy + 8*t
        H = set(K + L)
        eps = int(3 in H)
        repeated = sum((-2-b) % p in H for b in H)
        assert (repeated - 1 - 2*eps) % 2 == 0
        u = (repeated - 1 - 2*eps)//2
        rest = parent_energy - (3*n*n-3*n) - n*(4*eps + 12*u)
        assert rest >= 0 and rest % (24*n) == 0
        out.update({"child_energy": child_energy, "parent_energy": parent_energy, "T": t,
                    "epsilon": eps, "repeated_entry_orbits": u, "four_distinct_orbits": rest//(24*n)})
    return out


def scan(n, count):
    p, checked, collisions, max_fiber, end = n**4 // 4 + 1, 0, [], 0, None
    while checked < count:
        if isprime(p):
            g = generator(p, n)
            mult = Counter(primitive_labels(p, n, g))
            maximum = max(mult.values())
            max_fiber = max(max_fiber, maximum)
            if maximum > 1:
                collisions.append(profile(p, n, literal=len(collisions) < 2))
            checked += 1
            end = p
        p += n
    return {"parent_order": n, "eligible_primes_checked": checked, "start": n**4 // 4 + 1,
            "last_prime": end, "maximum_fiber": max_fiber, "collisions": collisions}


def main():
    old_path = ROOT / "results/kernel_discriminant.json"
    old = json.loads(old_path.read_text())
    integers = {n: primitive_polynomial(n) for n in (4, 8, 16, 32, 64, 128)}
    historical = []
    for entry in old["small_orders"]:
        n = entry["n"]
        for row in entry["all_eligible_exceptional_primes"]:
            historical.append(profile(row["p"], n, literal=True, integer_poly=integers[n]))
    witnesses = [profile(6700417, 64, True, integers[64]), profile(67403009, 128, True, integers[128])]
    scans = []
    for n in (256, 512, 1024):
        result = scan(n, 2048)
        scans.append(result)
        print(json.dumps({k: v if k != "collisions" else len(v) for k, v in result.items()}), flush=True)
    script = ROOT / "experiments/parallel2_subgroup_2026_09_04.py"
    result = {
        "status": "Exact fiber and gcd criteria proved; triple-free quartic hypothesis is not proved by these scans.",
        "integer_primitive_polynomials": {str(n): [str(x) for x in poly] for n, poly in integers.items()},
        "historical_exception_profiles": historical,
        "quartic_witnesses": witnesses,
        "quartic_scans": scans,
        "source_sha256": {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in (script, old_path)},
    }
    (ROOT / "results/parallel2_subgroup_2026_09_04.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"historical_profiles": len(historical), "historical_maximum_fiber": max(x["maximum_fiber"] for x in historical),
                      "scanned_quartic_primes": sum(x["eligible_primes_checked"] for x in scans),
                      "new_balanced_collision_fields": sum(len(x["collisions"]) for x in scans)}), flush=True)


if __name__ == "__main__":
    main()
