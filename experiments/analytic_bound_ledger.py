#!/usr/bin/env python3
"""Exact algebra for analytic-bounds-and-amplification.md.

No numerical subgroup bound or proof of Paley is produced. The generic
amplification inequalities are hypotheses; their realizability is not tested.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def eliminate(r, s, ell, er, es, e2ell):
    # Multiply |X|, |Y|, |Z|^(1/2), Delta3^4.
    power_n = 2*r + 2*s + 2*ell - er - es - e2ell/2
    powers = [F(2*r), F(2*s-2), F(2*ell-2), F(3)]
    assert all(x >= 0 for x in powers)
    # Delta3 >= Delta2^(2*ell), Delta2 >= Delta1^s,
    # Delta1 >= Delta^r. These substitutions are in the lower bound.
    for j, multiplier in [(3, 2*ell), (2, s), (1, r)]:
        powers[j-1] += multiplier*powers[j]
        powers[j] = F(0)
    assert powers == [F(8*ell*r*s), 0, 0, 0]
    return power_n, powers[0], (power_n-4)/powers[0]


def energy_floor(j):
    return F(max(j, 2*j-4))


def ideal_saving(r, s, ell):
    _, _, saving = eliminate(r, s, ell, energy_floor(r),
                            energy_floor(s), energy_floor(2*ell))
    independent_formula = F(min(r, 4)+min(s, 4)+min(ell, 2)-4,
                            8*ell*r*s)
    assert saving == independent_formula
    return saving


def amplification_audit():
    a, b, saving = eliminate(3, 3, 1, F(4), F(4), F(49, 20))
    assert (a, b, saving) == (F(191, 40), 72, F(31, 2880))
    assert F(2689, 2880)+F(4, 72) == 1-saving
    ideal_source = eliminate(3, 3, 1, F(3), F(3), F(2))
    assert ideal_source == (7, 72, F(1, 24))

    # This is an exhaustive proof after the clamping lemma in the note:
    # positive maxima have r,s<=4 and ell<=2. It is not a guessed cutoff.
    reduced = [(r, s, ell, ideal_saving(r, s, ell))
               for r, s, ell in product(range(1, 5), range(1, 5), range(1, 3))]
    best = max(row[3] for row in reduced)
    winners = [list(row[:3]) for row in reduced if row[3] == best]
    assert best == F(1, 16) and winners == [[1, 4, 1], [4, 1, 1]]
    restricted = [row for row in reduced if min(row[:2]) >= 2]
    restricted_best = max(row[3] for row in restricted)
    restricted_winners = [list(row[:3]) for row in restricted
                          if row[3] == restricted_best]
    assert restricted_best == F(3, 64)
    assert restricted_winners == [[2, 4, 1], [4, 2, 1]]
    # Sanity checks of the clamp on larger inputs supplement its proof.
    checks = 0
    for r, s, ell in product([1, 2, 3, 4, 5, 17, 101], repeat=3):
        current = ideal_saving(r, s, ell)
        if current > 0:
            assert current <= ideal_saving(min(r, 4), min(s, 4), min(ell, 2))
        checks += 1
    return {"published_source_case": {"power_n": str(a), "power_delta": str(b),
                                       "quartic_saving": str(saving),
                                       "amplitude_exponent": str(1-saving)},
            "ideal_source_case_saving": str(ideal_source[2]),
            "abstract_ledger_best_saving": str(best), "attainers": winners,
            "r_and_s_at_least_two_best_saving": str(restricted_best),
            "restricted_attainers": restricted_winners,
            "all_reduced_cases": [{"r": r, "s": s, "ell": ell, "saving": str(value)}
                                  for r, s, ell, value in reduced],
            "clamp_sanity_checks": checks,
            "scope": "Conditional exponent ledger only; not a general impossibility theorem."}


def centered_formula_audit(e2):
    # Conditional on the displayed subgroup moment formula in the note.
    def saving(k):
        return (F(k+9, 2)-e2-4) / 2**(k+1)

    rows = [{"k": k, "r": 2**k, "saving": str(saving(k))} for k in range(1, 12)]
    # Let a=2*e2-1. Then B(k)=(k-a)/(4*2^k), and
    # B(k+1)-B(k)=(a+1-k)/(8*2^k). This proves the global maximum.
    a = 2*e2-1
    for k in range(1, 12):
        assert saving(k+1)-saving(k) == (a+1-k)/(8*2**k)
    best = max(saving(k) for k in range(1, 12))
    attainers = [k for k in range(1, 12) if saving(k) == best]
    assert a+1 < 11  # Beyond the checked candidates it strictly decreases.
    if e2 == F(49, 20):
        assert best == F(11, 1280) and attainers == [5]
    else:
        assert e2 == 2 and best == F(1, 64) and attainers == [4, 5]
    return {"energy_two_exponent": str(e2), "best_saving": str(best),
            "attaining_k": attainers, "fixed_k_rows": rows,
            "scope": "Algebraic implication of a published formula; not its independent proof."}


def convolve(left, right):
    p = len(left)
    return [sum(left[x]*right[(z-x) % p] for x in range(p)) for z in range(p)]


def origin_test():
    # f=delta_0-1/p is a mean-zero convolution idempotent.
    # Exact rational checks are independent of the Fourier derivation.
    examples = []
    for p in [3, 5, 17, 97]:
        f = [F(p-1, p)] + [F(-1, p)]*(p-1)
        assert sum(f) == 0 and convolve(f, f) == f
        assert all(f[x] == f[(x*g) % p] for x in range(p) for g in range(1, p))
        norm1 = sum(abs(x) for x in f)
        norm2_squared = sum(x*x for x in f)
        assert norm1 == 2*F(p-1, p) and norm2_squared == F(p-1, p)
        invalid_bound = norm1**2/(p-1)
        if p > 4:
            assert norm2_squared > invalid_bound
        examples.append({"p": p, "norm1": str(norm1),
                         "all_positive_order_T_values": str(norm2_squared),
                         "norm_bound_rhs_in_equation_57": str(invalid_bound),
                         "equation_57_fails": norm2_squared > invalid_bound})
    assert 2**12*2**4 == 65536
    return {"examples": examples,
            "k_two_required_inequality": "1 <= 65536*C*(log p)^4*(1-1/p)^4/sqrt(p-1)",
            "right_side_limit_for_fixed_C": "0 as p tends to infinity",
            "scope": "Unqualified arbitrary-function statement as rendered in v2 HTML; not subgroup Theorem 3."}


def archive_audit():
    archive = ROOT / "sources/analytic-bounds-2026-09-04"
    manifest = json.loads((archive/"manifest.json").read_text())
    for row in manifest["files"]:
        contents = (archive/row["file"]).read_bytes()
        assert len(contents) == row["bytes"]
        assert sha256(contents).hexdigest() == row["sha256"]
    return {"manifest_sha256": sha256((archive/"manifest.json").read_bytes()).hexdigest(),
            "verified_files": [row["file"] for row in manifest["files"]]}


def main():
    result = {"status": "passed; Paley and prize unproved",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "amplification": amplification_audit(),
              "conditional_centered_formula": [centered_formula_audit(F(49, 20)),
                                               centered_formula_audit(F(2))],
              "origin_mass_check": origin_test(), "archive": archive_audit()}
    destination = ROOT / "results/analytic_bound_ledger.json"
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"],
                      "published_quartic_saving": "31/2880",
                      "abstract_ledger_maximum_saving": "1/16",
                      "conditional_centered_best_saving": "11/1280",
                      "source_files_verified": len(result["archive"]["verified_files"]),
                      "output": str(destination)}, indent=2))


if __name__ == "__main__":
    main()
