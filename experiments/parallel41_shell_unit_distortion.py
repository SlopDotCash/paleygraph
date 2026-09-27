#!/usr/bin/env python3
"""Check three specified unit multiples and an exact distortion threshold.

This is a deterministic certificate check, not a unit search or optimizer.
"""

from pathlib import Path
import importlib.util
import json


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "real_recovery", ROOT / "experiments/parallel41_shell_real_recovery.py"
)
arithmetic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(arithmetic)


def main():
    source = json.loads((ROOT / "results/parallel41_shell_real_recovery_2026_09_06.json").read_text())
    p, h, f = source["p"], source["h"], source["recovered_f"]
    N = len(h)
    unit = [0] * N
    unit[0], unit[1], unit[-1] = 1, 1, -1
    assert arithmetic.conjugate(unit) == unit
    assert arithmetic.norm_adjugate(unit)[0] == 1
    rows = []
    for power in (0, 1, 2):
        assert arithmetic.conjugate(h) == h
        assert arithmetic.norm_adjugate(h)[0] == p * p
        assert arithmetic.norm_adjugate(f)[0] == 1217 * p
        rows.append({
            "unit_power": power,
            "h": h,
            "f": f,
            "h_coefficient_energy": sum(x * x for x in h),
            "f_coefficient_energy": sum(x * x for x in f),
        })
        if power < 2:
            h = arithmetic.multiply(h, unit)
            f = arithmetic.multiply(f, unit)
    assert [row["h_coefficient_energy"] for row in rows] == [19, 53, 305]
    assert [row["f_coefficient_energy"] for row in rows] == [13, 43, 293]
    margin = 5**16 - 2**37
    assert margin == 15148937153 > 0
    result = {
        "scope": "Exact checks for three specified real-unit multiples and the integer threshold; the universal trace obstruction is proved in the note.",
        "p": p, "n": 2 * N, "N": N,
        "real_unit": unit,
        "complex_norm_unit": 1,
        "complex_norm_h_for_every_row": p * p,
        "complex_norm_f_for_every_row": 1217 * p,
        "rows": rows,
        "threshold_5_to_16_minus_2_to_37": margin,
        "no_optimization_claim": True,
    }
    path = ROOT / "results/parallel41_shell_unit_distortion_2026_09_06.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "h_energies": [19, 53, 305],
        "f_energies": [13, 43, 293],
        "threshold_margin": margin,
        "certificate": str(path),
    }, indent=2))


if __name__ == "__main__":
    main()
