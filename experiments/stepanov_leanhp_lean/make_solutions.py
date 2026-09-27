#!/usr/bin/env python3
"""Generate prove2me-format statement and solution files from the development files.

Theorems/Thm_<name>.lean : preamble + `theorem <name> ... := by sorry`   (what /submit-problem gets)
Solutions/Sol_<name>.lean: self-contained: imports + the development namespace + `theorem solution`
                           with exactly the same binders and type, closed by the development theorem.
Run from anywhere; writes next to this script.  Standard library only.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def body(fname, ns):
    src = open(os.path.join(HERE, fname)).read()
    imports = [ln for ln in src.splitlines() if ln.startswith("import ")]
    opens = [ln for ln in src.splitlines() if ln.startswith("open ")]
    m = re.search(r"(namespace %s\n.*?\nend %s\n)" % (ns, ns), src, re.S)
    assert m, (fname, ns)
    return imports, "section\n" + "\n".join(opens) + "\n\n" + m.group(1) + "\nend\n"


PREAMBLE = """import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

open Finset"""

PREAMBLE_REAL = """import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.Data.Real.Sqrt

set_option autoImplicit false

open Finset"""

# name -> (preamble, binders+type (text after the name, before ':='), dev file, namespace, proof term)
TARGETS = {
    "hanson_petridis_paley": (PREAMBLE, """(p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (A B : Finset (ZMod p))
    (hAB : ∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) :
    A.card * B.card ≤ (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card""",
        "StepanovHP.lean", "StepanovHP", "StepanovHP.hanson_petridis_paley p hp A B hAB"),
    "hanson_petridis_theorem_1_2": (PREAMBLE, """(p : ℕ) [Fact p.Prime] (d : ℕ) (hd : d ∣ p - 1) (hd' : d < p - 1)
    (A B : Finset (ZMod p)) (hAB : ∀ a ∈ A, ∀ b ∈ B, a + b = 0 ∨ (a + b) ^ d = 1) :
    A.card * B.card ≤ d + (B.filter (fun b => -b ∈ A)).card""",
        "StepanovHP.lean", "StepanovHP", "StepanovHP.hanson_petridis p d hd hd' A B hAB"),
    "paley_clique_number_bound": (PREAMBLE, """(p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    A.card * (A.card - 1) ≤ (p - 1) / 2""",
        "StepanovHP.lean", "StepanovHP", "StepanovHP.paley_clique_number_bound p hp A hA"),
    "paley_hp_sharp_example": (PREAMBLE, """(p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) :
    let A : Finset (ZMod p) := {0, 1}
    let B : Finset (ZMod p) := Finset.univ.filter (fun b => IsSquare b ∧ IsSquare (b + 1))
    B.card = (p + 3) / 4 ∧ (∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) ∧
      A.card * B.card = (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card""",
        "StepanovSharp.lean", "StepanovHPSharp", "StepanovHPSharp.paley_hp_sharp_example p hp"),
    # statements already on prove2me (private, CITED/Open), verbatim from research/sigma-lean-2026-09-05.md
    "paley_hanson_petridis_difference_bound": (PREAMBLE, """(p : ℕ) [Fact p.Prime] (d : ℕ)
    (hd : d ∣ p - 1) (hd' : d < p - 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a - a' = 0 ∨ (a - a') ^ d = 1) :
    A.card * (A.card - 1) ≤ d""",
        "StepanovHP.lean", "StepanovHP", "StepanovHP.difference_bound p d hd hd' A hA"),
    "paley_hanson_petridis_clique_number": (PREAMBLE_REAL, """(p : ℕ) [Fact p.Prime] (hp : p % 4 = 1)
    (A : Finset (ZMod p)) (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    (A.card : ℝ) ≤ (Real.sqrt (2 * (p : ℝ) - 1) + 1) / 2""",
        "StepanovHP.lean", "StepanovHP", "StepanovHP.clique_number_real p hp A hA"),
}


def main():
    os.makedirs(os.path.join(HERE, "Theorems"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "Solutions"), exist_ok=True)
    for name, (pre, sig, dev, ns, proof) in TARGETS.items():
        stmt = f"theorem {name} {sig} := by sorry\n"
        with open(os.path.join(HERE, "Theorems", f"Thm_{name}.lean"), "w") as fh:
            fh.write(pre + "\n\n" + stmt)
        imports, nsbody = body(dev, ns)
        pre_imports = [ln for ln in pre.splitlines() if ln.startswith("import ")]
        all_imports = []
        for ln in imports + pre_imports:
            if ln not in all_imports:
                all_imports.append(ln)
        sol = ("\n".join(all_imports) + "\n\nset_option autoImplicit false\n\n"
               + nsbody + "\nopen Finset\n\n"
               + f"theorem solution {sig} :=\n  {proof}\n\n#print axioms solution\n")
        with open(os.path.join(HERE, "Solutions", f"Sol_{name}.lean"), "w") as fh:
            fh.write(sol)
        print("wrote", name)


if __name__ == "__main__":
    main()
