#!/usr/bin/env python3
"""Check the local proof using an existing pinned mathlib dependency cache.

The cache is read only. All generated output stays in this workspace.
Pass --cache-workspace to use another Lean 4.29.1 project with the same mathlib.
"""
import argparse
import hashlib
import os
from pathlib import Path
import re
import subprocess

PIN = "5e932f97dd25535344f80f9dd8da3aab83df0fe6"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--cache-workspace", type=Path,
                    default=Path.home() / "TheLeaningOfEverything")
mode = parser.add_mutually_exclusive_group()
mode.add_argument("--subgroup", action="store_true",
                    help="Check the subgroup arithmetic certificate instead")
mode.add_argument("--quartic", action="store_true",
                  help="Check the quartic-window subgroup arithmetic certificate")
mode.add_argument("--sharp", action="store_true",
                  help="Check the analytic counterexample to the square-root-of-two constant")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
cache = args.cache_workspace / ".lake" / "packages"
mathlib = cache / "mathlib"
head = subprocess.check_output(["git", "-C", str(mathlib), "rev-parse", "HEAD"], text=True).strip()
if head != PIN:
    raise SystemExit(f"Expected mathlib {PIN}; found {head}")
dirty = subprocess.check_output(
    ["git", "-C", str(mathlib), "status", "--porcelain", "--untracked-files=no"], text=True)
if dirty.strip():
    raise SystemExit("The mathlib cache has tracked changes; use a clean pinned cache")
lean = Path.home() / ".elan/toolchains/leanprover--lean4---v4.29.1/bin/lean"
env = os.environ.copy()
env["LEAN_PATH"] = os.pathsep.join(str(p / ".lake/build/lib/lean")
                                   for p in sorted(cache.iterdir()) if p.is_dir())
version = subprocess.check_output([str(lean), "--version"], text=True).strip()
source_name = ("SharpConstantCounterexample.lean" if args.sharp else
               "SubgroupQuarticCounterexample.lean" if args.quartic else
               "SubgroupCounterexample.lean" if args.subgroup else "MomentObstruction.lean")
source = root / "research" / source_name
if re.search(r"\b(sorry|admit|axiom)\b", source.read_text()):
    raise SystemExit("Proof source contains a forbidden hole or axiom declaration token")
proc = subprocess.run([str(lean), str(source)],
                      env=env, capture_output=True, text=True, cwd=root)
allowed = {"propext", "Classical.choice", "Quot.sound"}
censuses = re.findall(r"depends on axioms: \[([^]]*)\]", proc.stdout)
audit_ok = len(censuses) == (6 if args.sharp else 4 if args.subgroup or args.quartic else 3) and all(
    {a.strip() for a in census.split(",") if a.strip()} <= allowed for census in censuses)
audit_ok = audit_ok and "sorryAx" not in proc.stdout + proc.stderr
digest = hashlib.sha256(source.read_bytes()).hexdigest()
log = (f"{version}\nmathlib HEAD: {head}\nsource SHA-256: {digest}\n"
       f"exit code: {proc.returncode}\naxiom audit passed: {audit_ok}\n"
       + proc.stdout + proc.stderr)
(root / "results").mkdir(exist_ok=True)
logname = ("sharp-constant-lean-verification.txt" if args.sharp else
           "quartic-subgroup-lean-verification.txt" if args.quartic else
           "subgroup-lean-verification.txt" if args.subgroup else "lean-verification.txt")
(root / "results" / logname).write_text(log)
print(log, end="")
raise SystemExit(proc.returncode or (0 if audit_ok else 1))
