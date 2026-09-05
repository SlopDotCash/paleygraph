#!/usr/bin/env python3
"""Compile the algebraic core and cross-check sparse counts by direct enumeration."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]
WORKSPACE=Path.home()/'prove2me_workspace'
def H(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def main():
    proof='prove2me/Check_paley_balanced_split_collision.lean'
    proof_hash=H(proof)
    assert not re.search(r'\b(sorry|admit|axiom|unsafe)\b|import Theorems\.',(ROOT/proof).read_text())
    started=time.monotonic()
    lean=subprocess.run([str(Path.home()/'.local/bin/lake'),'env','lean',str(ROOT/proof)],
        cwd=WORKSPACE,capture_output=True,text=True)
    log=lean.stdout+lean.stderr
    (ROOT/'results/parallel27_balanced_split_lean_2026_09_05.log').write_text(log)
    assert lean.returncode==0 and 'sorryAx' not in log
    assert H(proof)==proof_hash
    axioms={name:[a.strip() for a in values.split(',') if a.strip()]
        for name,values in re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",log)}
    assert len(axioms)==3
    assert all(set(v)<={'propext','Classical.choice','Quot.sound'} for v in axioms.values())
    toolchain=(WORKSPACE/'lean-toolchain').read_text().strip()
    assert toolchain=='leanprover/lean4:v4.30.0'
    mathlib=subprocess.run(['git','-C',str(WORKSPACE/'.lake/packages/mathlib'),'rev-parse','HEAD'],
        capture_output=True,text=True,check=True).stdout.strip()
    assert mathlib=='c5ea00351c28e24afc9f0f84379aa41082b1188f'
    print(json.dumps({'local_lean':'passed','statements':len(axioms)}),flush=True)

    sparse_path='results/parallel27_six_distinct_2026_09_05.json'
    sparse_hash=H(sparse_path)
    sparse=json.loads((ROOT/sparse_path).read_text())
    rows=sparse['cases']
    cpp='experiments/parallel27_direct_six.cpp'
    checked=[]
    with tempfile.TemporaryDirectory(prefix='paley-direct-six-') as temporary:
        executable=Path(temporary)/'direct-six'
        compiled=subprocess.run(['/usr/bin/clang++','-std=c++17','-O3','-Wall','-Wextra',
            str(ROOT/cpp),'-o',str(executable)],text=True,capture_output=True)
        assert compiled.returncode==0,compiled.stderr
        for expected in rows:
            case_start=time.monotonic()
            output=subprocess.run([str(executable)],input=f"{expected['p']} {expected['n']}\n",
                                  capture_output=True,text=True,check=True).stdout
            observed=json.loads(output)
            for k in ['p','n','E3','T6','J6','repeated_R6','distinct_balanced_R6','distinct_unbalanced_R6']:
                assert observed[k]==expected[k],(k,expected['p'],observed[k],expected[k])
            observed['seconds']=time.monotonic()-case_start
            checked.append(observed)
            print(json.dumps({'p':observed['p'],'n':observed['n'],'direct':'passed',
                              'seconds':observed['seconds']}),flush=True)
    assert sparse_hash==H(sparse_path)
    inputs=[proof,cpp,'experiments/parallel27_independent_check_2026_09_05.py',sparse_path,
            'results/parallel27_balanced_split_lean_2026_09_05.log']
    result={'status':'passed','checked_at_utc':datetime.now(timezone.utc).isoformat(),
        'seconds':time.monotonic()-started,'toolchain':toolchain,'mathlib_rev':mathlib,
        'axioms_by_theorem':axioms,'direct_cases':checked,
        'input_sha256':{p:H(p) for p in inputs},
        'scope':'Same-author independently implemented integer enumeration, plus Lean algebraic core. '
                'Full combinatorial counting formula is an ordinary proof; no uniform upper bound.'}
    (ROOT/'results/parallel27_independent_check_2026_09_05.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','direct_cases':len(checked),'lean_statements':len(axioms)}),flush=True)

if __name__=='__main__': main()
