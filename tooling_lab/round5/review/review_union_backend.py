#!/usr/bin/env python3
"""Independent column-by-column oracle for the Newton union backend.

This intentionally uses neither Newton identities, grouped sign power sums,
nor the backend's cover coefficients. It expands one literal column at a time.
"""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import random
import subprocess
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'third_moment'/'union_coefficients.cpp'


def literal_columns(columns,degree,max_union):
    state={(0,0,0,0):1}
    for col in columns:
        terms=[]
        for mask in range(1,8):
            exp=tuple((mask>>i)&1 for i in range(3))
            value=1
            for i in range(3):
                if exp[i]:value*=col[i]
            if value:terms.append((exp,value))
        nxt=state.copy()
        for (k,a,b,c),old in state.items():
            if k==max_union:continue
            for (x,y,z),value in terms:
                if a+x<=degree and b+y<=degree and c+z<=degree:
                    key=k+1,a+x,b+y,c+z
                    nxt[key]=nxt.get(key,0)+old*value
        state={key:value for key,value in nxt.items() if value}
    return [state.get((k,degree,degree,degree),0) for k in range(max_union+1)]


def main():
    started=time.perf_counter()
    source_hash=sha256(SOURCE.read_bytes()).hexdigest()
    binary=HERE/'review_union_backend'
    subprocess.run(['clang++','-std=c++17','-O2','-I/opt/homebrew/include',str(SOURCE),'-o',str(binary)],check=True)
    cases=[]
    for col in product((-1,0,1),repeat=3):
        for degree in (0,1):cases.append(([col],degree,min(1,3*degree)))
    rng=random.Random(95013)
    for q in (3,5,7,8):
        for degree in range(7):
            for trial in range(3):
                columns=[tuple(rng.choice((-1,0,1)) for _ in range(3)) for _ in range(q)]
                cases.append((columns,degree,min(q,3*degree)))
    # Nontrivial full-degree elementary products, zero columns, and repeated rows.
    for q in (6,7,8):
        for template in ([(1,1,1)],[(1,-1,1),(-1,1,-1)],[(1,1,-1),(0,0,0),(1,-1,1)]):
            columns=[template[i%len(template)] for i in range(q)]
            for degree in (3,6):cases.append((columns,degree,min(q,3*degree)))
    payload=[str(len(cases))]
    expected=[]
    for columns,d,k in cases:
        hist=Counter(columns)
        payload.append(f'{d} {k} {len(hist)}')
        payload.extend(' '.join(map(str,(*col,count))) for col,count in sorted(hist.items()))
        expected.append(literal_columns(columns,d,k))
    raw=subprocess.check_output([str(binary)],input='\n'.join(payload)+'\n',text=True)
    results=[json.loads(line) for line in raw.splitlines()]
    assert len(results)==len(cases)
    for i,(row,want) in enumerate(zip(results,expected)):
        assert row['case_id']==i and row['coefficients']==want,(i,cases[i],row,want)
    assert sha256(SOURCE.read_bytes()).hexdigest()==source_hash,'source changed during review'
    out={'status':'all literal column expansions agree exactly','case_count':len(cases),
         'degrees':list(range(7)),'populations':[1,3,5,6,7,8],
         'coefficient_equalities':sum(len(x) for x in expected),
         'independence':'column-by-column multiplication with seven possible nonempty row choices; no Newton identities or grouped cover coefficients',
         'source_sha256':{'../third_moment/union_coefficients.cpp':source_hash,
                          Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()},
         'elapsed_seconds':time.perf_counter()-started}
    (HERE/'review_union_backend.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
