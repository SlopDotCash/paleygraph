#!/usr/bin/env python3
"""Find actual same-cell-profile inputs before comparing conditional moments."""
import hashlib,itertools,json,sys,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
ROUND=HERE.parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROUND/'marked_moments'))


def profile(S,M):
    return tuple(sorted(Counter(tuple(int(S[x][m]) for m in M) for x in range(len(S))).items()))


def canonical_profile(S,M):
    options=[(profile(S,P),P) for P in itertools.permutations(M)]
    return min(options)


def prime_matrix(p):
    chi=[0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1) for x in range(p)]
    return [[chi[(x-y)%p] for y in range(p)] for x in range(p)]


def scan_saved_oracle():
    path=ROUND/'novelty/direct_conditional_oracle.json'
    data=json.loads(path.read_text());findings=[];censuses=[]
    for case in data['cases']:
        if case['n']<=6 or case['q']>17:continue
        S=prime_matrix(case['q']);groups={};canonical={};candidates=0
        for r in case['records']:
            M=tuple(r['marks'])
            if len(M)!=3:continue
            k=profile(S,M)
            if k in groups and groups[k]['second_moment']!=r['second_moment']:
                findings.append({'graph':case['graph'],'q':case['q'],'n':case['n'],'left':groups[k],'right':r,'matching':'literal ordered profile','profile':k})
                break
            groups[k]=r;candidates+=1
        censuses.append({'graph':case['graph'],'n':case['n'],'marks':3,'literal_profiles':len(groups),'checked_mark_sets':candidates,'difference_found':bool(findings and findings[-1]['graph']==case['graph'])})
    out={'oracle_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'findings':findings,'censuses':censuses}
    (HERE/'saved_oracle_scan.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':scan_saved_oracle()
