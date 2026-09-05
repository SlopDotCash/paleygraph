#!/usr/bin/env python3
"""Search actual same-profile marked inputs; retain a concrete witness or census."""
from collections import Counter,defaultdict
import argparse,hashlib,itertools,json,sys,time
from pathlib import Path
import numpy as np
from search_profiles import HERE,LAB,ROUND,profile,prime_matrix
from conditional_moments import matrix_counts,conditional_moments,canonical_pair


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonicalize(S,M):
    P=profile(S,M);best=None
    for perm in itertools.permutations(range(len(M))):
        candidate=tuple(sorted((tuple(u[i] for i in perm),size) for u,size in P))
        if best is None or candidate<best[0]:best=(candidate,tuple(M[i] for i in perm))
    return best


def rel_signature(record):
    c=Counter()
    for r in record['relations']:c[canonical_pair(r['left_pattern'],r['right_pattern'],r['sign'])]+=r['count']
    return tuple(sorted(c.items()))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--twins-only',action='store_true');parser.add_argument('--cross-twins',action='store_true');args=parser.parse_args()
    if args.cross_twins:args.twins_only=True
    prefix='cross_twins_' if args.cross_twins else ('twins_' if args.twins_only else '')
    started=time.monotonic();twin_path=LAB/'round3/pair_type_twins/results.json';twins=json.loads(twin_path.read_text())['sign_matrices']
    graphs=[] if args.twins_only else [(f'Paley{p}',prime_matrix(p),[(0,1,t) for t in range(2,p)]) for p in [29,37,41,53,61,73,89,101]]
    graphs.extend((name,S,[(0,a,b) for a,b in itertools.combinations(range(1,49),2)]) for name,S in zip(['Paley49','Peisert49'],twins))
    # These are actual graph records, not freely perturbed relation tables.
    pools={};signatures=defaultdict(set);census=[]
    for name,S,marks in graphs:
        S=np.array(S,dtype=np.int64);q=len(S);matrix_counts(S,[]);stats=Counter()
        for M in marks:
            key,P=canonicalize(S,M);poolkey=(q,len(M),key);stats['marked_inputs']+=1
            if poolkey not in pools:
                record=matrix_counts(S,list(P),validate_matrix=False)
                results={n:conditional_moments(record,n) for n in [7,8]}
                first={'graph':name,'marks':list(P),'record':record,'moments':results}
                pools[poolkey]=first;signatures[poolkey].add(rel_signature(record));stats['new_profiles']+=1;continue
            record=matrix_counts(S,list(P),validate_matrix=False);sig=rel_signature(record)
            if sig in signatures[poolkey]:stats['duplicate_canonical_relations']+=1;continue
            signatures[poolkey].add(sig);stats['new_relations_same_profile']+=1
            if args.cross_twins and pools[poolkey]['graph']==name:continue
            results={n:conditional_moments(record,n) for n in [7,8]};first=pools[poolkey]
            for n in [7,8]:
                if first['moments'][n]['second_moment']!=results[n]['second_moment']:
                    pair=[first,{'graph':name,'marks':list(P),'record':record,'moments':results}]
                    for row in pair:row['moments'][6]=conditional_moments(row['record'],6)
                    assert pair[0]['moments'][6]['second_moment']==pair[1]['moments'][6]['second_moment']
                    out={'status':'actual same-cell-profile different second moments found','n':n,'q':q,'marks_count':len(M),
                         'canonical_cell_profile':key,'pair':pair,'scanned_before_witness':census+[{'graph':name,**stats}],
                         'input_twin_sha256':sha(twin_path),'compiler_source_sha256':sha(ROUND/'marked_moments/conditional_moments.py'),
                         'source_sha256':sha(__file__),'elapsed_seconds':time.monotonic()-started}
                    (HERE/(prefix+'witness.json')).write_text(json.dumps(out,indent=2)+'\n')
                    # Export standalone actual matrices for a direct oracle.
                    matrices={g[0]:g[1] for g in graphs}
                    (HERE/(prefix+'witness_matrices.json')).write_text(json.dumps({'graphs':[{'name':r['graph'],'S':matrices[r['graph']]} for r in pair]},indent=2)+'\n')
                    print(json.dumps({'found':True,'q':q,'n':n,'left':[pair[0]['graph'],pair[0]['marks'],pair[0]['moments'][n]['second_moment']],
                                     'right':[pair[1]['graph'],pair[1]['marks'],pair[1]['moments'][n]['second_moment']], 'seconds':out['elapsed_seconds']}),flush=True)
                    return
        census.append({'graph':name,**stats});print(json.dumps(census[-1]),flush=True)
        (HERE/(prefix+'progress.json')).write_text(json.dumps(census,indent=2)+'\n')
    (HERE/(prefix+'no_triple_witness.json')).write_text(json.dumps({'status':'no triple ablation witness in specified scan','census':census,
        'elapsed_seconds':time.monotonic()-started,'source_sha256':sha(__file__)},indent=2)+'\n')


if __name__=='__main__':main()
