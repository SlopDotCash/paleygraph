#!/usr/bin/env python3
"""Reproduce round1 equality, extension controls, failure/break-even and larger exact cases."""
from itertools import combinations,product
from math import comb
from pathlib import Path
import json
import random
import time

from compressed_support import *
from extension_field_probe import subgroup
from stack_provenance import parity_catalog,census

HERE=Path(__file__).parent


def planted(F,dom,k,s,seed):
    rng=random.Random(seed)
    cs0=[rng.randrange(F.q) for _ in range(k)]
    cs1=[rng.randrange(F.q) for _ in range(k)]
    f0=[F.sum(F.mul(c,F.power(x,j)) for j,c in enumerate(cs0)) for x in dom]
    f1=[F.sum(F.mul(c,F.power(x,j)) for j,c in enumerate(cs1)) for x in dom]
    u0=[rng.randrange(F.q) for _ in dom]
    u1=[rng.randrange(F.q) for _ in dom]
    for i in range(s):
        u0[i]=f0[i]
    for i in range(len(dom)-s,len(dom)):
        u1[i]=F.sub(f1[i],u0[i])
    return u0,u1


def tiny_cover_audit():
    count=0
    for n in range(2,11):
        for k in range(1,n+1):
            for s in range(k,n+1):
                plan=cover_plan(n,k,s)
                verify_cover(plan)
                bases=[set(B) for block in plan['blocks'] for B in combinations(block,k)]
                for S in combinations(range(n),s):
                    assert any(B<=set(S) for B in bases)
                    count+=1
    return count


def main():
    start=time.perf_counter()
    out={'schema':'compressed-agreement-cover/v1','scope':'complete per-stack finite-affine census, no worst-case bound',
         'tiny_cover_supports_exhaustively_checked':tiny_cover_audit(),'round1_comparisons':[],
         'extension_comparisons':[],'scale_experiments':[],'costs':[]}
    matrix_checks=0
    for F,dom in [(PrimeField(17),domain(8,17)),(Field(3,[1,0,1]),[0,1,2,3,4,5,6,7,8])]:
        for k in [1,2,3]:
            for base in combinations(dom,k):
                assert interpolation_matrix(F,list(base),dom)==interpolation_matrix_reference(F,list(base),dom)
                matrix_checks+=1
    out['optimized_reference_lagrange_matrix_equalities']=matrix_checks
    old=json.loads((ROUND1/'stack_results.json').read_text())
    for r in old['records']:
        p=r['p']; dom=domain(r['n'],p)
        result=run_census(PrimeField(p),dom,r['k'],r['s'],r['u0'],r['u1'])
        assert node_set(result)==node_set(r)
        assert result['finite_bad_scalar_count']==r['finite_bad_scalar_count']
        assert result['whole_field_correlated']==r['whole_field_correlated']
        out['round1_comparisons'].append({'configuration':r['configuration'],'sample':r['sample'],
            'full_subset_oracle_count':comb(r['n'],r['s']),'full_node_equality':True,'result':result})
        print('round1',r['configuration'],r['sample'],'bases',result['stats']['bases_processed'],
              'bad',result['finite_bad_scalar_count'],'seconds',result['elapsed_seconds'],flush=True)
    ext=json.loads((ROUND1/'extension_results.json').read_text())
    for r in ext['experiments']:
        F=Field(r['p'],r['modulus'])
        result=run_census(F,r['domain'],r['k'],r['s'],r['u0'],r['u1'])
        assert node_set(result)==node_set(r)
        assert result['bad_scalars']==r['bad_scalars']
        assert result['whole_field_correlated']==bool(r['common_subsets'])
        out['extension_comparisons'].append({'name':r['name'],'full_node_equality':True,'result':result})
        print('extension',r['name'],'bad',result['finite_bad_scalar_count'],flush=True)
    # Measured failure: enumerate ALL k-point bases. It is complete but needlessly costly.
    r=old['records'][4]; F=PrimeField(r['p']); dom=domain(r['n'],r['p'])
    naive=run_census(F,dom,r['k'],r['s'],r['u0'],r['u1'],mode='full_k')
    t=time.perf_counter(); cat=parity_catalog(dom,r['k'],r['s'],r['p'])
    full=census(r['u0'],r['u1'],dom,r['k'],r['s'],r['p'],cat)
    full_seconds=time.perf_counter()-t
    assert node_set(naive)==node_set(r)
    out['measured_naive_failure']={'n':r['n'],'k':r['k'],'s':r['s'],
        'full_k_basis_count':naive['stats']['bases_processed'],'full_k_seconds':naive['elapsed_seconds'],
        'full_s_subset_count':len(cat),'full_s_seconds':round(full_seconds,6),
        'all_modes_same_nodes':True}
    configs=[(PrimeField(65537),32,8,20),(PrimeField(65537),64,4,40),
             (PrimeField(65537),64,4,12),
             (PrimeField(65537),128,8,80),(Field(3,[1,0,0,0,1,1,1]),56,4,35)]
    for F,n,k,s in configs:
        dom=domain(n,F.q) if F.e==1 else subgroup(F,n)
        u0,u1=planted(F,dom,k,s,12000+n)
        result=run_census(F,dom,k,s,u0,u1)
        assert 0 in result['bad_scalars'] and 1 in result['bad_scalars']
        verify_cover(result['cover_certificate'])
        # A coordinate rotation chooses a different complete cover. Equality is
        # a useful independent covering choice, not an independent decoding implementation.
        rot=7%n; perm=list(range(rot,n))+list(range(rot))
        second=run_census(F,[dom[i] for i in perm],k,s,[u0[i] for i in perm],[u1[i] for i in perm])
        restored=[]
        for node in second['nodes']:
            cw=[0]*n
            for i,original in enumerate(perm): cw[original]=node['codeword'][i]
            restored.append({'scalar':node['scalar'],'codeword':cw,
                             'agreement_support':sorted(perm[i] for i in node['agreement_support'])})
        assert node_set(result)==node_set({'nodes':restored})
        out['scale_experiments'].append({'seed':12000+n,'certificate_complete':True,
            'full_subset_count_not_enumerated':comb(n,s),'rotated_cover_same_all_nodes':True,
            'result':result,'rotated_cover_seconds':second['elapsed_seconds']})
        print('scale',F.q,n,k,s,'bases',result['stats']['bases_processed'],'bad',result['finite_bad_scalar_count'],
              'nodes',len(result['nodes']),'seconds',result['elapsed_seconds'],flush=True)
    for n,k,s in [(16,8,11),(20,10,14),(32,8,20),(32,8,16),(32,8,12),
                  (64,4,40),(128,8,80),(1024,64,410)]:
        plan=cover_plan(n,k,s)
        out['costs'].append({'n':n,'k':k,'s':s,'full_s_subsets':comb(n,s),
            'naive_full_k_bases':comb(n,k),'naive_full_k_plus_1_bases':comb(n,k+1),
            'anchor_bases':cover_plan(n,k,s,'anchor')['basis_count'],
            'partition_bases':plan['basis_count'],'partition_sizes':[len(b) for b in plan['blocks']]})
    out['elapsed_seconds']=round(time.perf_counter()-start,3)
    dest=HERE/'results.json'; dest.write_text(json.dumps(out,indent=2)+'\n')
    print('Saved',dest,'in',out['elapsed_seconds'],'seconds',flush=True)


if __name__=='__main__':main()
