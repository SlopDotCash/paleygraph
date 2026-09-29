#!/usr/bin/env python3
"""Independent certificate/ledger audit, without production algebra imports.

Reads all saved masks, but does not independently recompute every large mask.
The saved complete barycentric/Newton replay supplies that layer explicitly.
"""
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb,isqrt
from pathlib import Path
import random
import struct

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
PRODUCER=LAB/'round6/cover_barriers'


def digest(path):return sha256(path.read_bytes()).hexdigest()


def packed(n,k,M):
    return n//M*min(M,k-1)+min(n%M,k-1)


def partitions(n,largest):
    if n==0:
        yield ()
    else:
        for first in range(min(n,largest),0,-1):
            for rest in partitions(n-first,first):yield (first,)+rest


def evaluate(cs,x,p):return sum(c*pow(x,j,p) for j,c in enumerate(cs))%p


def lagrange(base,word,dom,p,x):
    result=0
    for i in base:
        numerator=word[i];denominator=1
        for j in base:
            if i!=j:
                numerator=numerator*(x-dom[j])%p
                denominator=denominator*(dom[i]-dom[j])%p
        result=(result+numerator*pow(denominator,-1,p))%p
    return result


def main():
    result_path=PRODUCER/'results.json';d=json.loads(result_path.read_text())
    check_path=PRODUCER/'verification.json';saved=json.loads(check_path.read_text())
    assert saved['status']=='passed' and saved['results_sha256']==digest(result_path)
    bindings={}
    for record in (d,saved):
        for key,folder in (('source_sha256',PRODUCER),('artifact_sha256',PRODUCER),('dependency_sha256',LAB)):
            for name,h in record.get(key,{}).items():
                path=folder/name;assert digest(path)==h,path
                bindings[str(path.relative_to(LAB))]=h
    formula_checks=0
    for n in range(1,21):
        all_parts=list(partitions(n,n))
        for k in range(1,n+1):
            for M in range(1,n+1):
                expected=min(sum(min(x,k-1) for x in parts) for parts in all_parts if max(parts)<=M)
                assert packed(n,k,M)==expected
                formula_checks+=1

    # Direct full-codebook maxima versus independent base interpolation,
    # including every prime-field length3 pencil at k1.
    tiny=0
    for p,n,k,count in ((3,3,1,729),(5,5,2,32)):
        rng=random.Random(924+p)
        book=[tuple(evaluate(cs,x,p) for x in range(n)) for cs in product(range(p),repeat=k)]
        pencils=product(range(p),repeat=2*n) if p==3 else ([rng.randrange(p) for _ in range(2*n)] for _ in range(count))
        for pencil in pencils:
            words=[pencil[:n],pencil[n:]]
            masks=[[sum(1<<i for i,(a,b) in enumerate(zip(cw,w)) if a==b) for cw in book] for w in words]
            direct=[max(mask.bit_count() for mask in row) for row in masks]
            direct.append(max((a&b).bit_count() for a in masks[0] for b in masks[1]))
            found=[0,0,0]
            for base in combinations(range(n),k):
                pair=[sum(1<<i for i,x in enumerate(range(n)) if lagrange(base,w,list(range(n)),p,x)==w[i]) for w in words]
                for j,mask in enumerate(pair+[pair[0]&pair[1]]):found[j]=max(found[j],mask.bit_count())
            assert found==direct
            tiny+=1

    ledger_rows=0;sampled_scalar_masks=0;maximum_witnesses=0;attainments=0
    groups={}
    for row in d['hard_stacks']:groups.setdefault(row['configuration'],[]).append(row)
    for name,records in groups.items():
        records.sort(key=lambda r:r['sample'])
        n,k,p=records[0]['n'],records[0]['k'],records[0]['field']['p']
        assert all(p%j for j in range(2,isqrt(p)+1))
        count=comb(n,k)
        generation=json.loads((PRODUCER/f'{name}.enumerated.json').read_text())
        replay=json.loads((PRODUCER/f'{name}.verified.json').read_text())
        assert generation['basis_count']==replay['basis_count']==count
        assert generation['samples']==replay['samples']
        assert generation['method']=='barycentric enumeration' and replay['method']=='Newton replay'
        data=(PRODUCER/records[0]['ledger']['path']).read_bytes()
        assert len(data)==count*len(records)*12
        iterator=struct.iter_unpack('<III',data)
        maxima=[[0]*3 for _ in records]
        for index,base in enumerate(combinations(range(n),k)):
            base_mask=sum(1<<i for i in base)
            for j,row in enumerate(records):
                masks=next(iterator)
                assert masks[2]==masks[0]&masks[1]
                assert all(mask<(1<<n) and mask&base_mask==base_mask for mask in masks)
                for mode,mask in enumerate(masks):maxima[j][mode]=max(maxima[j][mode],mask.bit_count())
                if index in (0,1,count//3,count//2,count-2,count-1):
                    for mode in (0,1):
                        word=row['u0' if mode==0 else 'u1'];dom=row['domain']
                        literal=sum(1<<i for i,x in enumerate(dom) if lagrange(base,word,dom,p,x)==word[i])
                        assert literal==masks[mode]
                        sampled_scalar_masks+=1
                ledger_rows+=1
        for row,got in zip(records,maxima):
            assert got==row['exact_maxima_u0_u1_joint']
            for mode,witness in enumerate(row['maximum_witnesses']):
                words=[row['u0'],row['u1']] if mode==2 else [row['u0' if mode==0 else 'u1']]
                coefficients=witness['coefficients']
                assert all(len(cs)<=k and all(0<=c<p for c in cs) for cs in coefficients)
                support=[i for i,x in enumerate(row['domain']) if all(evaluate(cs,x,p)==w[i] for cs,w in zip(coefficients,words))]
                assert support==witness['support'] and len(support)==got[mode]
                maximum_witnesses+=1
                partition=row['actual_cap_minimum_attainment_u0_u1_joint'][mode]
                used=[]
                for piece in partition['pieces']:
                    assert len(piece['coefficients'])==len(words)
                    assert all(len(cs)<=k for cs in piece['coefficients'])
                    for i in piece['coordinates']:
                        assert all(evaluate(cs,row['domain'][i],p)==w[i] for cs,w in zip(piece['coefficients'],words))
                    used+=piece['coordinates']
                assert sorted(used)==list(range(n))
                cap=sum(min(len(piece['coordinates']),k-1) for piece in partition['pieces'])
                assert cap==packed(n,k,got[mode])==partition['cap']
                assert row['barriers_u0_u1_joint'][mode]['minimum']==cap
                attainments+=1
            assert packed(n,k,got[2])>=row['s']
    assert ledger_rows==d['hard_stack_bases_exhausted']==893464
    rate_checks=0
    for k in range(3,70):
        for delta in range(k-2):
            threshold=2*k-(delta+1)//2
            assert packed(4*k,k,threshold)<2*k+delta
            assert packed(4*k,k,threshold-1)>=2*k+delta
            rate_checks+=1
    for row in d['rate_quarter']:
        n,k,s,M=[row[key] for key in ('n','k','s','minimum_M_not_blocked_by_abstract_cap')]
        assert packed(n,k,M)<s
        if M>1:assert packed(n,k,M-1)>=s
    report={'status':'passed','integer_partition_checks':formula_checks,'complete_tiny_prime_pencil_oracles':tiny,
            'saved_ledger_sample_base_rows_read':ledger_rows,'independent_scalar_mask_samples':sampled_scalar_masks,
            'exact_maximum_witnesses':maximum_witnesses,'actual_optimal_cap_attainments':attainments,
            'joint_family_obstructions':len(d['hard_stacks']),'rate_quarter_boundary_checks':rate_checks,
            'scope':'All saved masks checked for structure/maxima and source-bound to completed barycentric/Newton replay. Six bases per configuration independently re-evaluated; no third complete large interpolation census claimed.',
            'source_sha256':{**bindings,str(result_path.relative_to(LAB)):digest(result_path),
                            str(check_path.relative_to(LAB)):digest(check_path),
                            str(Path(__file__).relative_to(LAB)):digest(Path(__file__))}}
    (HERE/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
