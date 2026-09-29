#!/usr/bin/env python3
"""Independent tiny-field oracles for bounded monic piece discovery."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from copy import deepcopy
from hashlib import sha256
from itertools import product
from pathlib import Path
import importlib.util
import json
import random
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'piece_discovery'/'piece_discovery.py'
FIELD_ORACLE=HERE.parents[1]/'round4'/'novelty'/'review_support_batches.py'


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def trim(a):
    while len(a)>1 and not a[-1]:a.pop()
    return a


def padd(F,a,b):
    out=[0]*max(len(a),len(b))
    for i in range(len(out)):out[i]=F.add(a[i] if i<len(a) else 0,b[i] if i<len(b) else 0)
    return trim(out)


def pmul(F,a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=F.add(out[i+j],F.mul(x,y))
    return trim(out)


def subst(F,relation,h):
    out=[0];power=[1]
    for a in relation:
        out=padd(F,out,pmul(F,a,power));power=pmul(F,power,h)
    return out


def factor_product(F,branches):
    out=[[1]]
    for h in branches:
        nxt=[[0] for _ in range(len(out)+1)]
        for j,a in enumerate(out):
            nxt[j]=padd(F,nxt[j],pmul(F,a,[F.neg(x) for x in h]))
            nxt[j+1]=padd(F,nxt[j+1],a)
        out=nxt
    return out


def rank(F,rows):
    """Independent full reduced elimination, rather than forward echelon."""
    a=[list(row) for row in rows];height=len(a);width=len(a[0]);r=0
    for c in range(width):
        pivot=next((i for i in range(r,height) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        inv=F.inv(a[r][c]);a[r]=[F.mul(v,inv) for v in a[r]]
        for i in range(height):
            if i!=r and a[i][c]:
                value=a[i][c];a[i]=[F.sub(x,F.mul(value,y)) for x,y in zip(a[i],a[r])]
        r+=1
        if r==height:break
    return r


def dot(F,a,b):
    out=0
    for x,y in zip(a,b):out=F.add(out,F.mul(x,y))
    return out


def matrix(F,dom,word,k,r):
    pairs=[(i,j) for j in range(r) for i in range((r-j)*(k-1)+1)]
    rows=[];rhs=[]
    for x,y in zip(dom,word):
        xp=[1];yp=[1]
        for _ in range(r*(k-1)):xp.append(F.mul(xp[-1],x))
        for _ in range(r):yp.append(F.mul(yp[-1],y))
        rows.append([F.mul(xp[i],yp[j]) for i,j in pairs]);rhs.append(F.neg(yp[r]))
    return rows,rhs


def audit_record(F,result,oracle,candidate,counts):
    dom,word,k,s=result['domain'],result['word'],result['k'],result['s']
    book=[];masks=set()
    for cs in product(range(F.q),repeat=k):
        cw=tuple(F.eval(cs,x) for x in dom)
        support=tuple(i for i,v in enumerate(cw) if v==word[i])
        book.append((cw,support));masks.add(sum(1<<i for i in support))
    distances=[len(dom)+1]*(1<<len(dom));distances[0]=0
    for covered in range(len(distances)):
        for mask in masks:
            total=covered|mask
            distances[total]=min(distances[total],distances[covered]+1)
    minimum=distances[-1]
    assert result['polynomial_cover_size_lower_bound']<=minimum
    if result['minimum_polynomial_cover_size_certified'] is not None:
        assert result['minimum_polynomial_cover_size_certified']==minimum
        counts['exact_minimum_cover_sizes_checked']+=1
    counts['minimum_cover_lower_bounds_checked']+=1
    for attempt in result['attempts']:
        if 'linear_certificate' not in attempt:continue
        r=attempt['piece_bound'];A,b=matrix(F,dom,word,k,r);cert=attempt['linear_certificate']
        mr=rank(F,A);ar=rank(F,[row+[v] for row,v in zip(A,b)])
        assert cert['rank']==mr and (cert['status']=='consistent')==(mr==ar)
        assert len(cert['kernel_basis'])==len(A[0])-mr
        for v in cert['kernel_basis']:assert all(dot(F,row,v)==0 for row in A)
        if cert['status']=='inconsistent':
            w=cert['inconsistency_witness'];i=w['row'];indices=w['basis_rows'];weights=w['coefficients']
            assert all(A[i][j]==dot(F,weights,[A[z][j] for z in indices]) for j in range(len(A[0])))
            residue=F.sub(b[i],dot(F,weights,[b[z] for z in indices]))
            assert residue!=0 and residue==w['rhs_residual'];counts['dual_inconsistency_checks']+=1
        else:
            assert all(dot(F,row,cert['particular'])==v for row,v in zip(A,b))
            relation=attempt['relation_y_coefficients']
            assert relation[-1]==[1] and len(relation)==r+1
            assert all(len(a)<=((r-j)*(k-1)+1) for j,a in enumerate(relation))
            assert all(F.eval([F.eval(a,x) for a in relation],y)==0 for x,y in zip(dom,word))
            factors=attempt['factor_recovery']
            for h in factors['branches']:
                assert len(h)==k and subst(F,relation,h)==[0]
                counts['full_branch_identities']+=1
            if factors.get('factor_product_verified'):
                assert factor_product(F,factors['branches'])==relation
                counts['complete_factor_products']+=1
        counts['independent_linear_rank_checks']+=1
    if result['complete_static_list_certified']:
        wanted={(cw,support) for cw,support in book if len(support)>=s}
        actual={(tuple(r['codeword']),tuple(r['agreement_support'])) for r in result['static_certificate']['complete_static_list']}
        assert actual==wanted
        counts['complete_static_list_oracles']+=1
    if result['complete_scalar_fibers_certified']:
        record=result['scalar_certificate'];nodes=set()
        for z in range(F.q):
            if z==record['anchor']['scalar']:
                a=record['anchor'];nodes.add((z,tuple(a['codeword']),tuple(a['agreement_support'])))
            else:
                for line in record['nonanchor_fibers']:
                    cw=tuple(F.add(a,F.mul(z,b)) for a,b in zip(line['intercept_codeword'],line['slope_codeword']))
                    support=tuple(i for i,v in enumerate(cw) if v==F.add(result['u0'][i],F.mul(z,word[i])))
                    assert list(support)==line['agreement_support_for_every_scalar_in_domain']
                    nodes.add((z,cw,support))
        assert nodes==oracle.complete_nodes(F,dom,k,s,result['u0'],word)
        assert len(nodes)==record['complete_symbolic_node_count']
        counts['complete_scalar_fiber_oracles']+=1;counts['complete_nodes_compared']+=len(nodes)
        counts['whole_field_cases']+=record['bad_scalar_count']==F.q
    counts['cases']+=1;counts['status:'+result['status']]+=1


def main():
    started=time.perf_counter();digest=sha256(SOURCE.read_bytes()).hexdigest()
    candidate=load('reviewed_piece_discovery',SOURCE)
    oracle=load('independent_field_oracle',FIELD_ORACLE)
    counts=Counter();records=[];rng=random.Random(493005);samples={}
    for q in (3,5,9):
        F=oracle.OracleField(q)
        CF=candidate.Field(3,[1,0,1]) if q==9 else candidate.PrimeField(q)
        if q==3:
            n=3;cases=[(list(w),k,s) for w in product(range(q),repeat=n) for k in (1,2) for s in range(k,n+1)]
        elif q==5:
            n=4;cases=[(list(w),2,3) for w in product(range(q),repeat=n)]
        else:
            n=7;cases=[]
            # Every assignment of coordinates to two fixed affine pieces.
            pieces=([1,1],[3,4])
            for labels in product((0,1),repeat=n):
                cases.append(([F.eval(pieces[j],i) for i,j in enumerate(labels)],2,4))
            cases += [([rng.randrange(q) for _ in range(n)],2,4) for _ in range(40)]
        for case_index,(word,k,s) in enumerate(cases):
            dom=list(range(n));u0=None
            if case_index%5==0:
                cs=[rng.randrange(q) for _ in range(k)];zstar=rng.randrange(q)
                u0=[F.sub(F.eval(cs,x),F.mul(zstar,y)) for x,y in zip(dom,word)]
            result=candidate.discover(CF,dom,k,s,word,u0,max_pieces=3,center_budget=24,split_budget=64)
            audit_record(F,result,oracle,candidate,counts)
            samples.setdefault((q,result['status']),result)
            counts[f'field_{q}_cases']+=1
        records.append({'q':q,'n':n,'cases':len(cases)})
    # Distinct polynomial branches can still have no simple split center.
    F=oracle.OracleField(3);CF=candidate.PrimeField(3)
    collision=factor_product(F,[[0,1],[0,2],[1]])
    got=candidate.recover_branches(CF,collision,2)
    assert got['status']=='no_simple_split_center_within_budget' and got['branches']==[]
    # A valid local jet need not be a polynomial root: Y²=1+X.
    CF=candidate.PrimeField(5)
    nonpoly=candidate.recover_branches(CF,[[4,4],[0],[1]],2)
    assert nonpoly['status']=='some_simple_branches_nonpolynomial' and nonpoly['branches']==[]
    assert not nonpoly['factor_product_verified']
    for (q,status),record in samples.items():
        CF=candidate.Field(3,[1,0,1]) if q==9 else candidate.PrimeField(q)
        assert candidate.verify_export(CF,record)
        counts['fresh_export_readbacks']+=1
    scalar=next(v for v in samples.values() if v['complete_scalar_fibers_certified'])
    inconsistent=next(v for v in samples.values() if any(a.get('linear_certificate',{}).get('status')=='inconsistent' for a in v['attempts']))
    mutations=[]
    bad=deepcopy(scalar);bad['minimum_polynomial_cover_size_certified']=99;mutations.append(('false_minimum_cover',bad))
    bad=deepcopy(scalar);bad['scalar_certificate']['complete_symbolic_node_count']+=1;mutations.append(('false_node_count',bad))
    bad=deepcopy(scalar)
    attempt=next(a for a in bad['attempts'] if 'relation_y_coefficients' in a)
    attempt['relation_y_coefficients'][-1]=[0];mutations.append(('nonmonic_relation',bad))
    bad=deepcopy(scalar)
    attempt=next(a for a in bad['attempts'] if a.get('recovered_maximal_supports'))
    attempt['recovered_maximal_supports'][0].clear();mutations.append(('false_support',bad))
    bad=deepcopy(inconsistent)
    attempt=next(a for a in bad['attempts'] if a.get('linear_certificate',{}).get('status')=='inconsistent')
    attempt['linear_certificate']['inconsistency_witness']['rhs_residual']=0;mutations.append(('false_dual_residual',bad))
    rejected=[]
    for name,bad in mutations:
        q=bad['field']['q'];CF=candidate.Field(3,[1,0,1]) if q==9 else candidate.PrimeField(q)
        try:candidate.verify_export(CF,bad)
        except AssertionError:rejected.append(name)
        else:raise AssertionError(('corrupt export accepted',name))
    assert sha256(SOURCE.read_bytes()).hexdigest()==digest,'source changed during review'
    out={'status':'all successful bounded discoveries match independent complete lists',
         'counts':dict(counts),'field_cases':records,
         'mutated_exports_rejected':rejected,
         'intentional_incompleteness_controls':{'three_distinct_polynomial_factors_without_simple_center':got,
                                                'nonpolynomial_Hensel_jets_rejected':nonpoly},
         'source_sha256':{'../piece_discovery/piece_discovery.py':digest,
                          Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest(),
                          '../../round4/novelty/review_support_batches.py':sha256(FIELD_ORACLE.read_bytes()).hexdigest()},
         'elapsed_seconds':time.perf_counter()-started}
    (HERE/'review_piece_discovery.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
