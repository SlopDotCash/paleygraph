#!/usr/bin/env python3
"""Independent algebraic replay of saved n1024 discovered-cover certificates.

Uses plain field arithmetic and root-count uniqueness, not the discovery
elimination, factor recovery, or frozen list verifier.
"""
import sys
sys.dont_write_bytecode=True
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import time

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def main():
    started=time.perf_counter()
    helper=load('independent_piece_review',HERE/'review_piece_discovery.py')
    oracle=load('independent_field_oracle',helper.FIELD_ORACLE)
    saved=HERE.parent/'piece_discovery'/'results.json'
    data=json.loads(saved.read_text());reports=[]
    for entry in data['large']:
        out=entry['result'];F=oracle.OracleField(out['field']['q'])
        assert out['field']['e']==1
        dom,word,k,s=out['domain'],out['word'],out['k'],out['s'];n=len(dom)
        assert n==out['n'] and len(set(dom))==n
        assert all(0<=x<F.q for x in dom+word)
        pieces=out['discovered_pieces'];r=len(pieces)
        polys=[p['coefficients'] for p in pieces]
        assert len({tuple(h) for h in polys})==r
        assert all(len(h)==k and all(0<=x<F.q for x in h) for h in polys)
        supports=[[i for i,x in enumerate(dom) if F.eval(h,x)==word[i]] for h in polys]
        assert sorted(i for p in pieces for i in p['coordinates'])==list(range(n))
        assert all(set(p['coordinates'])<=set(supp) for p,supp in zip(pieces,supports))
        assert all(len(supp)>r*(k-1) for supp in supports)
        product_relation=helper.factor_product(F,polys)
        # Any competing monic relation of Y-degree j<=r and weighted degree
        # j(k-1) vanishes identically on each of these r polynomials: it has
        # more zeros than its X-degree. Distinct linear Y factors then force
        # j>=r, and monicity forces this product at j=r. The homogeneous
        # ansatz has Y-degree<j, giving zero nullity for all tested j<=r.
        checked_equations=0
        for attempt in out['attempts']:
            j=attempt['piece_bound'];cert=attempt['linear_certificate']
            assert j<=r
            expected_unknowns=(k-1)*j*(j+1)//2+j
            assert cert['rank']==cert['unknowns']==expected_unknowns
            assert cert['kernel_basis']==[] and cert['free_columns']==[]
            if j<r:
                assert cert['status']=='inconsistent'
            else:
                assert cert['status']=='consistent'
                assert attempt['relation_y_coefficients']==product_relation
                position=0
                for ydegree in range(j):
                    a=product_relation[ydegree]
                    for xdegree in range((j-ydegree)*(k-1)+1):
                        assert cert['particular'][position]==(a[xdegree] if xdegree<len(a) else 0)
                        position+=1
                assert position==expected_unknowns
                for x,y in zip(dom,word):
                    assert F.eval([F.eval(a,x) for a in product_relation],y)==0
                    checked_equations+=1
        assert out['minimum_polynomial_cover_size_certified']==r
        assert out['polynomial_cover_size_lower_bound']==r
        static=out['static_certificate']
        cap=sum(min(len(p['coordinates']),k-1) for p in pieces)
        assert static['outside_piece_agreement_cap']==cap and s>cap
        qualifying={(tuple(h),tuple(supp)) for h,supp in zip(polys,supports) if len(supp)>=s}
        actual={(tuple(p['coefficients']),tuple(p['agreement_support'])) for p in static['complete_static_list']}
        assert actual==qualifying
        for candidate in static['piece_candidates']:
            h=candidate['coefficients'];cw=[F.eval(h,x) for x in dom]
            supp=[i for i,v in enumerate(cw) if v==word[i]]
            assert candidate['codeword']==cw and candidate['agreement_support']==supp
            assert candidate['qualifies']==(len(supp)>=s)
        pencil=out['scalar_certificate'];anchor=pencil['anchor'];zstar=anchor['scalar']
        fcs=anchor['coefficients'];f=[F.eval(fcs,x) for x in dom]
        assert anchor['codeword']==f and anchor['agreement_support']==list(range(n))
        assert all(a==F.add(u,F.mul(zstar,w)) for a,u,w in zip(f,out['u0'],word))
        assert {(tuple(v['direction_coefficients']),tuple(v['agreement_support_for_every_scalar_in_domain'])) for v in pencil['nonanchor_fibers']}==qualifying
        for line in pencil['nonanchor_fibers']:
            h=line['direction_coefficients'];cw=[F.eval(h,x) for x in dom]
            assert line['scalar_domain']=={'kind':'all_except','excluded':[zstar]}
            assert line['direction_codeword']==line['slope_codeword']==cw
            assert line['slope_coefficients']==h
            assert line['intercept_coefficients']==[F.sub(a,F.mul(zstar,b)) for a,b in zip(fcs,h)]
            assert line['intercept_codeword']==[F.sub(a,F.mul(zstar,b)) for a,b in zip(f,cw)]
        count=len(qualifying)
        assert pencil['list_size_at_anchor']==1
        assert pencil['list_size_at_every_nonanchor_scalar']==count
        assert pencil['complete_symbolic_node_count']==1+(F.q-1)*count
        assert pencil['bad_scalar_count']==(F.q if count else 1)
        reports.append({'name':entry['name'],'q':F.q,'n':n,'k':k,'s':s,
                        'independently_certified_minimum_piece_count':r,
                        'maximal_support_sizes':[len(x) for x in supports],
                        'weighted_degree':r*(k-1),'monic_coefficient_equations_checked':checked_equations,
                        'rank_certification':'root-count factor divisibility gives unique monic relation and zero homogeneous kernel; no discovery solver reuse',
                        'complete_static_list_size':count,'symbolic_node_count':1+(F.q-1)*count})
    result={'status':'both large discovered covers and symbolic lists independently verified',
            'records':reports,'source_sha256':{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest(),
                'review_piece_discovery.py':sha256((HERE/'review_piece_discovery.py').read_bytes()).hexdigest(),
                '../../round4/novelty/review_support_batches.py':sha256(helper.FIELD_ORACLE.read_bytes()).hexdigest()},
            'results_sha256':sha256(saved.read_bytes()).hexdigest(),'elapsed_seconds':time.perf_counter()-started}
    (HERE/'review_piece_scale.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
