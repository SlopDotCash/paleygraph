#!/usr/bin/env python3
"""Bounded independent replay of discovery parameterization and cap failures."""
import sys
sys.dont_write_bytecode=True
from hashlib import sha256
from pathlib import Path
import json
import time
from review_tracks import load,SOURCE,ORACLE,audit

HERE=Path(__file__).resolve().parent


def check_partition(F,result):
    dom=result['domain'];u0,u1=result['u0'],result['u1'];tracks=result['verified_tracks'];k=result['k']
    assert sorted(i for t in tracks for i in t['coordinates'])==list(range(len(dom)))
    for t in tracks:
        for i in t['coordinates']:
            assert F.eval(t['intercept_coefficients'],dom[i])==u0[i]
            assert F.eval(t['slope_coefficients'],dom[i])==u1[i]
    cap=sum(min(len(t['coordinates']),k-1) for t in tracks)
    assert cap==result['joint_root_count_cap']
    return cap


def main():
    started=time.perf_counter();candidate=load('discovery_review_candidate',SOURCE)
    oracle=load('discovery_review_field',ORACLE)
    helper_path=HERE.parents[1]/'round5'/'review'/'review_piece_discovery.py'
    helper=load('independent_monic_matrix',helper_path)
    result_path=SOURCE.parent/'results.json';saved=json.loads(result_path.read_text());reports=[]
    first=saved['iterations'][0];base=first['initial'];F=oracle.OracleField(17);CF=candidate.PrimeField(17)
    for reassign,want in ((False,first['initial']),(True,first['iteration'])):
        actual=candidate.discover_pencil(CF,base['domain'],base['k'],base['s'],base['u0'],base['u1'],reassign=reassign)
        cap=check_partition(F,actual)
        assert actual['status']==want['status'] and actual['verified_tracks']==want['verified_tracks']
        if reassign:
            checked=audit(F,actual['certificate'],True)
            assert checked['nodes']==34 and cap==2
        else:assert cap==actual['s']==3 and not actual['complete']
        reports.append({'case':'ambiguity_fragment','reassignment':reassign,'cap':cap,'status':actual['status'],'fresh_discovery_replayed':True})
    second=saved['iterations'][1];initial,final=second['initial'],second['iteration']
    F=oracle.OracleField(257);F.inv=lambda a:pow(a,-1,257)
    assert check_partition(F,final)==28
    tracks=final['verified_tracks'];k=final['k'];dom=final['domain'];z=1
    polys=[tuple(F.add(a,F.mul(z,b)) for a,b in zip(t['intercept_coefficients'],t['slope_coefficients'])) for t in tracks]
    assert len(set(polys))==len(polys)==4
    word=[F.add(a,F.mul(z,b)) for a,b in zip(final['u0'],final['u1'])]
    supports=[[i for i,x in enumerate(dom) if F.eval(h,x)==word[i]] for h in polys]
    assert all(len(support)>3*(k-1) for support in supports)
    # Every putative three-polynomial cover must contain each of these four
    # polynomials: otherwise its three differences cover at most3(k-1) of
    # the corresponding matching coordinates. Thus minimum slice cover is4.
    bad_marginal=initial['marginal_discoveries'][1]
    assert bad_marginal['word']==word and bad_marginal['polynomial_cover_size_lower_bound']==4
    duals=0
    for attempt in bad_marginal['attempts']:
        r=attempt['piece_bound'];assert r in (1,2,3)
        A,b=helper.matrix(F,dom,word,k,r);linear=attempt['linear_certificate']
        assert linear['status']=='inconsistent';w=linear['inconsistency_witness']
        index=w['row'];basis=w['basis_rows'];weights=w['coefficients']
        assert all(A[index][j]==helper.dot(F,weights,[A[i][j] for i in basis]) for j in range(len(A[0])))
        residue=F.sub(b[index],helper.dot(F,weights,[b[i] for i in basis]))
        assert residue!=0 and residue==w['rhs_residual'];duals+=1
    assert duals==3 and not initial['complete']
    reports.append({'case':'choice_of_marginals','default_joint_cap':28,'default_symbolic_nodes':final['certificate']['complete_symbolic_node_count'],
                    'slice_scalar':1,'independent_exact_minimum_slice_cover':4,'slice_support_sizes':[len(x) for x in supports],
                    'independent_dual_inconsistency_witnesses':duals,'discovery_budget':3,
                    'fresh_large_discovery_replayed':False})
    negative=next(x['result'] for x in saved['negative_controls'] if x['name']=='nine-joint-regions-cap-failure')
    F=oracle.OracleField(257);F.inv=lambda a:pow(a,-1,257)
    assert check_partition(F,negative)==63 and negative['s']==50 and not negative['complete']
    try:candidate.compile_certificate(candidate.PrimeField(257),negative['domain'],negative['k'],negative['s'],negative['u0'],negative['u1'],negative['verified_tracks'])
    except candidate.UnavailableCertificate:pass
    else:raise AssertionError('non-strict cap accepted')
    reports.append({'case':'nine_joint_regions','verified_track_count':len(negative['verified_tracks']),'cap':63,'threshold':50,'complete_list_claim_rejected':True})
    out={'status':'discovery repair and failure claims independently verified','records':reports,
         'source_sha256':{'pencil_tracks.py':sha256(SOURCE.read_bytes()).hexdigest(),Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest(),
                          '../../round5/review/review_piece_discovery.py':sha256(helper_path.read_bytes()).hexdigest()},
         'results_sha256':sha256(result_path.read_bytes()).hexdigest(),'elapsed_seconds':time.perf_counter()-started}
    (HERE/'review_discovery.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
