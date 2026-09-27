#!/usr/bin/env python3
"""Test supplied-space pruning against complete small and structural large oracles."""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import perf_counter
from pruning_sampler import prepare,run,trials_for

HERE=Path(__file__).resolve().parent


def main():
    smallpath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';small=json.loads(smallpath.read_text())
    smallcert=HERE/'hard_rank3.certificate.json';c=json.loads(smallcert.read_text());prepared=prepare(c,failure_bits=46);data=c['input'];p=data['p'];domain=data['domain'];basis=data['basis'];a=data['origin']
    # Independent full-space oracle from explicit polynomial coefficients.
    all_words=[]
    for parameters in product(range(p),repeat=3):
        coefficients=[((a[j] if j<len(a) else 0)+sum(parameters[i]*(basis[i][j] if j<len(basis[i]) else 0) for i in range(3)))%p for j in range(data['k'])]
        word=[sum(v*pow(x,j,p) for j,v in enumerate(coefficients))%p for x in domain]
        all_words.append((parameters,word))
    small_runs=[]
    for scalar in range(p):
        received=[(u+scalar*v)%p for u,v in zip(small['u0'],small['u1'])]
        truth=[list(parameters) for parameters,word in all_words if sum(a==b for a,b in zip(word,received))>=data['s']]
        start=perf_counter();result=run(prepared,received);seconds=perf_counter()-start
        assert [r['parameters'] for r in result['outputs']]==truth
        small_runs.append({'scalar':scalar,'run':result,'oracle_parameters':truth,'seconds':seconds})
    assert len(small_runs)*(1<<40)<=1<<46
    print(json.dumps({'phase':'small_complete','received_words':len(small_runs),'ambient_space_members':len(all_words),'trials_per_word':prepared['budget']['trials'],'outputs':[len(r['run']['outputs']) for r in small_runs]}),flush=True)
    largepath=HERE.parent/'round6/pencil_tracks/results.json';large=next(x['result'] for x in json.loads(largepath.read_text())['large'] if x['name']=='skew-two')
    largecert=HERE/'full_length.certificate.json';bigc=json.loads(largecert.read_text());bigprepared=prepare(bigc,failure_bits=42);bigdata=bigc['input'];p=bigdata['p'];k=bigdata['k'];domain=bigdata['domain'];large_runs=[]
    for scalar in (0,1):
        received=[(u+scalar*v)%p for u,v in zip(large['u0'],large['u1'])];pieces=[];supports=[]
        for track in large['verified_tracks']:
            a,b=track['intercept_coefficients'],track['slope_coefficients']
            poly=[((a[j] if j<len(a) else 0)+scalar*(b[j] if j<len(b) else 0))%p for j in range(k)]
            values=[sum(v*pow(x,j,p) for j,v in enumerate(poly))%p for x in domain]
            support=[i for i,(v,w) in enumerate(zip(values,received)) if v==w];assert len(support)>=bigdata['s']
            pieces.append(poly);supports.append(support)
        assert len(pieces)==2 and pieces[0]!=pieces[1] and set(supports[0])|set(supports[1])==set(range(bigdata['n']))
        assert 2*(k-1)<bigdata['s']
        # Any degree-<k polynomial different from both pieces has at most k-1
        # agreement roots on each piece, hence fewer than s agreements overall.
        start=perf_counter();result=run(bigprepared,received);seconds=perf_counter()-start
        assert sorted(r['coefficients'] for r in result['outputs'])==sorted(pieces)
        large_runs.append({'scalar':scalar,'run':result,'oracle_polynomials':pieces,'oracle_agreement_supports':supports,
                           'other_degree_bounded_polynomial_agreement_cap':2*(k-1),'seconds':seconds})
    assert len(large_runs)*(1<<40)<=1<<42
    # Compare the certified input-specific probability with the universal
    # common-root bound, keeping every other sampling assumption fixed.
    from math import comb
    from fractions import Fraction
    baseline=Fraction(comb(bigdata['s']-bigdata['k']+3,3),comb(bigdata['n'],3))
    baseline_budget=trials_for(p,[baseline.numerator,baseline.denominator],42)
    assert Fraction(len(small_runs),1<<46)+Fraction(len(large_runs),1<<42)<=Fraction(1,1<<40)
    out={'status':'passed','scope':'Actual supplied-space pruning. All17 small received scalars checked against exhaustive space membership; two large scalars checked against a complete two-polynomial root-bound oracle. This is not discovery of the needed clusters or a new general decoder.',
         'small_runs':small_runs,'large_runs':large_runs,'small_full_space_members_enumerated':len(all_words),
         'bundle_failure_probability_at_most_under_uniform_sampling':[1,1<<40],
         'large_certificate_trial_budget':bigprepared['budget'],'large_universal_baseline_trial_budget':baseline_budget,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_pruning.py','pruning_sampler.py','certificate_verifier.py')},
         'input_sha256':{'hard_rank3.certificate.json':sha256(smallcert.read_bytes()).hexdigest(),'full_length.certificate.json':sha256(largecert.read_bytes()).hexdigest(),
                         '../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(smallpath.read_bytes()).hexdigest(),
                         '../round6/pencil_tracks/results.json':sha256(largepath.read_bytes()).hexdigest()}}
    (HERE/'pruning_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({'status':'passed','small_received_words':len(small_runs),'large_received_words':len(large_runs),
                      'large_trials':bigprepared['budget']['trials'],'large_baseline_trials':baseline_budget['trials'],
                      'large_output_counts':[len(r['run']['outputs']) for r in large_runs]}),flush=True)


if __name__=='__main__':main()
