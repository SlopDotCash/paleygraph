#!/usr/bin/env python3
"""Source-bound rational readback; shares model helpers, never invokes HiGHS."""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from trace_envelopes import problem, verify

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]


def digest(path):return sha256(path.read_bytes()).hexdigest()


def main():
    result_path=HERE/'results.json'
    result=json.loads(result_path.read_text())
    assert digest(HERE/'trace_envelopes.py')==result['source_sha256']
    for name,value in result['input_sha256'].items():
        assert digest(LAB/name)==value
    moments=json.loads((LAB/'round5/trace_invariants/results.json').read_text())['comparisons']
    inequalities=0;equalities=0;mutations=0
    for row in result['cases']:
        path=LAB/row['inventory'];assert digest(path)==row['inventory_sha256']
        moment=next(x for x in moments if (x['q'],x['n'])==(row['q'],row['n']))
        data=problem(moment,json.loads(path.read_text()),row['model'])
        assert row['support_points']==[list(point) for point in data['points']]
        for name in ('lower_certificate','upper_certificate'):
            cert=row[name]
            value=verify(data,cert)
            inequalities+=len(data['points'])
            equalities+=len(data['rows']) if cert['primal'] is not None else 0
            bound=data['constant']+cert['sign']*value
            assert bound==Fraction(*row['lower' if cert['sign']==-1 else 'upper'])
        assert Fraction(*row['lower'])<=Fraction(*row['actual'])<=Fraction(*row['upper'])
        width=Fraction(*row['upper'])-Fraction(*row['lower'])
        assert width==Fraction(*row['width'])
        assert width**2/Fraction(*row['variance'])**3==Fraction(*row['standardized_width_squared'])
        # Deliberate invalid proposals must not become certified results.
        for which in ('dual','primal'):
            damaged=deepcopy(row['upper_certificate'])
            if which=='dual':damaged['dual_multipliers'][0]=[-10**100,1]
            else:damaged['primal'][0]['weight']=[-1,1]
            try:verify(data,damaged)
            except AssertionError:mutations+=1
            else:raise AssertionError(f'accepted corrupted {which} certificate')
    output={'status':'passed','scope':'shared-model exact rational readback; no optimizer and no independent-model claim',
            'case_count':len(result['cases']),'pointwise_dual_inequalities':inequalities,
            'primal_equalities':equalities,'corrupted_certificates_rejected':mutations,
            'results_sha256':digest(result_path),'source_sha256':{name:digest(HERE/name)
                for name in ('trace_envelopes.py','verify_results.py')}}
    (HERE/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
