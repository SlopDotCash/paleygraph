#!/usr/bin/env python3
"""Independent acceptance using direct small subsets and the old twin census."""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    source=HERE/'third_moment_preflight.json'
    records=json.loads(source.read_text())
    old=HERE.parents[1]/'round3/pair_type_twins/results.json'
    twins=json.loads(old.read_text())['census']
    out=[]
    for row in records['cases']:
        q=row['q']
        if q in (13,17):
            chi=[0]+[1 if pow(x,(q-1)//2,q)==1 else -1 for x in range(1,q)]
            total=0
            for C in combinations(range(q),6):
                value=sum(prod(chi[(y-x)%q] for x in C) for y in range(q))
                total+=value**3
            expectation=Fraction(total,comb(q,6))
            method='fresh complete six-subset oracle using literal integer row products'
        else:
            assert row['graph'] in ('Paley49','Peisert49')
            histogram=twins['paley' if row['graph']=='Paley49' else 'peisert']
            total=sum(int(value)**3*count for value,count in histogram.items())
            expectation=Fraction(total,twins['six_sets'])
            method='independently verified full round3 histogram, all13983816 six-subsets'
        assert expectation==Fraction(*row['third_moment'])
        assert total==row['sum_T6_cubed']
        assert row['six_column_subsets_enumerated']==0
        out.append({'graph':row['graph'],'third_moment':[expectation.numerator,expectation.denominator],
                    'method':method})
    result={'status':'passed','cases':out,
            'preflight_sha256':sha256(source.read_bytes()).hexdigest(),
            'prior_twin_census_sha256':sha256(old.read_bytes()).hexdigest(),
            'review_source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'note':'The proposed method uses no six-set enumeration; this independent small acceptance checker deliberately does.'}
    (HERE/'root_acceptance_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Third-moment preflight acceptance passed on all four cases.')


if __name__=='__main__':main()
