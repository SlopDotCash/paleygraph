#!/usr/bin/env python3
"""Seven measured statistics for the degree-six global third moment.

The formal logarithm proves universal theta^6 and zero theta^5 terms.
The remaining two quartics require ten powers; conference identities leave
six measured low powers plus the total sixth power. These are sufficient
statistics from a trusted inventory, not a graph-realizability certificate.
"""
from fractions import Fraction
from hashlib import sha256
from math import factorial
from pathlib import Path
import time
from folded_moment import compiler,fold_inventory,unfold_inventory,polynomial_from_values,encode

STAT_NAMES=('mono_theta1','mono_theta2','mono_theta3','mixed_theta3','mono_theta4','mixed_theta4','total_theta6')


def universal_sixth_union():
    coefficients=[Fraction(1)]
    for _ in range(6):
        new=[Fraction() for _ in range(len(coefficients)+3)]
        for j,c in enumerate(coefficients):
            for k,v in ((1,1),(2,-3),(3,2)):new[j+k]+=v*c
        coefficients=new
    return [c/factorial(6) for c in coefficients]


def seven_statistics(record):
    folded=record if 'family_power_sums' in record else fold_inventory(record)
    q=folded['q']
    assert len(folded['family_power_sums'])==2
    assert all(type(row['monochromatic']) is bool for row in folded['family_power_sums'])
    powers={row['monochromatic']:row['powers'] for row in folded['family_power_sums']}
    assert set(powers)=={False,True} and all(len(v)==7 and all(isinstance(x,int) for x in v) for v in powers.values())
    if 'records' in folded:
        compiler.validate_inventory({'q':q,'records':unfold_inventory(folded)})
        for mono in (False,True):
            assert powers[mono]==[sum(row['count']*row['theta']**j for row in folded['records'] if row['monochromatic']==mono)
                                  for j in range(7)]
    a,b=powers[True],powers[False]
    assert a[0]==(q-5)//4 and b[0]==3*(q-1)//4
    assert b[1]==3*a[1]-6
    assert a[2]+b[2]==(q-3)*(q+1)
    values=(a[1],a[2],a[3],b[3],a[4],b[4],a[6]+b[6])
    return {'q':q,'statistics':dict(zip(STAT_NAMES,values)),
            'normalization_contract':folded['normalization_contract'],
            'ordered_distinct_row_weight':q*(q-1),'normalized_parameter_count':q-2,
            'trust_boundary':'Exact sufficient statistics from a separately certified inventory; this object does not certify realizability.'}


def seven_stat_global_third_moment(record,n):
    start=time.perf_counter()
    record=record if 'statistics' in record else seven_statistics(record)
    q=record['q'];assert q>=13 and q%4==1 and 6<=n<=q
    assert record['ordered_distinct_row_weight']==q*(q-1) and record['normalized_parameter_count']==q-2
    stats=record['statistics'];assert set(stats)==set(STAT_NAMES) and all(isinstance(v,int) for v in stats.values())
    powers={True:[(q-5)//4,stats['mono_theta1'],stats['mono_theta2'],stats['mono_theta3'],stats['mono_theta4']],
            False:[3*(q-1)//4,3*stats['mono_theta1']-6,(q-3)*(q+1)-stats['mono_theta2'],stats['mixed_theta3'],stats['mixed_theta4']]}
    c6=compiler.inclusion_moment(universal_sixth_union(),q,n)
    histograms=compiler.repeated_histograms(q);plans=[]
    for mono in (False,True):
        edges=(1,1,1) if mono else (1,1,-1)
        sign=1 if mono else -1
        taus=compiler.admissible_tau_nodes(q,edges)[:6]
        plans.append((mono,len(histograms),[sign*t for t in taus]))
        histograms.extend(compiler.distinct_histogram(q,edges,t) for t in taus)
    raw,backend=compiler.union_coefficients_batch(histograms,6,max_union=min(n,18))
    values=[compiler.inclusion_moment(row['coefficients'],q,n) for row in raw]
    assert values[1]==values[2]
    all_equal=q*values[0];twice=3*q*(q-1)*values[1]
    normalized=c6*stats['total_theta6'];details=[]
    for mono,index,nodes in plans:
        reduced=[value-c6*theta**6 for value,theta in zip(values[index:index+len(nodes)],nodes)]
        poly=polynomial_from_values(nodes[:5],reduced[:5])
        if len(nodes)>5:assert sum((c*nodes[-1]**j for j,c in enumerate(poly)),Fraction())==reduced[-1]
        contribution=sum((c*powers[mono][j] for j,c in enumerate(poly)),Fraction())
        normalized+=contribution
        details.append({'monochromatic':mono,'theta_nodes':nodes,'quartic_coefficients':[encode(c) for c in poly],
                        'extra_node_checked':len(nodes)>5,'synthetic_nodes_claimed_realizable':False,
                        'polynomial_scope':'all theta by degree bound' if len(nodes)>=5 else 'all admissible theta values only',
                        'normalized_low_power_contribution':encode(contribution)})
    distinct=q*(q-1)*normalized;total=all_equal+twice+distinct
    return {'q':q,'n':n,'degree':6,'third_moment':encode(total),'third_moment_float':float(total),
            'all_equal_rows_contribution':encode(all_equal),'exactly_two_equal_rows_contribution':encode(twice),
            'distinct_rows_contribution':encode(distinct),'universal_theta6_coefficient':encode(c6),
            'theta5_coefficient':[0,1],'family_plans':details,'measured_integer_statistics':7,
            'coefficient_evaluations':len(histograms),'maximum_coefficient_evaluations':15,
            'backend_seconds':backend,'elapsed_seconds':time.perf_counter()-start,
            'maximum_intermediate_bits':max(x['maximum_intermediate_bits'] for x in raw),
            'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in
                (Path(__file__),Path(__file__).with_name('folded_moment.py'),compiler.SOURCE if hasattr(compiler,'SOURCE') else Path(compiler.__file__))}}
