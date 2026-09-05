#!/usr/bin/env python3
"""Two symmetry families for the existing exact third-moment compiler.

Known row-permutation/complement symmetry; no historical originality claim.
The scalar theta encodes the family modulo eight. Two families are retained
because a single polynomial need not represent both residue classes.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'third_moment'/'third_moment.py'
spec=importlib.util.spec_from_file_location('trace_invariants_compiler',SOURCE)
compiler=importlib.util.module_from_spec(spec);spec.loader.exec_module(compiler)


def encode(x):return [x.numerator,x.denominator]


def theta_family(q,theta):
    residue=theta%8
    if residue==(15-q)%8:return True
    assert residue==(q-7)%8, (q,theta)
    return False


def fold_inventory(record):
    q,rows=compiler.validate_inventory(record)
    folded=Counter()
    for row in rows:
        a,b,c=row['edges'];theta=a*b*c*row['tau']
        mono=a==b==c
        assert theta_family(q,theta)==mono
        folded[mono,theta]+=row['count']
    records=[{'monochromatic':m,'theta':t,'count':count} for (m,t),count in sorted(folded.items())]
    powers=[{'monochromatic':m,'powers':[sum(count*t**j for (mm,t),count in folded.items() if mm==m)
                                         for j in range(7)]} for m in (False,True)]
    return {'q':q,'records':records,'family_power_sums':powers,
            'normalized_parameter_count':q-2,'ordered_distinct_row_weight':q*(q-1),
            'normalization_contract':record.get('normalization',record.get('normalization_contract',
                'Requires independently verified normalized-pair transitivity under row/column permutations and complementation.')),
            'invariant':'theta=e01*e02*e12*tau; monochromatic iff theta mod8=(15-q) mod8',
            'input_inventory_rows':len(rows),'folded_inventory_rows':len(records)}


def unfold_inventory(folded):
    """Recover a balanced normalized inventory, given the transitivity contract.

    Row permutations and complementation distribute every mixed equivalence
    class equally across the three e01=+1 mixed patterns. This requires the
    same normalized-pair transitivity used by the global moment formula.
    """
    q=folded['q'];rows=[]
    for row in folded['records']:
        theta,count=row['theta'],row['count']
        assert theta_family(q,theta)==row['monochromatic']
        if row['monochromatic']:rows.append({'edges':[1,1,1],'tau':theta,'count':count})
        else:
            assert count%3==0
            for edges in ((1,-1,-1),(1,-1,1),(1,1,-1)):
                rows.append({'edges':list(edges),'tau':edges[0]*edges[1]*edges[2]*theta,'count':count//3})
    return sorted(rows,key=lambda row:(row['edges'],row['tau']))


def lagrange_polynomial(nodes,i):
    coeff=[Fraction(1)]
    for j,node in enumerate(nodes):
        if i==j:continue
        out=[Fraction() for _ in range(len(coeff)+1)]
        for k,v in enumerate(coeff):out[k]-=node*v;out[k+1]+=v
        coeff=[v/(nodes[i]-node) for v in out]
    return coeff


def polynomial_from_values(nodes,values):
    result=[Fraction() for _ in nodes]
    for i,value in enumerate(values):
        for j,c in enumerate(lagrange_polynomial(nodes,i)):result[j]+=value*c
    return result


def folded_global_third_moment(record,n,degree=6):
    start=time.perf_counter()
    folded=record if 'family_power_sums' in record else fold_inventory(record)
    q=folded['q'];assert q>=13 and degree in (0,2,4,6) and degree<=n<=q
    assert folded['normalized_parameter_count']==q-2
    assert folded['ordered_distinct_row_weight']==q*(q-1)
    assert len(folded['family_power_sums'])==2
    assert all(type(row['monochromatic']) is bool for row in folded['family_power_sums'])
    powers={row['monochromatic']:row['powers'] for row in folded['family_power_sums']}
    assert set(powers)=={False,True} and all(len(v)==7 and all(isinstance(x,int) for x in v) for v in powers.values())
    assert powers[True][0]==(q-5)//4 and powers[False][0]==3*(q-1)//4
    if 'records' in folded:
        # Validate the full support, not merely its residue, before using the
        # small-q finite-domain interpolation fallback.
        compiler.validate_inventory({'q':q,'records':unfold_inventory(folded)})
        assert all(theta_family(q,row['theta'])==row['monochromatic'] and row['count']>0 for row in folded['records'])
        for mono in (False,True):
            assert powers[mono]==[sum(row['count']*row['theta']**j for row in folded['records'] if row['monochromatic']==mono)
                                  for j in range(7)]
    histograms=compiler.repeated_histograms(q);plans=[]
    for mono in (False,True):
        edges=(1,1,1) if mono else (1,1,-1)
        sign=1 if mono else -1
        tau_nodes=compiler.admissible_tau_nodes(q,edges)[:degree+2]
        # When fewer than d+1 admissible values exist, ALL possible values
        # are used. The polynomial only needs to agree on this finite domain.
        index=len(histograms)
        histograms.extend(compiler.distinct_histogram(q,edges,tau) for tau in tau_nodes)
        plans.append((mono,index,[sign*t for t in tau_nodes]))
    raw,backend_seconds=compiler.union_coefficients_batch(histograms,degree,max_union=min(n,3*degree))
    moments=[compiler.inclusion_moment(row['coefficients'],q,n) for row in raw]
    assert moments[1]==moments[2]
    repeated=q*moments[0]+3*q*(q-1)*moments[1]
    distinct=Fraction();details=[]
    for mono,index,nodes in plans:
        values=moments[index:index+len(nodes)]
        basis=nodes[:degree+1];basis_values=values[:degree+1]
        poly=polynomial_from_values(basis,basis_values)
        extra=len(nodes)>degree+1
        if extra:assert sum((c*nodes[-1]**j for j,c in enumerate(poly)),Fraction())==values[-1]
        value=sum((c*powers[mono][j] for j,c in enumerate(poly)),Fraction())
        distinct+=q*(q-1)*value
        details.append({'monochromatic':mono,'theta_nodes':nodes,'polynomial_in_theta':[encode(c) for c in poly],
                        'extra_node_checked':extra,'synthetic_nodes_claimed_realizable':False,
                        'polynomial_scope':'all theta by degree bound' if len(basis)==degree+1 else 'all admissible theta values only',
                        'normalized_family_sum':encode(value)})
    total=repeated+distinct
    return {'q':q,'n':n,'degree':degree,'third_moment':encode(total),'third_moment_float':float(total),
            'all_equal_rows_contribution':encode(q*moments[0]),'exactly_two_equal_rows_contribution':encode(3*q*(q-1)*moments[1]),
            'distinct_rows_contribution':encode(distinct),'coefficient_evaluations':len(histograms),
            'maximum_generic_coefficient_evaluations':2*(degree+2)+3,'family_plans':details,
            'power_sum_input_integers':14,'parameter_rows_consumed_by_moment_backend':0,
            'maximum_intermediate_bits':max(row['maximum_intermediate_bits'] for row in raw),
            'backend_integer_products':sum(row['integer_products'] for row in raw),
            'backend_seconds':backend_seconds,'elapsed_seconds':time.perf_counter()-start,
            'source_sha256':{str(p.relative_to(HERE.parent)):sha256(p.read_bytes()).hexdigest() for p in
                (Path(__file__),SOURCE,SOURCE.parent/'union_coefficients.cpp',SOURCE.parent/'union_coefficients')}}
