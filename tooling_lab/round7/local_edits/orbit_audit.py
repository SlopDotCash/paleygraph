#!/usr/bin/env python3
"""Test whether apparent sufficiency merely identifies the full arithmetic orbit."""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import local_edits as tool

HERE=Path(__file__).resolve().parent


def canonical(c,q,squares):
    # A lexicographically least translate contains zero, so translate a member
    # of C to zero and minimize over every nonzero square multiplier.
    return min(tuple(sorted(a*(x-pivot)%q for x in c)) for a in squares for pivot in c)


def main():
    q=17;squares={x*x%q for x in range(1,q)};rows=[]
    for n in (6,7,8):
        groups=defaultdict(set);population=defaultdict(int);seen_orbits=set()
        for tail in combinations(range(2,q),n-2):
            c=(0,1)+tail;r=tool.from_prime(q,list(c))
            hist=tuple(map(tuple,r['signed_row_count_histogram']))
            summaries=tuple(sorted((x['insertion_sum'],x['insertion_norm_squared'],r['internal_contraction'][a][a],sum(r['internal_contraction'][a])) for a,x in enumerate(r['deletion_records'])))
            orbit=canonical(c,q,squares)
            groups[(hist,summaries)].add(orbit);population[(hist,summaries)]+=1;seen_orbits.add(orbit)
        # Every n-set must contain a square difference to reach the normalized
        # domain by a square affine map. Check this finite coverage explicitly.
        covered=0
        for c in combinations(range(q),n):
            assert any((a-b)%q in squares for a,b in combinations(c,2))
            covered+=1
        assert covered==comb(q,n)
        counts=defaultdict(int)
        for orbits in groups.values():counts[len(orbits)]+=1
        full_canonical={c:canonical(c,q,range(1,q)) for c in seen_orbits}
        full_counts=defaultdict(int)
        for orbits in groups.values():full_counts[len({full_canonical[c] for c in orbits})]+=1
        row={'q':q,'n':n,'normalized_inputs':sum(population.values()),'all_sets_normalization_coverage_checked':covered,
             'square_affine_group_size':q*len(squares),'square_affine_orbits':len(seen_orbits),
             'deletion_summary_fibres':len(groups),'fibres_by_number_of_distinct_orbits':dict(sorted(counts.items())),
             'all_summary_fibres_identify_one_full_orbit':all(len(x)==1 for x in groups.values()),
             'full_affine_orbits':len(set(full_canonical.values())),
             'fibres_by_number_of_full_affine_orbits':dict(sorted(full_counts.items())),
             'all_summary_fibres_contained_in_one_full_affine_orbit':all(len({full_canonical[c] for c in x})==1 for x in groups.values()),
             'canonical_orbit_representatives':[list(c) for c in sorted(seen_orbits)]}
        rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='canonical_orbit_representatives'}),flush=True)
    out={'status':'passed','scope':'Finite square-affine orbit saturation audit at Paley17. An orbit-identifying feature is sufficient for every invariant here, not evidence for a general cheap moment closure.',
         'cases':rows,'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('orbit_audit.py','local_edits.py')}}
    (HERE/'orbit_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
