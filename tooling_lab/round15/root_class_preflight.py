#!/usr/bin/env python3
"""A no-zero polynomial space whose cheapest unconstrained profile is impossible."""
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import sys
from projective_profile import optimize,pack

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round14'))
from arc_profile import profile


def main():
    sourcepath=HERE.parent/'round14/d3_s128.certificate.json';source=json.loads(sourcepath.read_text())['input'];p=source['p'];domain=source['domain'];n=len(domain);k=64;s=96
    roots=domain[:62];P=[1]
    for root in roots:
        new=[0]*(len(P)+1)
        for j,a in enumerate(P):new[j]=(new[j]-root*a)%p;new[j+1]=(new[j+1]+a)%p
        P=new
    basis=[[1],P,[0]+P];assert len(P)==63 and len(basis[2])==64
    # Independent product evaluation checks the coefficient expansion and
    # proves there are exactly62 roots on the actual domain.
    columns=[];groups=defaultdict(list)
    for i,x in enumerate(domain):
        product=1
        for root in roots:product=product*(x-root)%p
        value=sum(a*pow(x,j,p) for j,a in enumerate(P))%p;assert value==product
        assert (value==0)==(x in roots)
        v=(1,value,x*value%p);columns.append(v);groups[v].append(i)
        if value:assert v[2]*pow(v[1],-1,p)%p==x
    classes=list(groups.values());assert sorted(len(g) for g in classes)==[1]*962+[62]
    assert not any(not any(v) for v in columns)
    a,b,c=columns[0],columns[62],columns[63]
    rank_witness=(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))%p
    assert rank_witness
    naive=profile(n,s,3);constrained=optimize([len(g) for g in classes],0,s,3);best=constrained['best']
    assert naive['queries']==70070 and naive['queried_blocks']==47 and naive['unqueried_coordinates']==1
    assert 62>47+1 and best['queries']==136155 and best['queried_blocks']==33 and best['unqueried_coordinates']==29
    blocks=pack(classes,[],constrained);owner={i:j for j,B in enumerate(blocks) for i in B}
    assert sorted(owner)==list(range(n)) and len(owner)==sum(len(B) for B in blocks)
    for group in classes:
        assigned=[owner[i] for i in group if owner[i]<best['queried_blocks']]
        assert len(assigned)==len(set(assigned))
    assert sum(min(len(B),2) for B in blocks)==95<s and 64<96 and s*s<(k-1)*n
    data={'p':p,'n':n,'k':k,'s':s,'dimension':3,'domain':domain,'origin':[],'basis':basis}
    out={'status':'parallel_preflight_passed','scope':'Synthetic degree<64 rank-three polynomial space on the saved actual1024-point domain. Its unconstrained70070-query arc profile is impossible because of a62-point projective class. A136155-query relaxed profile has a valid class allocation; higher-rank query independence is not yet checked.',
         'input':data,'root_coordinates':list(range(62)),'projective_classes':classes,
         'full_rank_witness':{'coordinates':[0,62,63],'determinant_mod_p':rank_witness},
         'unconstrained_profile':naive,'parallel_constrained_profile':constrained,'allocated_blocks':blocks,
         'all_triple_ranks_checked':False,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['root_class_preflight.py','projective_profile.py','../round14/arc_profile.py']},
         'input_sha256':{'../round14/d3_s128.certificate.json':sha256(sourcepath.read_bytes()).hexdigest()}}
    (HERE/'root_class_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({'status':out['status'],'p':p,'n':n,'k':k,'s':s,'projective_classes':len(classes),
                      'largest_class':62,'unconstrained_queries':70070,'constrained_profile_queries':136155,'higher_rank_checks_pending':True}),flush=True)


if __name__=='__main__':main()
