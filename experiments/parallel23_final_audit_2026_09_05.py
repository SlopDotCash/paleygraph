#!/usr/bin/env python3
"""Reconcile the four pass23 lanes, reviews, sources and preserved prior work."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import ast
import json
import re

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/parallel23_pass_audit_2026_09_05.json'
def H(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def J(p): return json.loads((ROOT/p).read_text())
central=['README.md','research/frontier.md','research/source-audit.md',
         'research/checkpoint-2026-09-04.md','sources/manifest.json']
prior=J('results/parallel23_prior_state_2026_09_05.json')
preserved=[]
for p,h in prior['prior_sha256'].items():
    if p in central: continue
    assert H(p)==h, ('prior artifact changed',p)
    preserved.append(p)
old=json.loads(prior['prior_manifest_text'])
assert sha256(prior['prior_manifest_text'].encode()).hexdigest()==prior['prior_sha256']['sources/manifest.json']
manifest=J('sources/manifest.json')
assert len(old)==27 and len(manifest)==32 and manifest[:27]==old
source_checks=[]
for s in manifest[27:]:
    assert (ROOT/s['path']).stat().st_size==s['bytes']
    assert H(s['path'])==s['sha256']
    source_checks.append({'path':s['path'],'sha256':s['sha256']})
extra=J('results/parallel23_prize_additional_sources_2026_09_05.json')
assert extra['commit']=='e65197892890b8fd9b0dc05b8980273cf1d595cc'
for s in extra['files']:
    assert H(s['path'])==s['sha256']
    assert any(t.get('path')==s['path'] and t['sha256']==s['sha256'] for t in manifest[27:])

lanes={
 'classical':'results/parallel23_classical_upper_2026_09_05.json',
 'subgroup':'results/parallel23_subgroup_upper_2026_09_05.json',
 'spectral':'results/parallel23_spectral_coupling_2026_09_05.json',
 'prize':'results/parallel23_prize_bridge_2026_09_05.json',
}
input_checks=[]
for lane,p in lanes.items():
    d=J(p)
    for name,h in d['input_sha256'].items():
        assert H(name)==h,(p,name)
        input_checks.append({'lane':lane,'path':name,'sha256':h})

review_checks=[]
reviews=[]
for lane in lanes:
    p='research/parallel23-'+lane+'-independent-review-2026-09-05.md'
    text=(ROOT/p).read_text()
    found=[]
    for line in text.splitlines():
        cols=line.split('|')
        if len(cols)<4: continue
        h=re.fullmatch(r'\s*`?([0-9a-f]{64})`?\s*',cols[2])
        if not h: continue
        cell=cols[1].strip().strip('`')
        link=re.search(r'\]\(([^)]+)\)',cell)
        if link:
            name=str((ROOT/Path(p).parent/link.group(1)).resolve().relative_to(ROOT))
        else:
            name=cell
        assert (ROOT/name).is_file(),(p,name)
        assert H(name)==h.group(1),(p,name,'review hash mismatch')
        found.append(name)
        review_checks.append({'review':p,'path':name,'sha256':h.group(1)})
    assert len(found)>=3,(p,found)
    assert lanes[lane] in found,(p,'results not pinned')
    reviews.append(p)

classical=J(lanes['classical'])
assert len(classical['small_exhaustive_slices'])==19
for f in classical['small_exhaustive_slices']:
    c=f['variance_certificate']
    assert not c['sidon_certificate_applicable'] and c['sidon_bad_fraction_upper'] is None
for f in classical['actual_prime_slice_samples']:
    c=f['probability_certificate']
    assert c['sidon_certificate_applicable'] and c['sidon_bad_fraction_upper'] is not None
assert classical['counts']=={
 'actual_quartic_energy_checks':5,'balanced_population_moment_domination':95,
 'exact_insertion_identities':15042,'exhaustive_slice_sets':3551,
 'full_character_transform_identities':1780,'large_prime_probability_certificates':5,
 'relation_count_fields':4,'sampled_actual_sidon_sets':287,
 'sixth_variance_bound_checks':19,'slice_poincare_checks':38,
}

artifacts=set()
for folder in ['research','experiments','results']:
    artifacts.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).glob('parallel23*') if p.is_file() and p!=OUT)
artifacts.update(s['path'] for s in manifest[27:])
syntax=[]
for p in sorted(artifacts):
    if p.endswith('.py'):
        ast.parse((ROOT/p).read_text())
        syntax.append(p)
links=[];excluded=[]
for doc in sorted(p for p in artifacts if p.endswith('.md'))+central[:-1]:
    for line in (ROOT/doc).read_text().splitlines():
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',line):
            if target.startswith(('http://','https://','#','mailto:','codex://')): continue
            if target=='u' and 'RawHyp' in line and doc=='research/source-audit.md':
                excluded.append({'source':doc,'target':target,'reason':'Historical mathematical function evaluation.'})
                continue
            path=(ROOT/Path(doc).parent/target.split('#')[0]).resolve()
            assert path.exists() or path==OUT,(doc,target)
            links.append({'source':doc,'target':target})

audit={
 'status':'Passed four-lane partial-result and artifact audit. Full goal remains unproved.',
 'audited_at_utc':datetime.now(timezone.utc).isoformat(),
 'previous_turn_classification':'progress',
 'input_hash_checks':input_checks,'separate_agent_review_input_checks':review_checks,
 'source_hash_checks':source_checks,'prior_artifacts_preserved':preserved,
 'source_manifest_prior_entries_preserved':27,'source_manifest_entries':32,
 'classical_reporting_repair':{
    'out_of_range_slices':19,'in_range_sample_certificates':5,
    'verifier_rerun_by_root':True,'all_check_counts_unchanged':True,
    'review':'research/parallel23-classical-independent-review-2026-09-05.md',
    'scope':'Out-of-range conditional Sidon bounds now null with false applicability. Proof unchanged; reviewer reconstructs old hashes after reversing reporting changes.',
 },
 'verification_counts_by_lane':{k:J(v).get('counts',J(v).get('checks')) for k,v in lanes.items()},
 'syntax_checks':syntax,'local_link_checks':links,'mathematical_notation_excluded':excluded,
 'artifact_sha256':{p:H(p) for p in sorted(artifacts)},
 'central_file_sha256':{p:H(p) for p in central},
 'separate_agent_reviews':reviews,
 'review_limitations':[
    'Agent review is not human refereeing or formal proof verification.',
    'Classical upper bound controls almost all Sidon inputs, not the exceptional ones.',
    'Subgroup upper bound covers the product-balanced part, not total opposite-free energy.',
    'Spectral result imports a known moment limit; no independent Hecke computation or all-vector upper bound.',
    'Prize count identities do not estimate production maxima or constitute a certificate.',
 ],
 'goal_status':'active and unachieved',
 'open_obligations':[
    'Uniform exceptional-input control sufficient for full classical Paley.',
    'Unbalanced subgroup relations and uniform square-root cancellation including exceptional primes.',
    'Remaining spectral operator with nonzero coupling retained.',
    'Arbitrary-word remainder-fibre and ratio bounds meeting the actual official certificate.',
 ],
}
OUT.write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'audit':str(OUT),'input_hash_checks':len(input_checks),
 'review_hash_checks':len(review_checks),'prior_artifacts_preserved':len(preserved),
 'artifacts':len(artifacts),'syntax_checks':len(syntax),'local_links':len(links),
 'manifest_entries':len(manifest),'goal_status':audit['goal_status']},indent=2))

