#!/usr/bin/env python3
"""Bind the current round2 source/results and check portable artifact links."""
from datetime import datetime,timezone
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import sys

HERE=Path(__file__).resolve().parent


def digest(p):return sha256(p.read_bytes()).hexdigest()


def main():
    scripts=sorted(HERE.rglob('*.py'))
    for p in scripts:compile(p.read_text(),str(p),'exec')
    pairs=[('observables/results.json','observables/critical_audit.py','source_sha256'),
           ('observables/cohort_review.json','observables/cohort_review.py','script_sha256'),
           ('observables/exchange_spectrum_results.json','observables/exchange_spectrum.py','script_sha256'),
           ('observables/stress_exchange_results.json','observables/stress_exchange.py','script_sha256'),
           ('observables/review_exchange_results.json','observables/review_exchange.py','review_sha256'),
           ('observables/review_exchange_results.json','observables/exchange_spectrum.py','generator_sha256'),
           ('observables/exchange_grouped_results.json','observables/exchange_grouped.py','script_sha256'),
           ('observables/exchange_grouped_results.json','observables/review_exchange.py','review_backend_sha256'),
           ('novelty/exchange_grouped_review.json','observables/exchange_grouped.py','source_sha256'),
           ('novelty/exchange_grouped_review.json','novelty/review_exchange_grouped.py','review_sha256'),
           ('novelty/arithmetic_lift_results.json','novelty/arithmetic_lift.py','script_sha256'),
           ('novelty/arithmetic_lift_validation.json','novelty/verify_arithmetic_lift.py','verifier_sha256'),
           ('novelty/compressed_review.json','proximity/compressed_support.py','reviewed_source_sha256'),
           ('novelty/compressed_review.json','novelty/review_compressed.py','reviewer_sha256'),
           ('novelty/spectral_exact_review.json','spectral/exact_ntt.cpp','cpp_source_sha256'),
           ('novelty/spectral_exact_review.json','novelty/review_spectral_exact.py','reviewer_sha256'),
           ('novelty/million_samples_review.json','novelty/review_million_samples.py','reviewer_sha256'),
           ('spectral/exact_validation.json','spectral/exact_ntt.cpp','cpp_source_sha256'),
           ('spectral/exact_validation.json','spectral/verify_exact.py','python_source_sha256'),
           ('spectral/runtime.json','spectral/fft_transport.py','source_sha256'),
           ('spectral/extra_transport_check.json','spectral/extra_transport_check.py','source_sha256'),
           ('spectral/million_transport_check.json','spectral/million_transport_check.py','source_sha256')]
    checks=[]
    for result,source,key in pairs:
        record=json.loads((HERE/result).read_text());h=digest(HERE/source)
        assert record[key]==h,(result,source,'source/result hash mismatch')
        checks.append({'result':result,'source':source,'sha256':h})
    witness_count=0
    for p in (HERE/'spectral').glob('results_*.json'):
        record=json.loads(p.read_text())
        assert len(record['sides'])==2 and 'elapsed_seconds' in record
        for side in record['sides']:
            w=side['witness'];assert digest(HERE/'spectral'/w['file'])==w['sha256'];witness_count+=1
    for file in ['extra_transport_check.json','million_transport_check.json']:
        record=json.loads((HERE/'spectral'/file).read_text())
        assert digest(HERE/'spectral'/record['file'])==record['sha256']
        assert record['autocorrelation_computed'] and record['count_translation_ge_3_5']==1
        witness_count+=1
    assert witness_count==14
    broken=[]
    for p in HERE.rglob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https://','http://','#','mailto:')):continue
            dest=p.parent/target.split('#')[0]
            if dest==HERE/'manifest.json':continue
            if not dest.exists():broken.append((str(p.relative_to(HERE)),target))
    assert not broken,broken
    log=HERE/'verification.log'
    assert 'Round2 verification passed.' in log.read_text()
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p!=HERE/'manifest.json'}
    inputs=[HERE.parent/'observables/results_incidence.json',HERE.parent/'proximity/stack_results.json',
            HERE.parent/'proximity/extension_results.json',
            Path('/Users/shawwalters/proximityprize/scripts/probes/probe_syz39_sylvester_badprime_structure.py')]
    output={'created_utc':datetime.now(timezone.utc).isoformat(),
            'status':'current round2 source bindings and finite artifact checks passed; no conjecture or priority certification',
            'runtime':{'python':sys.executable,'version':platform.python_version(),
                       'packages':{x:importlib.metadata.version(x) for x in ['numpy','scipy','sympy','matplotlib']}},
            'python_syntax_checks':len(scripts),'source_bindings':checks,'witness_hashes':witness_count,
            'markdown_local_links_passed':True,
            'read_only_inputs':[{'path':str(p),'sha256':digest(p)} for p in inputs],
            'artifacts':artifacts,
            'note':'Nested lane manifests are earlier snapshots; this manifest pins current bytes after combined verification.'}
    (HERE/'manifest.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'passed':True,'python_files':len(scripts),'source_bindings':len(checks),
                      'witnesses':witness_count,'artifacts':len(artifacts)},indent=2))


if __name__=='__main__':main()
