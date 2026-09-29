#!/usr/bin/env python3
"""Pin artifacts and inspected input documents; check syntax and stored source hashes."""
from datetime import datetime,timezone
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import sys

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent


def digest(path):return sha256(path.read_bytes()).hexdigest()


def main():
    scripts=sorted(ROOT.rglob('*.py'))
    for path in scripts:compile(path.read_text(),str(path),'exec')
    source_checks=[]
    pairs=[('observables/results_'+s+'.json','observables/fiber_audit.py','script_sha256') for s in ['toy','scale','holdout']]
    pairs += [('observables/results_trades.json','observables/moment_trades.py','script_sha256'),
              ('observables/results_incidence.json','observables/incidence_refinement.py','script_sha256'),
              ('spectral/results.json','spectral/transport_microscope.py','source_sha256'),
              ('spectral/refinement.json','spectral/refine_transport.py','source_sha256'),
              ('spectral/robustness.json','spectral/robustness.py','source_sha256'),
              ('spectral/validation.json','spectral/verify_witnesses.py','source_sha256')]
    for result,source,key in pairs:
        value=json.loads((ROOT/result).read_text())[key]
        assert value==digest(ROOT/source),(result,source,'source hash mismatch')
        source_checks.append({'result':result,'source':source,'sha256':value})
    broken=[]
    for path in ROOT.rglob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if target.startswith(('https://','http://','#','mailto:')):continue
            target=target.split('#')[0]
            if not (path.parent/target).exists():broken.append([str(path.relative_to(ROOT)),target])
    assert not broken,broken
    log=ROOT/'verification.log'
    assert log.exists() and 'All requested lab commands passed.' in log.read_text()
    inputs=[PROJECT/'README.md',PROJECT/'research/parallel25-pass-summary-2026-09-05.md',
            PROJECT/'research/parallel22-all-orders-obstruction-2026-09-05.md',
            PROJECT/'research/parallel21-squarefree-moments-2026-09-05.md',
            Path('/Users/shawwalters/proximityprize/docs/kb/deltastar-sw1-nec-2026-09-05.md'),
            Path('/Users/shawwalters/proximityprize/docs/kb/NEW_MATHEMATICS-2026-06-15.md')]
    artifacts={str(p.relative_to(ROOT)):{'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest.json'}
    manifest={'status':'finite artifact checks passed; no conjecture, production or priority certification',
              'created_utc':datetime.now(timezone.utc).isoformat(),
              'root_runtime':{'python':sys.executable,'version':platform.python_version(),
                              'packages':{x:importlib.metadata.version(x) for x in ['numpy','scipy','sympy','matplotlib']}},
              'spectral_discovery_runtime':{'python':'/opt/homebrew/opt/python@3.14/bin/python3.14',
                                            'numpy':'2.4.4','scipy':'1.17.1',
                                            'verified_live_by_root':True},
              'compiled_python_files':len(scripts),'markdown_local_links_passed':True,
              'current_source_hash_checks':source_checks,
              'inspected_source_documents':[{'path':str(p),'sha256':digest(p)} for p in inputs],
              'artifacts':artifacts}
    (ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':'passed','python_files':len(scripts),'artifacts':len(artifacts),
                      'source_hash_checks':len(source_checks),'local_links':'passed'}))


if __name__=='__main__':main()
