#!/usr/bin/env python3
"""Reproduce or verify the research lab without network or external services."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
VERIFY=[['novelty/replay_observable_certificates.py'],
        ['novelty/review_incidence.py'],['novelty/review_proximity.py'],
        ['spectral/verify_witnesses.py'],['proximity/verify_artifacts.py'],
        ['proximity/extension_field_probe.py']]
FULL=[['observables/fiber_audit.py','--suite','toy'],
      ['observables/fiber_audit.py','--suite','scale'],
      ['observables/fiber_audit.py','--suite','holdout'],
      ['observables/moment_trades.py'],['observables/incidence_refinement.py'],
      ['spectral/transport_microscope.py'],['spectral/refine_transport.py'],
      ['spectral/robustness.py'],['proximity/deformation_microscope.py'],
      ['proximity/stack_provenance.py']]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=['verify','full'],default='verify',nargs='?')
    args=ap.parse_args()
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1'
    for command in (FULL+VERIFY if args.mode=='full' else VERIFY):
        print('Running '+' '.join(command),flush=True)
        subprocess.run([sys.executable,str(ROOT/command[0]),*command[1:]],cwd=ROOT,env=env,check=True)
    print('All requested lab commands passed. This verifies finite artifacts, not either conjecture.')


if __name__=='__main__':main()
