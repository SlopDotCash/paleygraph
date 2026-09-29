#!/usr/bin/env python3
"""Replay round2 or regenerate its bounded experiments. No external services."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
VERIFY=[['observables/cohort_review.py'],['observables/exchange_spectrum.py'],
        ['observables/stress_exchange.py'],['observables/review_exchange.py'],
        ['observables/exchange_grouped.py'],['novelty/review_exchange_grouped.py'],['novelty/verify_arithmetic_lift.py'],
        ['novelty/review_compressed.py'],['spectral/verify_exact.py'],
        ['spectral/extra_transport_check.py'],['spectral/million_transport_check.py'],
        ['novelty/review_spectral_exact.py'],['novelty/review_million_samples.py']]
EXPERIMENTS=[['observables/critical_audit.py'],['novelty/arithmetic_lift.py'],
             ['proximity/run_experiments.py'],['spectral/fft_transport.py'],
             ['spectral/completion_audit.py'],['spectral/scaled_operator_check.py']]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=['verify','full','million'],default='verify',nargs='?')
    args=ap.parse_args();env=os.environ.copy()
    env.update({'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    subprocess.run(['clang++','-std=c++17','-O3','-Wall','-Wextra',str(HERE/'spectral/exact_ntt.cpp'),
                    '-o',str(HERE/'spectral/exact_ntt')],check=True)
    if args.mode=='million':
        commands=[['spectral/fft_transport.py','--primes','1000033','--k','2','--tol','1e-6']]+VERIFY
    else:commands=(EXPERIMENTS if args.mode=='full' else [])+VERIFY
    for command in commands:
        print('Running '+' '.join(command),flush=True)
        subprocess.run([sys.executable,str(HERE/command[0]),*command[1:]],cwd=HERE,env=env,check=True)
    print('Round2 verification passed. Exact finite artifacts and declared statistical checks only.',flush=True)


if __name__=='__main__':main()
