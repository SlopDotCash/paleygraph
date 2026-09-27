#!/usr/bin/env python3
"""Plot exact proved candidate caps, separately from observed list sizes."""
from hashlib import sha256
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def main():
    source = HERE/'pullback_summary.json'; heldout = HERE/'heldout_summary.json'
    rows = json.loads(source.read_text())['cases']; geometry = json.loads(heldout.read_text())['cases']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), constrained_layout=True)
    x = [r['erasures'] for r in rows]
    axes[0].plot(x, [r['baseline_universal_cap'] for r in rows], 'o-', color='#ba5039', label='Separate digit bounds')
    axes[0].plot(x, [r['universal_candidate_box_cap'] for r in rows], 'o-', color='#286f96', label='Centered pullback bounds')
    axes[0].set_title('Consecutive erasures: retain cancellations')
    axes[0].set_xlabel('Number of erased coordinates'); axes[0].set_ylim(.6, 8e11)
    axes[0].legend(loc='upper left', frameon=False, fontsize=9)
    axes[0].annotate('41,472', (32, 41472), xytext=(29, 250000), ha='center', color='#286f96')
    caps = [r['universal_candidate_box_cap'] for r in geometry]
    axes[1].bar(range(3), [float(v) for v in caps], color=['#286f96', '#aa7b4f', '#937494'])
    axes[1].set_xticks(range(3), ['Consecutive', 'Alternating', 'Seeded random'])
    axes[1].set_title('32 erasures: geometry still matters'); axes[1].set_ylim(1e3, 1e37)
    for j, value in enumerate(caps): axes[1].text(j, float(value)*6, f'{float(value):.2g}', ha='center', fontsize=9)
    for ax in axes:
        ax.set_yscale('log'); ax.set_ylabel('Uniform upper bound on candidate combinations')
        ax.grid(axis='y', alpha=.2); ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Certified erasure decoding at p = 2,013,265,921; N = 64', fontsize=13)
    fig.savefig(HERE/'candidate_caps.png', dpi=180); fig.savefig(HERE/'candidate_caps.pdf'); plt.close(fig)
    out = {'status': 'produced', 'scope': 'Exact universal candidate-box caps from rational certificates; not list-size observations, timings, or bounds for every erasure pattern.',
           'source_sha256': {'plot_caps.py': sha256(Path(__file__).read_bytes()).hexdigest()},
           'input_sha256': {source.name: sha256(source.read_bytes()).hexdigest(), heldout.name: sha256(heldout.read_bytes()).hexdigest()}}
    (HERE/'plot_metadata.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
