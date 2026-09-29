#!/usr/bin/env python3
"""Scientific plots of exact verified finite residual widths."""
from hashlib import sha256
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def main():
    inputs = {}; data = {}
    for name in ['p65537_canonical', 'p65537_saved', 'p65537_saved_reordered', 'symbolic_N32_k641']:
        path = HERE/(name+'.json'); data[name] = json.loads(path.read_text())['graph']['widths']
        inputs[path.name] = sha256(path.read_bytes()).hexdigest()
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for name, label, style, color in [('p65537_saved', 'Saved order', '-', '#c64c32'),
                                      ('p65537_canonical', 'Canonical encoding', '-', '#256d9b'),
                                      ('p65537_saved_reordered', 'Saved encoding, reordered', '--', '#2d8a56')]:
        axs[0].plot(data[name], style, label=label, color=color, linewidth=2)
    axs[0].set_title('Same 65,537 scalar encodings'); axs[0].legend(loc='upper left', frameon=False, fontsize=9)
    axs[0].annotate('1,597', (10, 1597), xytext=(10, 2500), ha='center', color='#c64c32')
    axs[0].set_ylim(.8, 6000)
    axs[1].plot(data['symbolic_N32_k641'], color='#73539a', linewidth=2)
    axs[1].set_title('All 6,700,417 encodings; cofactor 641')
    axs[1].annotate('Peak width 2,565 = 4 × 641 + 1', (16, 2565), xytext=(16, 4600), ha='center')
    axs[1].set_ylim(.8, 9000)
    for ax in axs:
        ax.set_yscale('log'); ax.set_xlabel('Number of digits read'); ax.set_ylabel('Nonempty residual states')
        ax.grid(True, alpha=.18); ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Exact state requirements depend on the digit representation', fontsize=13)
    fig.savefig(HERE/'state_widths.png', dpi=180); fig.savefig(HERE/'state_widths.pdf'); plt.close(fig)
    out = {'status': 'produced', 'scope': 'Exact fixed-order layered widths; implicit rejecting state excluded. No timing or asymptotic claim.',
           'source_sha256': {'plot_widths.py': sha256(Path(__file__).read_bytes()).hexdigest()}, 'input_sha256': inputs}
    (HERE/'plot_metadata.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
