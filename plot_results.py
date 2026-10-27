"""Export reproducible PNG/SVG figures and plotted CSV source data."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
CHARTS = ROOT / 'charts'


def read_csv(name):
    with (RESULTS / name).open(newline='') as file:
        return list(csv.DictReader(file))


def prepare():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titlesize': 13, 'axes.titleweight': 'bold',
        'axes.labelcolor': '#3b4455', 'xtick.color': '#4e5967', 'ytick.color': '#4e5967',
        'grid.color': '#dbe2ea', 'savefig.facecolor': 'white', 'figure.facecolor': 'white',
    })


def storage(ax, agg, pooled):
    items = [
        ('Standalone\n5-bit', agg['sum_packed_standalone_bytes']),
        ('Naive guarded\n+93 references', agg['shared_reference_guarded_total_bytes']),
        ('Naive lossless\n+93 references', agg['shared_reference_lossless_total_bytes']),
        ('Pooled guarded\n16 references', pooled['guarded']['bytes']),
        ('Pooled lossless\n15 references', pooled['lossless']['bytes']),
    ]
    labels = [x[0] for x in items]
    vals = [x[1] for x in items]
    colors = ['#6c7c96', '#e7a44b', '#e6c07a', '#258e81', '#5875af']
    bars = ax.bar(range(len(vals)), vals, color=colors, width=.71)
    ax.set_xticks(range(len(vals)), labels, fontsize=8.5)
    ax.tick_params(axis='x', length=0)
    ax.set_ylim(0, max(vals) * 1.23)
    ax.set_ylabel('Total stored bytes')
    ax.set_title('Reference costs: measured payload and pool bytes')
    ax.yaxis.grid(True, alpha=.6)
    ax.set_axisbelow(True)
    for bar, val in zip(bars, vals):
        ax.annotate(f'{val:,}', (bar.get_x() + bar.get_width()/2, val),
                    ha='center', va='bottom', fontsize=9, fontweight='bold', xytext=(0, 5), textcoords='offset points')


def feature(ax, agg):
    groups = ['Real UniProt homologs', 'Synthetic single-site perturbation']
    unguarded = [agg['reference_site_recall_macro'] * 100,
                 agg['synthetic_reference_site_recall_macro'] * 100]
    guarded = [100, 100]
    x = np.arange(2)
    w = .34
    bars1 = ax.bar(x - w/2, unguarded, w, label='Reference only', color='#ce8c48')
    bars2 = ax.bar(x + w/2, guarded, w, label='Feature guarded', color='#258e81')
    ax.set_xticks(x, groups, fontsize=9)
    ax.set_ylim(0, 112)
    ax.set_ylabel('Mean exact site-feature agreement (%)')
    ax.set_title('Preserving annotated residues ≠ testing function')
    ax.yaxis.grid(True, alpha=.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, loc='lower left')
    for bars in (bars1, bars2):
        for bar in bars:
            val = bar.get_height()
            ax.annotate(f'{val:.1f}%', (bar.get_x() + bar.get_width()/2, val),
                        ha='center', va='bottom', fontsize=9, xytext=(0, 4), textcoords='offset points')


def whole_corpus(ax, agg, pooled):
    items = [
        ('ASCII + LZMA\n(lossless)', agg['ascii_corpus_lzma_xz_bytes'], True),
        ('Pooled guard + gzip\n(site-exact, lossy)', pooled['guarded']['gzip_bytes'], False),
        ('ASCII + gzip\n(lossless)', agg['ascii_corpus_gzip_bytes'], True),
        ('Pooled guard\n(site-exact, lossy)', pooled['guarded']['bytes'], False),
        ('Pooled lossless\n(lossless)', pooled['lossless']['bytes'], True),
        ('Packed 5-bit\n(lossless)', agg['sum_packed_standalone_bytes'], True),
    ]
    with (RESULTS / 'whole_corpus_baselines.csv').open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['method', 'bytes', 'sequence_lossless'])
        writer.writeheader()
        writer.writerows({'method': label.replace('\n', ' '), 'bytes': size, 'sequence_lossless': exact}
                         for label, size, exact in items)
    labels = [a for a, _, _ in items]
    values = [v for _, v, _ in items]
    colors = ['#5d76a8' if lossless else '#258e81' for _, _, lossless in items]
    bars = ax.bar(range(len(values)), values, color=colors, width=.73)
    ax.set_title('Whole-corpus compression changes the ranking')
    ax.set_ylabel('Bytes for all 100 sequences')
    ax.set_xticks(range(len(labels)), labels, fontsize=8.5)
    ax.set_ylim(0, max(values) * 1.19)
    ax.yaxis.grid(True, alpha=.6)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax.annotate(f'{value:,}', (bar.get_x() + bar.get_width()/2, value),
                    ha='center', va='bottom', fontsize=9, fontweight='bold', xytext=(0, 5), textcoords='offset points')


def tradeoff(ax, rows, agg):
    names = ['reference_only', 'guarded', 'optional25', 'optional50', 'optional75', 'lossless']
    labels = ['No edits', 'Guard', '+25%', '+50%', '+75%', 'Lossless']
    data = []
    total_residues = sum(int(row['length_aa']) for row in rows)
    for name, label in zip(names, labels):
        num_matches = sum(float(row[f'{name}_identity']) * int(row['length_aa']) for row in rows)
        site_sum = sum(float(row[f'{name}_site_recall']) * int(row['site_features']) for row in rows)
        data.append({
            'arm': name, 'label': label, 'payload_bytes': agg['payload_bytes'][name],
            'sequence_identity_percent': round(100 * num_matches / total_residues, 6),
            'site_feature_recall_percent': round(100 * site_sum / agg['site_feature_annotations'], 6),
        })
    with (RESULTS / 'rate_distortion.csv').open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    x = [r['payload_bytes'] for r in data]
    y = [r['sequence_identity_percent'] for r in data]
    ax.plot(x, y, '-o', lw=2.5, color='#5875af', markersize=7)
    ax.set_title('Sequence fidelity versus conditional payload')
    ax.set_xlabel('Payload bytes; reference storage excluded')
    ax.set_ylabel('Residue agreement with target (%)')
    ax.yaxis.grid(True, alpha=.6)
    ax.set_axisbelow(True)
    ax.set_xlim(min(x) - 130, max(x) + 350)
    ax.set_ylim(min(y) - 1, 101)
    for j, (xx, yy, label) in enumerate(zip(x, y, labels)):
        ax.annotate(label, (xx, yy), xytext=(6, -17 if j % 2 == 0 else 10),
                    textcoords='offset points', fontsize=8.5)


def family(ax, family_rows):
    names = [r['family'].replace('_', ' ') for r in family_rows]
    ratios = [100 * int(r['natural_site_challenge_cases']) / int(r['proteins']) for r in family_rows]
    site_recall = [100 * float(r['site_recall_reference_weighted']) for r in family_rows]
    yy = np.arange(len(names))
    ax.barh(yy, ratios, color='#cc8b49', height=.6, label='Targets needing site corrections')
    ax.set_yticks(yy, names, fontsize=8.5)
    ax.set_xlim(0, max(max(ratios) + 15, 50))
    ax.set_xlabel('Fraction of proteins with real annotated-site mismatch (%)')
    ax.set_title('Functional-site changes are uncommon in these homologs')
    ax.xaxis.grid(True, alpha=.6)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    for i, (rate, recall) in enumerate(zip(ratios, site_recall)):
        ax.text(rate + 1.5, i, f'{rate:.0f}%  •  site match {recall:.1f}%', va='center', fontsize=8)


def save(fig, name, dpi=200):
    fig.savefig(CHARTS / f'{name}.png', dpi=dpi, bbox_inches='tight')
    fig.savefig(CHARTS / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)


def run():
    prepare()
    CHARTS.mkdir(parents=True, exist_ok=True)
    agg = json.loads((RESULTS / 'aggregate.json').read_text())
    pooled = json.loads((RESULTS / 'pooled_aggregate.json').read_text())
    rows = read_csv('per_protein.csv')
    families = read_csv('by_family.csv')

    fig, ax = plt.subplots(figsize=(10.5, 5.3), layout='constrained')
    storage(ax, agg, pooled)
    save(fig, 'storage_comparison')

    fig, ax = plt.subplots(figsize=(11, 5.3), layout='constrained')
    whole_corpus(ax, agg, pooled)
    save(fig, 'whole_corpus_baselines')

    fig, ax = plt.subplots(figsize=(8.8, 4.9), layout='constrained')
    feature(ax, agg)
    save(fig, 'site_preservation')

    fig, ax = plt.subplots(figsize=(9.4, 5.1), layout='constrained')
    tradeoff(ax, rows, agg)
    save(fig, 'rate_distortion')

    fig, ax = plt.subplots(figsize=(9.6, 5.3), layout='constrained')
    family(ax, families)
    save(fig, 'family_challenges')

    fig, axes = plt.subplots(2, 2, figsize=(16, 11.5), layout='constrained')
    storage(axes[0, 0], agg, pooled)
    feature(axes[0, 1], agg)
    tradeoff(axes[1, 0], rows, agg)
    family(axes[1, 1], families)
    fig.suptitle('UniProt 100-protein guarded compression — audited results', fontsize=19, fontweight='bold')
    save(fig, 'overview', dpi=125)
    print('Exported 6 figures in PNG and SVG, plus plotted rate and baseline CSV files')


if __name__ == '__main__':
    run()
