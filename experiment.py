"""Run an auditable, actual-byte 100-protein UniProt sequence experiment.

Usage: python experiment.py
All results are computed from bundled frozen UniProt sequence/feature snapshots.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import lzma
from pathlib import Path
import statistics

from build_dataset import SITE_TYPES, feature_positions, protected_positions
from codec import ALPHABET, decode, decode_standalone, encode, encode_standalone

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data' / 'uniprot_100.json'
RESULTS = ROOT / 'results'


def exact_feature_matches(candidate: str, target: str, features: list[dict]) -> int:
    return sum(all(candidate[i] == target[i] for i in feature_positions(feature)) for feature in features)


def identity(a: str, b: str) -> float:
    assert len(a) == len(b)
    return sum(x == y for x, y in zip(a, b)) / len(a)


def stats(values: list[float | int]) -> dict:
    values = sorted(values)
    return {
        'mean': round(statistics.mean(values), 6),
        'median': round(statistics.median(values), 6),
        'min': min(values), 'max': max(values),
        'q1': round(values[(len(values) - 1) // 4], 6),
        'q3': round(values[(3 * (len(values) - 1)) // 4], 6),
    }


def mutate_for_stress(seq: str, position: int, target_aa: str) -> str:
    # This is SYNTHETIC, not a UniProt record or evidence of actual function.
    mutation = next(c for c in ALPHABET if c not in (seq[position], target_aa))
    return seq[:position] + mutation + seq[position + 1 :]


def run() -> tuple[list[dict], dict]:
    data = json.loads(DATA.read_text())
    rows = []
    reference_ids = set()
    all_arms = ['reference_only', 'guarded', 'optional25', 'optional50', 'optional75', 'lossless']
    payloads = {name: 0 for name in all_arms}
    stress_cases = []

    for case in data['cases']:
        target = data['entries'][case['target']]
        ref = data['entries'][case['reference']]
        seq, rseq = target['sequence'], ref['sequence']
        features = [f for f in target['features'] if f['type'] in SITE_TYPES]
        signals = [f for f in target['features'] if f['type'] == 'Signal']
        assert features
        protected = protected_positions(target['features'], len(seq))
        sites = set().union(*(feature_positions(f) for f in features))
        raw_differences = sum(a != b for a, b in zip(seq, rseq))
        guarded_differences = sum(seq[i] != rseq[i] for i in protected)
        ref_site_matches = exact_feature_matches(rseq, seq, features)
        reference_ids.add(case['reference'])

        outputs = {}
        blobs = {}
        recipes = {
            'reference_only': (set(), 0.0),
            'guarded': (protected, 0.0),
            'optional25': (protected, 0.25),
            'optional50': (protected, 0.50),
            'optional75': (protected, 0.75),
            'lossless': (protected, 1.0),
        }
        for label, (mask, fraction) in recipes.items():
            blob = encode(seq, rseq, mask, fraction)
            result = decode(blob, rseq)
            outputs[label] = result
            blobs[label] = blob
            payloads[label] += len(blob)
            assert all(result[p] == seq[p] for p in mask)
            if label != 'reference_only':
                assert exact_feature_matches(result, seq, features) == len(features)

        assert outputs['reference_only'] == rseq
        assert outputs['lossless'] == seq
        standalone = encode_standalone(seq)
        assert decode_standalone(standalone) == seq
        packed_ref = encode_standalone(rseq)
        assert decode_standalone(packed_ref) == rseq

        case_row = {
            'family': case['family'], 'target': case['target'], 'reference': case['reference'],
            'target_organism': target['organism'], 'reference_organism': ref['organism'],
            'length_aa': len(seq),
            'reference_identity': identity(seq, rseq),
            'mismatches': raw_differences,
            'protected_residues': len(protected),
            'protected_fraction': len(protected) / len(seq),
            'protected_mismatches': guarded_differences,
            'site_residues': len(sites),
            'site_mismatches': sum(seq[p] != rseq[p] for p in sites),
            'site_features': len(features),
            'experimentally_evidenced_sites': sum('ECO:0000269' in f['evidence_codes'] for f in features),
            'signal_features': len(signals),
            'reference_site_matches': ref_site_matches,
            'reference_site_recall': ref_site_matches / len(features),
            'guarded_site_recall': exact_feature_matches(outputs['guarded'], seq, features) / len(features),
            'reference_signal_matches': exact_feature_matches(rseq, seq, signals),
            'guarded_signal_matches': exact_feature_matches(outputs['guarded'], seq, signals),
            'ascii_bytes': len(seq), 'gzip_bytes': len(gzip.compress(seq.encode(), mtime=0)),
            'standalone_5bit_bytes': len(standalone),
            'reference_5bit_bytes': len(packed_ref),
            'feature_mask_bitmap_bytes': (len(seq) + 7) // 8,
        }
        for label in all_arms:
            case_row[f'{label}_bytes'] = len(blobs[label])
            case_row[f'{label}_identity'] = identity(seq, outputs[label])
            case_row[f'{label}_site_recall'] = exact_feature_matches(outputs[label], seq, features) / len(features)
            case_row[f'{label}_patches'] = sum(a != b for a, b in zip(outputs[label], rseq))

        # Controlled synthetic test explicitly isolated from the genuine data.
        unchanged_sites = sorted(p for p in sites if seq[p] == rseq[p])
        site_to_corrupt = unchanged_sites[0] if unchanged_sites else min(sites)
        altered_reference = mutate_for_stress(rseq, site_to_corrupt, seq[site_to_corrupt])
        altered_codec = encode(seq, altered_reference, protected)
        altered_output = decode(altered_codec, altered_reference)
        assert all(altered_output[i] == seq[i] for i in protected)
        assert exact_feature_matches(altered_output, seq, features) == len(features)
        stress_cases.append({
            'target': case['target'], 'position_1based': site_to_corrupt + 1,
            'synthetic_reference_site_recall': exact_feature_matches(altered_reference, seq, features) / len(features),
            'synthetic_guarded_site_recall': exact_feature_matches(altered_output, seq, features) / len(features),
            'synthetic_guarded_bytes': len(altered_codec),
        })
        rows.append(case_row)

    assert len(rows) == 100
    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / 'per_protein.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (RESULTS / 'synthetic_site_stress.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(stress_cases[0]))
        writer.writeheader()
        writer.writerows(stress_cases)

    grouped_rows = []
    for family in sorted({r['family'] for r in rows}):
        subs = [r for r in rows if r['family'] == family]
        sites_total = sum(r['site_features'] for r in subs)
        grouped_rows.append({
            'family': family, 'proteins': len(subs), 'mean_length_aa': round(statistics.mean(r['length_aa'] for r in subs), 2),
            'mean_reference_identity': round(statistics.mean(r['reference_identity'] for r in subs), 5),
            'total_site_features': sites_total,
            'natural_site_challenge_cases': sum(r['site_mismatches'] > 0 for r in subs),
            'site_recall_reference_weighted': round(sum(r['reference_site_matches'] for r in subs) / sites_total, 5),
            'site_recall_guarded_weighted': 1.0,
            'total_packed_bytes': sum(r['standalone_5bit_bytes'] for r in subs),
            'total_guarded_bytes': sum(r['guarded_bytes'] for r in subs),
            'total_lossless_bytes': sum(r['lossless_bytes'] for r in subs),
            'reference_pool_packed_bytes': sum(len(encode_standalone(data['entries'][i]['sequence'])) for i in {r['reference'] for r in subs}),
        })
    with (RESULTS / 'by_family.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(grouped_rows[0]))
        writer.writeheader()
        writer.writerows(grouped_rows)

    packed_total = sum(r['standalone_5bit_bytes'] for r in rows)
    ascii_corpus = ('\n'.join(data['entries'][case['target']]['sequence'] for case in data['cases']) + '\n').encode('ascii')
    packed_corpus = b''.join(encode_standalone(data['entries'][case['target']]['sequence']) for case in data['cases'])
    ref_pool_bytes = sum(len(encode_standalone(data['entries'][i]['sequence'])) for i in reference_ids)
    site_features = sum(r['site_features'] for r in rows)
    site_correct = sum(r['reference_site_matches'] for r in rows)
    aggregate = {
        'dataset': DATA.name,
        'dataset_sha256': hashlib.sha256(DATA.read_bytes()).hexdigest(),
        'proteins': len(rows), 'families': len(grouped_rows),
        'distinct_reference_accessions': len(reference_ids),
        'site_feature_annotations': site_features,
        'site_features_with_experimental_ECO_0000269': sum(r['experimentally_evidenced_sites'] for r in rows),
        'natural_challenged_proteins': sum(r['site_mismatches'] > 0 for r in rows),
        'natural_challenged_site_residues': sum(r['site_mismatches'] for r in rows),
        'reference_site_recall_micro': site_correct / site_features,
        'guarded_site_recall_micro': 1.0,
        'reference_site_recall_macro': statistics.mean(r['reference_site_recall'] for r in rows),
        'guarded_site_recall_macro': 1.0,
        'reference_identity': stats([r['reference_identity'] for r in rows]),
        'guarded_identity': stats([r['guarded_identity'] for r in rows]),
        'protected_fraction': stats([r['protected_fraction'] for r in rows]),
        'guarded_patches': stats([r['guarded_patches'] for r in rows]),
        'lossless_patches': stats([r['lossless_patches'] for r in rows]),
        'standalone_bytes': stats([r['standalone_5bit_bytes'] for r in rows]),
        'guarded_bytes': stats([r['guarded_bytes'] for r in rows]),
        'lossless_bytes': stats([r['lossless_bytes'] for r in rows]),
        'sum_ascii_bytes': sum(r['ascii_bytes'] for r in rows),
        'ascii_corpus_newline_delimited_bytes': len(ascii_corpus),
        'ascii_corpus_gzip_bytes': len(gzip.compress(ascii_corpus, mtime=0)),
        'ascii_corpus_lzma_xz_bytes': len(lzma.compress(ascii_corpus)),
        'packed_corpus_gzip_bytes': len(gzip.compress(packed_corpus, mtime=0)),
        'packed_corpus_lzma_xz_bytes': len(lzma.compress(packed_corpus)),
        'sum_gzip_bytes': sum(r['gzip_bytes'] for r in rows),
        'sum_packed_standalone_bytes': packed_total,
        'sum_unique_reference_packed_bytes': ref_pool_bytes,
        'sum_feature_mask_bitmap_bytes': sum(r['feature_mask_bitmap_bytes'] for r in rows),
        'payload_bytes': payloads,
        'shared_reference_guarded_total_bytes': payloads['guarded'] + ref_pool_bytes,
        'shared_reference_lossless_total_bytes': payloads['lossless'] + ref_pool_bytes,
        'shared_reference_guarded_vs_standalone_ratio': (payloads['guarded'] + ref_pool_bytes) / packed_total,
        'shared_reference_lossless_vs_standalone_ratio': (payloads['lossless'] + ref_pool_bytes) / packed_total,
        'synthetic_site_stress_proteins': len(stress_cases),
        'synthetic_reference_site_recall_macro': statistics.mean(r['synthetic_reference_site_recall'] for r in stress_cases),
        'synthetic_guarded_site_recall_macro': statistics.mean(r['synthetic_guarded_site_recall'] for r in stress_cases),
        'cost_caveats': [
            'Reference selection/accession index side information is not included in wire sizes.',
            'Functional-site annotation metadata and protected masks are not included in wire sizes.',
            'The real cohort is deterministic but not an independently sampled UniProt population.',
            'Exact sequence at annotated sites is not laboratory proof of biological function.',
            'Reference-conditioned streams require the exact correct reference at decode time.',
        ],
    }
    (RESULTS / 'aggregate.json').write_text(json.dumps(aggregate, indent=2, sort_keys=True) + '\n')
    return rows, aggregate


if __name__ == '__main__':
    records, totals = run()
    for key in ['proteins', 'families', 'site_feature_annotations', 'natural_challenged_proteins',
                'reference_site_recall_micro', 'sum_packed_standalone_bytes',
                'sum_unique_reference_packed_bytes', 'payload_bytes',
                'shared_reference_guarded_total_bytes', 'shared_reference_lossless_total_bytes']:
        print(f'{key}: {totals[key]}')
