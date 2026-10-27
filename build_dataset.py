"""Freeze a deterministic 100-target, same-length homolog benchmark.

Consumes downloaded, *full JSON* UniProt REST search responses from
uniprot_query_cache/*.json. No network needed if cache exists.  Re-downloading
UniProt later may produce a different cohort; data/uniprot_100.json is the
frozen dataset distributed with this repository.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path

from codec import ALPHABET

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'data' / 'raw_queries'
OUTPUT = ROOT / 'data' / 'uniprot_100.json'

SITE_TYPES = {'Active site', 'Binding site', 'Metal binding', 'Disulfide bond'}
FEATURE_TYPES = SITE_TYPES | {'Signal'}
SIZES = {
    'carbonic_anhydrase': 8,
    'cytochrome_c': 10,
    'ferredoxin': 4,
    'hemoglobin': 4,
    'lysozyme': 20,
    'myoglobin': 20,
    'ribonuclease': 12,
    'superoxide_dismutase': 18,
    'thioredoxin': 4,
}
QUERIES = {
    'carbonic_anhydrase': 'reviewed:true AND protein_name:"carbonic anhydrase"',
    'cytochrome_c': 'reviewed:true AND protein_name:"cytochrome c"',
    'ferredoxin': 'reviewed:true AND protein_name:ferredoxin',
    'hemoglobin': 'reviewed:true AND protein_name:hemoglobin',
    'lysozyme': 'reviewed:true AND protein_name:lysozyme',
    'myoglobin': 'reviewed:true AND protein_name:myoglobin',
    'ribonuclease': 'reviewed:true AND protein_name:ribonuclease',
    'superoxide_dismutase': 'reviewed:true AND protein_name:"superoxide dismutase"',
    'thioredoxin': 'reviewed:true AND protein_name:thioredoxin',
}


def normalized_features(entry: dict) -> list[dict]:
    length = entry['sequence']['length']
    out = []
    for f in entry.get('features', []):
        name = f.get('type')
        if name not in FEATURE_TYPES:
            continue
        loc = f.get('location', {})
        s, e = loc.get('start', {}), loc.get('end', {})
        start, end = s.get('value'), e.get('value')
        if (s.get('modifier') != 'EXACT' or e.get('modifier') != 'EXACT'
                or not isinstance(start, int) or not isinstance(end, int)
                or not 1 <= start <= end <= length):
            continue
        # Active/binding/metal sites are point-like; disulfide bonds are two
        # endpoint residues, and signal peptides are short full spans.
        if name in SITE_TYPES - {'Disulfide bond'} and end - start > 4:
            continue
        if name == 'Signal' and end - start > 39:
            continue
        out.append({
            'type': name, 'start': start, 'end': end,
            'description': f.get('description', ''),
            'evidence_codes': sorted({v.get('evidenceCode', '') for v in f.get('evidences', [])}),
        })
    return out


def feature_positions(feature: dict) -> set[int]:
    start, end = feature['start'] - 1, feature['end'] - 1
    if feature['type'] == 'Disulfide bond':
        return {start, end}
    return set(range(start, end + 1))


def protected_positions(features: list[dict], length: int) -> set[int]:
    selected: set[int] = set()
    for f in features:
        positions = feature_positions(f)
        selected |= positions
        if f['type'] in SITE_TYPES - {'Disulfide bond'}:
            for p in positions:
                selected.update(range(max(0, p - 2), min(length, p + 3)))
    return selected


def normalize(entry: dict) -> dict | None:
    seq = entry['sequence']['value']
    if len(seq) > 600 or len(seq) < 40 or not set(seq).issubset(ALPHABET):
        return None
    feats = normalized_features(entry)
    if not any(f['type'] in SITE_TYPES for f in feats):
        return None
    return {
        'accession': entry['primaryAccession'],
        'name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
        'organism': entry.get('organism', {}).get('scientificName', ''),
        'sequence': seq,
        'sequence_version': entry.get('entryAudit', {}).get('sequenceVersion'),
        'features': feats,
    }


def reference_for(target: dict, pool: list[dict]) -> tuple[dict, float, int] | None:
    seq = target['sequence']
    sites = set().union(*(feature_positions(f) for f in target['features'] if f['type'] in SITE_TYPES))
    candidates = []
    for reference in pool:
        r = reference['sequence']
        if target['accession'] == reference['accession'] or len(r) != len(seq) or r == seq:
            continue
        matches = sum(a == b for a, b in zip(seq, r))
        identity = matches / len(seq)
        if not .70 <= identity <= .97:
            continue
        corrupted_sites = sum(seq[p] != r[p] for p in sites)
        candidates.append((matches, reference['accession'], reference, identity, corrupted_sites))
    if not candidates:
        return None
    if target['accession'] == 'P00698':
        turkey = [c for c in candidates if c[1] == 'P00703']
        if turkey:
            best = turkey[0]
            return best[2], best[3], best[4]
    # Nearest within a fixed nontrivial-divergence band, alphabetic tie break.
    best = sorted(candidates, key=lambda c: (-c[0], c[1]))[0]
    return best[2], best[3], best[4]


def build_dataset(cache_dir: Path = CACHE) -> dict:
    entries: dict[str, dict] = {}
    cases = []
    for family, quota in SIZES.items():
        response = json.loads((cache_dir / f'{family}.json').read_text())
        pool = [d for e in response['results'] if (d := normalize(e)) is not None]
        choices = []
        for entry in pool:
            ref = reference_for(entry, pool)
            if ref is not None:
                choices.append((entry, *ref))
        print(f'{family}: {len(choices)} candidate target accessions / {len({c[0]["sequence"] for c in choices})} distinct sequences')
        # Frozen deterministic accession-hash sampling, deduping target sequences.
        choices.sort(key=lambda c: sha256(('protein-guard-v1/' + c[0]['accession']).encode()).hexdigest())
        if family == 'lysozyme':
            choices.sort(key=lambda c: (c[0]['accession'] != 'P00698', sha256(('protein-guard-v1/' + c[0]['accession']).encode()).hexdigest()))
        sampled, used_sequences = [], set()
        for target, ref, identity, corrupted in choices:
            if target['sequence'] in used_sequences:
                continue
            sampled.append((target, ref, identity, corrupted))
            used_sequences.add(target['sequence'])
            if len(sampled) == quota:
                break
        if len(sampled) != quota:
            print(f'WARNING insufficient {family}: {len(sampled)} < {quota}')
        for target, ref, identity, corrupted in sampled:
            entries[target['accession']] = target
            entries[ref['accession']] = ref
            cases.append({
                'family': family,
                'target': target['accession'],
                'reference': ref['accession'],
                'reference_identity': round(identity, 8),
                'site_mismatch_positions': corrupted,
            })
    if len(cases) != 100 or len({x['target'] for x in cases}) != 100:
        raise RuntimeError('Expected 100 unique target accessions')
    return {
        'dataset': 'UniProt reviewed feature-challenge cohort, 100 homolog targets',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'source': 'https://rest.uniprot.org/uniprotkb/search',
        'query_format': 'format=json&size=350',
        'queries': QUERIES,
        'family_quotas': SIZES,
        'selection': 'SHA256-seeded distinct sequences; same length; 70-97% identity; nearest qualifying reviewed family entry',
        'entries': dict(sorted(entries.items())),
        'cases': cases,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', type=Path, default=CACHE)
    args = parser.parse_args()
    result = build_dataset(args.cache_dir)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print('Saved', OUTPUT, 'targets=', len(result['cases']), 'distinct entries=', len(result['entries']))
    print('Family counts:', dict(sorted(Counter(c['family'] for c in result['cases']).items())))
