"""A self-contained, actual-byte 100-protein archive with pooled references.

The restricted optimizer exhaustively tests *one reference per (family, length)*
from the frozen targets and their reviewed homologs. Each target may instead
use a standalone code. Results are in-sample corpus compression, NOT a test of
generalization or independent prediction. Target ordering/accession IDs and
functional annotation side information are external to both archive schemes.
"""

from __future__ import annotations

from collections import defaultdict
import csv
import gzip
import json
import lzma
from pathlib import Path
import struct

from build_dataset import SITE_TYPES, feature_positions, protected_positions
from codec import decode, decode_standalone, encode, encode_standalone

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'uniprot_100.json'
RESULTS = ROOT / 'results'
MAGIC = b'FPCA'
ARCHIVE_HEADER = struct.Struct('>4sHH')
RECORD_HEADER = struct.Struct('>cH')


def site_score(seq: str, target: str, features: list[dict]) -> int:
    return sum(all(seq[i] == target[i] for i in feature_positions(f)) for f in features)


def make_archive(dataset: dict, lossless: bool) -> tuple[bytes, list[dict], dict]:
    entries = dataset['entries']
    cases = dataset['cases']
    groups = defaultdict(list)
    for i, case in enumerate(cases):
        groups[(case['family'], len(entries[case['target']]['sequence']))].append(i)

    assignment = {}
    chosen_references = []
    selected_group_count = 0
    for group_key, members in sorted(groups.items()):
        references = {}
        for idx in members:
            case = cases[idx]
            for aid in (case['target'], case['reference']):
                seq = entries[aid]['sequence']
                references[seq] = min(aid, references.get(seq, aid))
        baseline_cost = sum(len(encode_standalone(entries[cases[i]['target']]['sequence'])) + 3 for i in members)
        best_cost = baseline_cost
        best_reference = None
        best_modes = None
        for rseq in sorted(references):
            cost = len(encode_standalone(rseq))
            proposed = {}
            for i in members:
                target = entries[cases[i]['target']]
                seq = target['sequence']
                protected = protected_positions(target['features'], len(seq))
                payload = encode(seq, rseq, protected, 1.0 if lossless else 0.0)
                packed = encode_standalone(seq)
                # Actual wire overhead: 1 mode + 2 length + 1 reference index.
                if len(payload) + 4 < len(packed) + 3:
                    proposed[i] = (b'R', payload)
                    cost += len(payload) + 4
                else:
                    proposed[i] = (b'S', packed)
                    cost += len(packed) + 3
            if cost < best_cost:
                best_cost, best_reference, best_modes = cost, rseq, proposed
        if best_reference is None:
            for i in members:
                assignment[i] = (b'S', None, encode_standalone(entries[cases[i]['target']]['sequence']))
        else:
            selected_group_count += 1
            ref_index = len(chosen_references)
            chosen_references.append((group_key, best_reference, references[best_reference]))
            for i, (mode, payload) in best_modes.items():
                assignment[i] = (mode, ref_index if mode == b'R' else None, payload)

    if len(chosen_references) > 255:
        raise ValueError('1-byte reference index cannot represent this pool')
    chunks = [ARCHIVE_HEADER.pack(MAGIC, len(chosen_references), len(cases))]
    chunks.extend(encode_standalone(seq) for _, seq, _ in chosen_references)
    output_rows = []
    for i, case in enumerate(cases):
        mode, index, payload = assignment[i]
        chunks.append(RECORD_HEADER.pack(mode, len(payload)))
        if mode == b'R':
            chunks.append(struct.pack('B', index))
        chunks.append(payload)
        output_rows.append({
            'family': case['family'], 'target': case['target'],
            'mode': mode.decode(),
            'reference_accession': chosen_references[index][2] if index is not None else '',
            'reference_pool_index': index if index is not None else '',
            'payload_bytes': len(payload),
            'wire_record_bytes': len(payload) + (4 if mode == b'R' else 3),
        })
    blob = b''.join(chunks)
    decoded = decode_archive(blob)
    if len(decoded) != len(cases):
        raise AssertionError('Archive decoder lost target records')
    for case, prediction in zip(cases, decoded):
        target = entries[case['target']]
        seq = target['sequence']
        protected = protected_positions(target['features'], len(seq))
        assert all(prediction[k] == seq[k] for k in protected)
        features = [f for f in target['features'] if f['type'] in SITE_TYPES]
        assert site_score(prediction, seq, features) == len(features)
        if lossless:
            assert prediction == seq
    return blob, output_rows, {
        'bytes': len(blob), 'reference_count': len(chosen_references),
        'reference_supported_targets': sum(r['mode'] == 'R' for r in output_rows),
        'standalone_targets': sum(r['mode'] == 'S' for r in output_rows),
        'group_count': len(groups), 'selected_reference_groups': selected_group_count,
        'all_feature_positions_preserved': True, 'all_sequences_lossless': lossless,
    }


def decode_archive(blob: bytes) -> list[str]:
    if len(blob) < ARCHIVE_HEADER.size:
        raise ValueError('Truncated archive')
    magic, nrefs, ntargets = ARCHIVE_HEADER.unpack_from(blob, 0)
    if magic != MAGIC:
        raise ValueError('Wrong archive type')
    offset = ARCHIVE_HEADER.size
    references = []
    for _ in range(nrefs):
        if offset + 2 > len(blob):
            raise ValueError('Truncated reference')
        (n,) = struct.unpack_from('>H', blob, offset)
        length = 2 + (5 * n + 7) // 8
        if offset + length > len(blob):
            raise ValueError('Truncated reference bytes')
        references.append(decode_standalone(blob[offset:offset + length]))
        offset += length
    predictions = []
    for _ in range(ntargets):
        if offset + RECORD_HEADER.size > len(blob):
            raise ValueError('Truncated protein record')
        mode, length = RECORD_HEADER.unpack_from(blob, offset)
        offset += RECORD_HEADER.size
        reference_index = None
        if mode == b'R':
            if offset >= len(blob):
                raise ValueError('Missing reference index')
            reference_index = blob[offset]
            offset += 1
            if reference_index >= nrefs:
                raise ValueError('Invalid reference index')
        elif mode != b'S':
            raise ValueError('Unknown record mode')
        if offset + length > len(blob):
            raise ValueError('Truncated protein payload')
        payload = blob[offset:offset + length]
        offset += length
        predictions.append(decode(payload, references[reference_index]) if mode == b'R' else decode_standalone(payload))
    if offset != len(blob):
        raise ValueError('Trailing archive bytes')
    return predictions


def run() -> dict:
    RESULTS.mkdir(parents=True, exist_ok=True)
    dataset = json.loads(DATA.read_text())
    packed_blob = b''.join(encode_standalone(dataset['entries'][case['target']]['sequence']) for case in dataset['cases'])
    (RESULTS / 'standalone_100.seq5').write_bytes(packed_blob)
    summary = {'standalone_100_bytes': len(packed_blob), 'records': len(dataset['cases'])}
    for label, lossless in [('guarded', False), ('lossless', True)]:
        blob, rows, metrics = make_archive(dataset, lossless)
        (RESULTS / f'pooled_{label}.fpca').write_bytes(blob)
        gz = gzip.compress(blob, mtime=0)
        xz = lzma.compress(blob)
        assert gzip.decompress(gz) == blob and lzma.decompress(xz) == blob
        (RESULTS / f'pooled_{label}.fpca.gz').write_bytes(gz)
        (RESULTS / f'pooled_{label}.fpca.xz').write_bytes(xz)
        with (RESULTS / f'pooled_{label}_assignments.csv').open('w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        summary[label] = metrics
        summary[label]['gzip_bytes'] = len(gz)
        summary[label]['lzma_xz_bytes'] = len(xz)
        summary[label]['ratio_to_standalone'] = len(blob) / len(packed_blob)
    (RESULTS / 'pooled_aggregate.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    return summary


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
