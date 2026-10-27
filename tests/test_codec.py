"""Unit, randomized, full-data and actual archive round-trip tests."""

from __future__ import annotations

import json
import gzip
import lzma
from pathlib import Path
import random
import struct
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from build_dataset import SITE_TYPES, feature_positions, protected_positions
from codec import ALPHABET, decode, decode_standalone, encode, encode_standalone, selected_patch_positions
from pool_archive import decode_archive


class CodecTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dataset = json.loads((ROOT / 'data' / 'uniprot_100.json').read_text())

    def test_single_residue_alphabet_exhaustive(self):
        for target in ALPHABET:
            for ref in ALPHABET:
                self.assertEqual(decode(encode(target, ref, {0}), ref), target)
                self.assertEqual(decode(encode(target, ref, set()), ref), ref)
                self.assertEqual(decode_standalone(encode_standalone(target)), target)

    def test_randomized_guarded_roundtrips(self):
        rnd = random.Random(20261026)
        for n in [1, 2, 3, 15, 147, 257, 511]:
            for _ in range(12):
                seq = ''.join(rnd.choices(ALPHABET, k=n))
                ref = ''.join(rnd.choices(ALPHABET, k=n))
                guard = set(rnd.sample(range(n), rnd.randrange(n + 1)))
                for fraction in (0, .25, .5, .75, 1):
                    output = decode(encode(seq, ref, guard, fraction), ref)
                    self.assertTrue(all(output[i] == seq[i] for i in guard))
                    if fraction == 1:
                        self.assertEqual(output, seq)

    def test_patch_proximity_and_required_exact(self):
        seq, ref = 'ACDEFGHIK', 'AAAAAAAAA'
        sites = {2, 6}
        patches = selected_patch_positions(seq, ref, sites, .0)
        self.assertEqual(patches, [2, 6])
        self.assertEqual(selected_patch_positions(seq, ref, sites, 1), list(range(1, 9)))
        self.assertEqual(decode(encode(seq, ref, sites), ref)[2], seq[2])

    def test_wrong_reference_and_malformed(self):
        target, ref = 'ACDEFGHIK', 'AAAAAAAAA'
        blob = encode(target, ref, {2})
        with self.assertRaises(ValueError):
            decode(blob, 'CCCCCCCCC')
        with self.assertRaises(ValueError):
            decode(blob[:13], ref)
        with self.assertRaises(ValueError):
            decode(b'BAD!' + blob[4:], ref)
        with self.assertRaises(ValueError):
            decode(blob + b'\0', ref)
        with self.assertRaises(ValueError):
            encode('AAAA', 'AAA', {1})
        with self.assertRaises(ValueError):
            encode('AXAA', 'AAAA', {1})
        with self.assertRaises(ValueError):
            encode('AAAA', 'AAAA', {9})
        with self.assertRaises(ValueError):
            encode('AAAA', 'AAAA', {1}, 1.1)

    def test_noncanonical_and_padding_guard(self):
        with self.assertRaises(ValueError):
            encode_standalone('BB')
        with self.assertRaises(ValueError):
            decode_standalone(b'\0\x01\xff')
        with self.assertRaises(ValueError):
            decode_standalone(b'\0\x01')
        blob = bytearray(encode('CC', 'AA', {0}))
        # 1 index bit + 5 symbol bits => two high-order padding bits.
        blob[-1] |= 0b10000000
        with self.assertRaises(ValueError):
            decode(bytes(blob), 'AA')

    def test_dataset_has_100_reviewed_targets_and_9_families(self):
        cases = self.dataset['cases']
        self.assertEqual(len(cases), 100)
        self.assertEqual(len(set(c['target'] for c in cases)), 100)
        self.assertEqual(len(set(c['family'] for c in cases)), 9)
        self.assertIn('P00698', set(c['target'] for c in cases))
        entries = self.dataset['entries']
        for c in cases:
            target, ref = entries[c['target']], entries[c['reference']]
            self.assertEqual(len(target['sequence']), len(ref['sequence']))
            self.assertTrue(any(f['type'] in SITE_TYPES for f in target['features']))

    def test_every_case_feature_guarded_and_lossless(self):
        entries = self.dataset['entries']
        for c in self.dataset['cases']:
            target, ref = entries[c['target']], entries[c['reference']]
            seq, rseq = target['sequence'], ref['sequence']
            guard = protected_positions(target['features'], len(seq))
            guarded = decode(encode(seq, rseq, guard), rseq)
            full = decode(encode(seq, rseq, guard, 1), rseq)
            self.assertTrue(all(guarded[k] == seq[k] for k in guard))
            for f in target['features']:
                if f['type'] in SITE_TYPES:
                    self.assertTrue(all(guarded[k] == seq[k] for k in feature_positions(f)))
            self.assertEqual(full, seq)

    def test_lysozyme_binding_site_real_change(self):
        entries = self.dataset['entries']
        seq = entries['P00698']['sequence']
        ref = entries['P00703']['sequence']
        self.assertEqual((seq[118], ref[118]), ('D', 'G'))
        guard = protected_positions(entries['P00698']['features'], len(seq))
        self.assertIn(118, guard)
        self.assertEqual(decode(encode(seq, ref, guard), ref)[118], 'D')

    def test_actual_pooled_archives(self):
        entries = self.dataset['entries']
        for mode in ['guarded', 'lossless']:
            path = ROOT / 'results' / f'pooled_{mode}.fpca'
            blob = path.read_bytes()
            self.assertEqual(gzip.decompress((ROOT / 'results' / f'pooled_{mode}.fpca.gz').read_bytes()), blob)
            self.assertEqual(lzma.decompress((ROOT / 'results' / f'pooled_{mode}.fpca.xz').read_bytes()), blob)
            outputs = decode_archive(blob)
            self.assertEqual(len(outputs), 100)
            for seq, c in zip(outputs, self.dataset['cases']):
                target = entries[c['target']]
                guard = protected_positions(target['features'], len(target['sequence']))
                self.assertTrue(all(seq[i] == target['sequence'][i] for i in guard))
                if mode == 'lossless':
                    self.assertEqual(seq, target['sequence'])
            with self.assertRaises(ValueError):
                decode_archive(blob[:-1])
        direct = (ROOT / 'results' / 'standalone_100.seq5').stat().st_size
        self.assertLess((ROOT / 'results' / 'pooled_guarded.fpca').stat().st_size, direct)
        self.assertLess((ROOT / 'results' / 'pooled_lossless.fpca').stat().st_size, direct)

    def test_archive_header_validation(self):
        blob = (ROOT / 'results' / 'pooled_guarded.fpca').read_bytes()
        with self.assertRaises(ValueError):
            decode_archive(b'No archive')
        with self.assertRaises(ValueError):
            decode_archive(b'BAD!' + blob[4:])

    def test_whole_corpus_baseline_accounted(self):
        data = json.loads((ROOT / 'results' / 'aggregate.json').read_text())
        pool = json.loads((ROOT / 'results' / 'pooled_aggregate.json').read_text())
        self.assertEqual(data['site_feature_annotations'], 511)
        self.assertEqual(data['natural_challenged_proteins'], 6)
        self.assertEqual(pool['guarded']['bytes'], 8107)
        self.assertEqual(pool['lossless']['bytes'], 9803)
        self.assertLess(data['ascii_corpus_lzma_xz_bytes'], pool['guarded']['bytes'])


if __name__ == '__main__':
    unittest.main()
