# Results and interpretation — frozen 100-protein UniProt benchmark

## Scientific question and outcome

**Question:** Can a reference-conditioned sparse-correction code preserve
known protein functional-site amino acids while using fewer bytes than storing
protein sequences individually?

**Answer:** Yes for exact annotated-site residues and the simple 5-bit packing
baseline, if multiple proteins share references. But the better baseline is
whole-corpus generic compression: lossless LZMA wins against our lossy guarded
method on this cohort. No biological activity or model screening accuracy was
measured, so function preservation beyond these residue-level conditions is
**unproven**.

## Experimental design

1. Query nine protein-name families of **reviewed UniProtKB** entries using
   `format=json&size=350` from `https://rest.uniprot.org/uniprotkb/search`.
   Search queries and source accession metadata are in `data/uniprot_100.json`.
2. Keep canonical sequences 40–600 amino acids long with at least one exact
   `Active site`, `Binding site`, `Metal binding`, or `Disulfide bond` feature.
3. For each candidate, find a *different* reviewed family entry with the same
   length and 70–97% identity; use the nearest qualifying homolog. An explicit
   continuity anchor pairs chicken lysozyme P00698 with turkey P00703.
4. Deterministically sample 100 **distinct target sequences** by a documented
   SHA-256 accession ordering, with predetermined family quotas. This is not
   a representative random sample; nine family search pages are capped at 350.
5. Protect each exactly annotated site, both endpoints of disulfide bonds,
   the full annotated short signal peptide, and ±2 residues around active,
   ligand-binding and metal-binding sites. Transmit every mismatch in this mask;
   optionally transmit 25%, 50%, 75% or 100% of the other mismatches by a
   nearest-protected-position heuristic.
6. Decode actual bitpacked Python messages, assert exact site residues and
   record lengths, mutations, identities, payload byte counts and feature
   annotations. The lossless arm sends every mismatch.
7. Optimize **at most one reference per (family, exact length)** over the
   frozen target/reference candidates. Every target can instead store its own
   5-bit sequence. Write actual self-contained binary `.fpca` archives; decode
   them without external reference sequences.
8. Separately alter one annotated-site residue in each of the 100 homolog
   references **synthetically** and re-run the site guard. Those perturbations
   are not naturally occurring UniProt mutants or functional assays.

## Exact feature measurement

A feature is counted as preserved only if *all* of its exactly annotated
residue coordinates match the target sequence. Disulfide features are scored
on their two endpoint amino acids, **not** on whether the bond actually forms.

| Measure | Result |
| --- | ---: |
| Proteins / families | 100 / 9 |
| Site-feature annotations counted | 511 |
| Annotations marked with direct experimental evidence `ECO:0000269` | 35 |
| Natural reference-site matches | 505 / 511 (98.83%; micro) |
| After mandatory protected patches | **511 / 511 (100%; micro)** |
| Targets with any natural site-feature mismatch | **6 / 100** |
| Distinct natural mismatched site residues | 6 |
| Mean reference residue identity | 89.84% (per-protein macro) |
| Mean guarded residue identity | 90.88% (per-protein macro) |
| Mean fraction of sequence protected | 14.60% |
| Mean guarded patches per protein | 1.8 |
| Mean lossless patches per protein | 18.53 |

The feature-mask guarantee is automatic given the mask. What is informative
is the *cost* of enforcing it and the fact that naturally homologous sequences
already conserve almost all of these curated sites. Most entries' site
annotations are **not explicitly tagged `ECO:0000269`**: reviewed status must
not be mistaken for laboratory verification of every residue annotation.

In the explicitly artificial stress condition, mean exact feature agreement
falls to **65.59%** when one site residue per reference is deliberately
changed; guarded decoding restores it to **100%**. This validates the codec's
reaction to constructed mismatches, *not an assay of preserved activity*.

For lysozyme P00698 compared with P00703, residue 119 is a real target D vs
reference G difference in an annotated binding site. It is corrected by the
guard, agreeing with the smaller prior single-protein demonstration.

## Actual serialized bytes

| Encoding for all 100 target sequences | Bytes | Full sequence reconstructed? |
| --- | ---: | :---: |
| ASCII sequences with newline delimiters, LZMA/XZ | **6,224** | Yes |
| Pooled guarded archive, then gzip | **6,828** | No |
| ASCII sequences with newline delimiters, gzip | **7,475** | Yes |
| Pooled guarded `.fpca`, **16 references** | **8,107** | No |
| Pooled lossless `.fpca`, **15 references** | **9,803** | Yes |
| Five-bit standalone `.seq5` concatenation | **11,519** | Yes |
| Naive guarded payloads + 93 unique packed references | 12,662 | No |
| Naive lossless payloads + 93 unique packed references | 15,444 | Yes |

The *naive* two totals are sums of real serialized protein payloads and real
packed references, but they omit an additional reference-index mapping. They
are **optimistic** costs, not completed archives. The pooled `.fpca` values
are the on-disk lengths of complete decoder-tested binaries: 8-byte archive
header, 5-bit encoded reference pool, and 3- or 4-byte record headers plus
actual payloads. The exact files are under `results/`.

Measured differences:

- Pooled guarded saves **3,412 B / 29.62%** versus standalone 5-bit packing.
- Pooled lossless saves **1,716 B / 14.90%** versus standalone 5-bit packing.
- Pooled guarded + gzip saves **647 B / 8.66%** versus whole-corpus ASCII + gzip,
  but is **lossy** while the gzip baseline is not.
- Whole-corpus **lossless** ASCII + LZMA is **604 B smaller than guarded+gzip**,
  and **1,883 B smaller than uncompressed pooled guarded**.

All generic compression figures use Python's standard `gzip.compress(mtime=0)`
and `lzma.compress()` settings. The ASCII corpus has 100 newline delimiters;
the standalone binary contains explicit two-byte per-sequence lengths. Neither
format stores original UniProt accession strings or functional annotations.

The pooled archive records target order but **not** target accession names,
full UniProt feature JSON, protection masks or the algorithm's fitted choice
of codebook assignment as an external reusable model. Its reference sequences
and per-record reference indices **are** included. The encoder needs annotated
features to decide which amino acids to protect; the decoder only needs the
stored patches and stored reference pool. Including annotation metadata in the
storage objective would make both biological claims and compression ratios
considerably less attractive.

## Reference sharing drives the improvement

The first scheme chooses the nearest qualifying reference *per protein*: **93
unique reference accessions** for 100 targets, making aggregate storage worse
than packing directly. Pooled optimization chooses shared references inside
the 50 (protein family, sequence length) groups; it uses only **16** codebook
sequences for the guarded case, routes 66 targets through the codebook, and
stores 34 directly. The lossless optimizer uses 15 references for 59 targets,
and direct storage for 41.

This is restricted *in-sample* optimization of known target sequences. There
is no held-out inference model, learned functional codebook, trained
compression model, insertions/deletions, sequence alignment, structure
prediction, measured reaction rates, or external functional screening score.

## Machine-readable exports

`results/per_protein.csv`, `by_family.csv`, `rate_distortion.csv`,
`whole_corpus_baselines.csv`, `synthetic_site_stress.csv`, both pooled
assignment CSVs, `aggregate.json` and `pooled_aggregate.json` retain
reproducible values. The frozen cohort includes full source accessions, entry
sequence versions, descriptions and the annotated feature evidence codes.

## Formal assurance: narrow but real

`Lean/FunctionGuard.lean` was compiled with Lean 4.20.0, without any `sorry`
or `admit`. Its main guarantee is pointwise exactness for all selected sites,
regardless of optional edits. A second theorem states an override-only codec
must correct each differing target coordinate to be exact there. These proofs
do not establish bitpacking correctness (tested independently in Python),
biochemical function, folding, gene ontology labels, or downstream screening
performance. **Eleven Python tests passed**, including all 100 real cases,
full archive roundtrips, malformed-stream rejection and randomized edits.

## Follow-up questions for actual function screening research

1. Can *sequence-derived* saliency choose a mask without using the answers
   (target annotations) at encoding time?
2. How does preservation behave on function-relevant mutations with **assay
   data** (e.g., enzyme kinetics, deep mutational scans), not only site labels?
3. Does a learned, hierarchical family codebook beat strong generic corpus
   compressors while accounting for all model bits and decode cost?
4. What happens on remote homologs, insertions/deletions, unseen families,
   hetero-oligomeric assemblies and epistatic networks?
5. Does compression preserve function classifiers and calibrated confidence
   on held-out, sequence-identity-separated protein sets?

This benchmark supports a coordinate-preserving codec as a **baseline and
formal construction**. It does not yet establish a new training-scale,
function-preserving compression method.
