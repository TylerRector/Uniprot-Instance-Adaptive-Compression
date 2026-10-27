# Uniprot Instance Adaptive Compression

Reference-based compression of 100 UniProt proteins across 9 families, with
Lean 4 proofs that selected functional-site residues remain exact.
This does not guarantee translation across biochemical activity.

| Encoding | Bytes | Exact full sequences? |
| --- | ---: | :---: |
| ASCII + LZMA/XZ | **6,224** | Yes |
| Guarded archive + gzip | 6,828 | No; annotated sites exact |
| ASCII + gzip | 7,475 | Yes |
| Guarded archive | 8,107 | No; annotated sites exact |
| Lossless archive | 9,803 | Yes |
| 5-bit packing | 11,519 | Yes |

Guarded correction preserves **511/511 annotated site features** (505/511
without correction). Lossless LZMA remains the smallest encoding.

Python uses only the standard library. Run from the repository root:

```bash
python3 -m unittest discover -s tests -v
python3 src/experiment.py
python3 src/pool_archive.py
lean Lean/FunctionGuard.lean
```

`src/`: Python code · `Lean/`: proofs · `data/`: frozen UniProt records ·
`results/`: measured outputs · [Detailed results](docs/RESULTS.md). MIT license.
