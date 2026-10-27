"""Exact, reference-conditioned and annotation-guarded protein patch codec.

This is an intentionally simple baseline, not a learned biological compressor.
The small wire format is real and losslessly decodable *relative to its chosen
reference for every explicitly transmitted position*.
"""

from __future__ import annotations

import hashlib
import math
import struct
from typing import Collection

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_CODE = {aa: i for i, aa in enumerate(ALPHABET)}
MAGIC = b"FPR1"
HEADER = struct.Struct(">4sH8sH")  # 16 bytes; 0-based positional patches follow


def validate_sequence(sequence: str) -> None:
    if not (1 <= len(sequence) <= 65535) or any(aa not in AA_TO_CODE for aa in sequence):
        raise ValueError("Expected 1..65535 canonical amino acids (no ambiguity or gaps)")


def reference_digest(reference: str) -> bytes:
    validate_sequence(reference)
    return hashlib.blake2s(reference.encode("ascii"), digest_size=8).digest()


def index_bits(length: int) -> int:
    return max(1, (length - 1).bit_length())


def selected_patch_positions(
    target: str,
    reference: str,
    protected: Collection[int],
    optional_fraction: float = 0.0,
) -> list[int]:
    """Pick mandatory mismatches plus optional ones near protected positions.

    Coordinates are zero-based. Fractions use floor to avoid exceeding budget.
    Optional selection is a deterministic proximity heuristic, not a trained
    or independently validated functional-risk estimator.
    """
    validate_sequence(target)
    validate_sequence(reference)
    if len(target) != len(reference):
        raise ValueError("Reference and target must have equal lengths")
    if not math.isfinite(optional_fraction) or not (0 <= optional_fraction <= 1):
        raise ValueError("optional_fraction must be in [0,1]")
    sites = frozenset(protected)
    if any(not isinstance(i, int) or i < 0 or i >= len(target) for i in sites):
        raise ValueError("Protected position is out of bounds")
    mismatches = [i for i, (a, b) in enumerate(zip(target, reference)) if a != b]
    required = [i for i in mismatches if i in sites]
    optional = [i for i in mismatches if i not in sites]
    if sites:
        optional.sort(key=lambda i: (min(abs(i - j) for j in sites), i))
    quota = len(optional) if optional_fraction == 1 else math.floor(len(optional) * optional_fraction)
    return sorted(required + optional[:quota])


def encode(
    target: str,
    reference: str,
    protected: Collection[int] = (),
    optional_fraction: float = 0.0,
) -> bytes:
    positions = selected_patch_positions(target, reference, protected, optional_fraction)
    bit_width = index_bits(len(target))
    packed = 0
    for index in positions:
        packed = (packed << bit_width) | index
        packed = (packed << 5) | AA_TO_CODE[target[index]]
    bit_count = len(positions) * (bit_width + 5)
    body = packed.to_bytes((bit_count + 7) // 8, "big")
    return HEADER.pack(MAGIC, len(target), reference_digest(reference), len(positions)) + body


def decode(blob: bytes, reference: str) -> str:
    validate_sequence(reference)
    if len(blob) < HEADER.size:
        raise ValueError("Truncated header")
    magic, length, digest, count = HEADER.unpack(blob[: HEADER.size])
    if magic != MAGIC or length != len(reference) or digest != reference_digest(reference):
        raise ValueError("Wrong format, target length or reference fingerprint")
    if count > length:
        raise ValueError("Too many patches")
    width = index_bits(length)
    total_bits = count * (width + 5)
    body = blob[HEADER.size :]
    if len(body) != (total_bits + 7) // 8:
        raise ValueError("Wrong patch payload length")
    bits = int.from_bytes(body, "big")
    if bits >> total_bits:
        raise ValueError("Nonzero trailing padding")
    output = list(reference)
    previous = -1
    for j in range(count):
        shift = (count - 1 - j) * (width + 5)
        value = (bits >> shift) & ((1 << (width + 5)) - 1)
        index, symbol = value >> 5, value & 31
        if index <= previous or index >= length or symbol >= len(ALPHABET):
            raise ValueError("Malformed, unsorted or out-of-range patch")
        output[index] = ALPHABET[symbol]
        previous = index
    return "".join(output)


def encode_standalone(sequence: str) -> bytes:
    """Two-byte length + five bits per residue; no external reference."""
    validate_sequence(sequence)
    bits = 0
    for aa in sequence:
        bits = (bits << 5) | AA_TO_CODE[aa]
    return struct.pack(">H", len(sequence)) + bits.to_bytes((5 * len(sequence) + 7) // 8, "big")


def decode_standalone(blob: bytes) -> str:
    if len(blob) < 2:
        raise ValueError("Truncated sequence")
    (n,) = struct.unpack(">H", blob[:2])
    if n == 0 or len(blob) != 2 + (5 * n + 7) // 8:
        raise ValueError("Malformed sequence")
    bits = int.from_bytes(blob[2:], "big")
    if bits >> (5 * n):
        raise ValueError("Nonzero padding")
    result = []
    for i in range(n):
        symbol = (bits >> (5 * (n - i - 1))) & 31
        if symbol >= len(ALPHABET):
            raise ValueError("Invalid amino acid")
        result.append(ALPHABET[symbol])
    return "".join(result)
