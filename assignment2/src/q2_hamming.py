"""Q2: Hamming distance between the two equal-length sequences in a FASTA file.

Usage:
    python q2_hamming.py <fasta_path>
"""
import sys

from fasta import read_fasta


def hamming(s, t):
    if len(s) != len(t):
        raise ValueError(f"sequences differ in length ({len(s)} vs {len(t)})")
    return sum(a != b for a, b in zip(s, t))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python q2_hamming.py <fasta_path>")
    (_, s), (_, t) = read_fasta(sys.argv[1])[:2]
    print(hamming(s, t))
