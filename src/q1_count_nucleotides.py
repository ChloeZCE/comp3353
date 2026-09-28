"""Q1: Count occurrences of A, C, G, T (case-insensitive) in a FASTA file.

Usage:
    python src/q1_count_nucleotides.py <fasta_path>

Reads the file line by line (no full-file load) so it scales to chromosome-sized
FASTA files such as data/Q1/chr22.fa.
"""
import sys
from collections import Counter


def count_nucleotides(path):
    counts = Counter()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            for ch in line.strip().upper():
                if ch in "ACGT":
                    counts[ch] += 1
    return counts["A"], counts["C"], counts["G"], counts["T"]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python q1_count_nucleotides.py <fasta_path>")
    a, c, g, t = count_nucleotides(sys.argv[1])
    print(a, c, g, t)
