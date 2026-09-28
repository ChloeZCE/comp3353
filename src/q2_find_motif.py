"""Q2: Count occurrences of a motif s as a substring of a target DNA sequence.

Usage:
    python src/q2_find_motif.py <fasta_path> <motif>

Counts overlapping, case-insensitive matches on the forward strand only.
Each FASTA record is searched independently (matches are not allowed to span
across two different records), and counts are summed across all records.
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from common import read_fasta_records, count_overlapping


def count_motif_in_file(path, motif):
    total = 0
    for _header, seq in read_fasta_records(path):
        total += count_overlapping(seq, motif)
    return total


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python q2_find_motif.py <fasta_path> <motif>")
    path, motif = sys.argv[1], sys.argv[2]
    print(count_motif_in_file(path, motif))
