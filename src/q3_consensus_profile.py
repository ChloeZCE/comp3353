"""Q3: Consensus string and profile matrix for aligned equal-length DNA sequences.

Variability metric: per-column Shannon entropy (bits) over the observed A/C/G/T
frequencies, H(i) = -sum_b p_b * log2(p_b). Entropy is used (rather than, say,
1 - max frequency) because it accounts for the *whole* distribution at a
column, not just the most common symbol -- a column split 13/12 between two
bases is clearly more variable than one split 24/1, but "count of the
non-majority symbols" alone does not distinguish them. Entropy is also the
natural building block for the mutual-information metric used in Q4, so using
it here keeps the two questions' variability measures consistent.
"""
import sys
import os
from math import log2

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from common import read_fasta_records

BASES = "ACGT"


def load_alignment(path):
    seqs = [seq.upper() for _header, seq in read_fasta_records(path)]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not all the same length: {lengths}")
    return seqs, lengths.pop()


def profile_matrix(seqs, length):
    profile = {b: [0] * length for b in BASES}
    for seq in seqs:
        for i, ch in enumerate(seq):
            if ch in profile:
                profile[ch][i] += 1
    return profile


def consensus_and_entropy(profile, length):
    consensus_chars = []
    entropies = []
    for i in range(length):
        counts = {b: profile[b][i] for b in BASES}
        total = sum(counts.values())
        best_base = max(BASES, key=lambda b: counts[b])
        consensus_chars.append(best_base if total > 0 else "N")

        entropy = 0.0
        if total > 0:
            for c in counts.values():
                if c > 0:
                    p = c / total
                    entropy -= p * log2(p)
        entropies.append(entropy)
    return "".join(consensus_chars), entropies


def top_variable_positions(entropies, n=3):
    order = sorted(range(len(entropies)), key=lambda i: (-entropies[i], i))
    return [i + 1 for i in order[:n]]  # 1-based


if __name__ == "__main__":
    seqs, length = load_alignment("data/Q3/BRCA_aligned.fa")
    profile = profile_matrix(seqs, length)
    consensus, entropies = consensus_and_entropy(profile, length)
    top3 = top_variable_positions(entropies, 3)

    print("Number of sequences:", len(seqs))
    print("Alignment length:", length)
    print("Consensus string:")
    print(consensus)
    for b in BASES:
        print(b + ":", " ".join(str(x) for x in profile[b]))
    print("Top 3 least-conserved (highest-entropy) 1-based positions:", top3)
    print("Their entropies (bits):", [round(entropies[i - 1], 4) for i in top3])
