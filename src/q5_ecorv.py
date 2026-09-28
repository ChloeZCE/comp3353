"""Q5: EcoRV restriction-site case study on chr22.

Usage:
    python src/q5_ecorv.py <chr22_fasta_path> [--outdir results]

Runs all four sub-questions (a-d) and writes a fragment-length histogram to
<outdir>/q5d_fragment_length_hist.png.
"""
import sys
import os
import re
import argparse
import statistics

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from common import read_fasta_records, reverse_complement, count_overlapping

SITE = "GATATC"


def load_single_sequence(path):
    records = list(read_fasta_records(path))
    if len(records) != 1:
        raise ValueError(f"Expected a single-record FASTA, found {len(records)}")
    return records[0][1]


def part_a_exact_sites(seq):
    return count_overlapping(seq, SITE)


def part_b_reverse_complement():
    return reverse_complement(SITE)


def part_c_relaxed_sites(seq):
    """GAX/ATC where X is any of A/C/G/T (mismatch tolerated at position 3)."""
    seq_u = seq.upper()
    pattern = re.compile(r"(?=(GA[ACGT]ATC))")
    total = sum(1 for _ in pattern.finditer(seq_u))
    exact = part_a_exact_sites(seq)
    return total, exact, total - exact


def part_d_fragments(seq, outdir):
    seq_u = seq.upper()
    site_len = len(SITE)
    cut_offset = 3  # GAT / ATC -> cut after the 3rd base

    # Find all non-overlapping exact site positions (GATATC cannot self-overlap).
    positions = []
    start = 0
    while True:
        idx = seq_u.find(SITE, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + site_len

    cut_points = [p + cut_offset for p in positions]
    boundaries = [0] + cut_points + [len(seq_u)]

    fragment_lengths = []
    for i in range(len(boundaries) - 1):
        fragment = seq_u[boundaries[i]:boundaries[i + 1]]
        length_no_n = len(fragment) - fragment.count("N")
        fragment_lengths.append(length_no_n)

    n_fragments = len(fragment_lengths)
    median_length = statistics.median(fragment_lengths)

    os.makedirs(outdir, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positive_lengths = [l for l in fragment_lengths if l > 0]
    import numpy as np
    log_bins = np.logspace(0, np.log10(max(positive_lengths)), 50)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(positive_lengths, bins=log_bins, color="#4C72B0", edgecolor="white")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Fragment length (bp, N's excluded, log scale)")
    ax.set_ylabel("Number of fragments (log scale)")
    ax.set_title(f"EcoRV digest of chr22: fragment length distribution (n={n_fragments})")
    fig.tight_layout()
    out_path = os.path.join(outdir, "q5d_fragment_length_hist.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    return n_fragments, median_length, out_path, fragment_lengths


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_path")
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    seq = load_single_sequence(args.fasta_path)

    print("=== (a) Exact EcoRV sites (GATATC) ===")
    n_exact = part_a_exact_sites(seq)
    print("Count:", n_exact)

    print("\n=== (b) Reverse complement of GATATC ===")
    rc = part_b_reverse_complement()
    print("Reverse complement:", rc)
    print("Palindromic:", rc == SITE)

    print("\n=== (c) Star-activity relaxed sites (GAX/ATC) ===")
    total, exact, relaxed_only = part_c_relaxed_sites(seq)
    print("Total GA[N]ATC matches (includes exact):", total)
    print("Of which exact (X=T):", exact)
    print("Of which relaxed-only (X != T):", relaxed_only)

    print("\n=== (d) Fragment lengths after cutting at every exact site ===")
    n_fragments, median_length, hist_path, _ = part_d_fragments(seq, args.outdir)
    print("Number of fragments:", n_fragments)
    print("Median fragment length (N's excluded):", median_length)
    print("Histogram saved to:", hist_path)
