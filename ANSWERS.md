# COMP3353 Bioinformatics — Assignment 1: Sequence Analysis

**Name:** Zheng Choi I
**University Number:** 3035987788
**Email:** u3598778@connect.hku.hk

All code below is implemented in Python 3 (using the standard library plus
`numpy`/`matplotlib` for Q4/Q5) and is also available, runnable, in this
submission's `src/` folder. All reported positions use **1-based indexing**.
`src/common.py` (reproduced once below) is a small shared helper module
imported by every question's script.

```python
"""Shared helpers for the COMP3353 Assignment 1 scripts."""

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

STANDARD_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta_records(path):
    """Yield (header, sequence) tuples without external deps, sequence as one string."""
    header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate_codon(codon):
    return STANDARD_CODON_TABLE.get(codon.upper(), "X")


def count_overlapping(seq, motif):
    """Count overlapping, case-insensitive occurrences of motif in seq."""
    seq_u = seq.upper()
    motif_u = motif.upper()
    count = 0
    start = 0
    while True:
        idx = seq_u.find(motif_u, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1
    return count
```

---

## Question 1: Counting DNA nucleotides

### Source code (`src/q1_count_nucleotides.py`)

```python
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
```

### Output

- On `data/Q1/input1.fa`: `13 13 17 17`  (A=13, C=13, G=17, T=17)
- On `data/Q1/chr22.fa`: `10382214 9160652 9246186 10370725`  (A=10382214, C=9160652, G=9246186, T=10370725)

---

## Question 2: Finding a motif in DNA

### Source code (`src/q2_find_motif.py`)

```python
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
```

Occurrences are counted **overlapping** (a sliding search that resumes one
position after each match start, so e.g. `AAA` in `AAAA` counts as 2), on the
forward strand only, matched case-insensitively.

### Output

- `{path: "data/Q1/input1.fa", s: "CGTAACC"}`: **4**
- `{path: "data/Q1/chr22.fa", s: "CGTAACC"}`: **206**

---

## Question 3: Consensus and profile matrix

### Variability metric

We measure variability at each alignment column using **Shannon entropy**
over the observed base frequencies:

H(i) = − Σ_b p_b · log2(p_b),  for b in {A, C, G, T}

where p_b is the fraction of the 25 aligned sequences carrying base b at
column i. We use entropy rather than a simpler count such as "number of
sequences differing from the majority base" because entropy captures the
*whole shape* of the column's distribution, not just the size of the
majority: a column split roughly 13/12 between two bases is intuitively far
more variable than one split 24/1, and entropy correctly rates the former as
close to its 2-bit maximum while a naive majority-mismatch count would not
distinguish them as sharply. Entropy is also the natural building block of
the mutual-information metric used in Question 4 (for two columns,
I(i,j) = H(i) + H(j) − H(i,j)), so using the same information-theoretic
quantity keeps the two questions' notions of "variability" consistent.

### Source code (`src/q3_consensus_profile.py`)

```python
"""Q3: Consensus string and profile matrix for aligned equal-length DNA sequences.

Usage:
    python src/q3_consensus_profile.py <fasta_path>

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
    if len(sys.argv) != 2:
        sys.exit("Usage: python q3_consensus_profile.py <fasta_path>")
    seqs, length = load_alignment(sys.argv[1])
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
```

### Output (on `data/Q3/BRCA_aligned.fa`, 25 sequences × 335 columns)

**Consensus string:**
```
GATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGAGATGTGGTCAATGGAAGAAACCACCAAGGTCCAAAGCGAGCAAGAGAATCCCAGGACAGAAAGGTAAAGCTCCCTCCCTCAAGTTGACAAAAATCTCACCCCACCACTCTGTATTCCACTCCCCTTTGCAGAGATGGGCCGCTTCATTTTGTAAGACTTATTACATACATACACAGTGCTAGATACTTTCACACAGGTTCTTTTTTCACTCTTCCATCCCAACCACATAAATAAGTATTGTCTCTACTTTATGAATGATAAAACTAAGAGATTTAGAGAGGCTGTGTA
```

**3 least-conserved (highest-entropy) 1-based positions:** 26, 158, 144
(entropies 1.164, 1.164, 1.0211 bits respectively — positions 26 and 158 tie
exactly and are both reported; ties are broken by lower column index first).
The full per-base profile matrix (counts of A/C/G/T at every column) is
computed by the script and available in `results/q3_output.txt`.

---

## Question 4: Multiple sequence alignment and covariation matrices

Alphabet: the 20 standard amino acids plus `-` for gaps; any non-standard
character is folded into the gap category, per the assignment's instructions.
Consensus at each column is simply the most frequent of these 21 symbols
(so a column can have `-` as its consensus if most sequences carry a gap
there). Mutual information is computed exactly as defined in the assignment,
skipping (a,b) terms with zero joint frequency (by the usual convention that
0·log2(0/x) = 0).

### Source code (`src/q4_covariation.py`)

```python
"""Q4: Consensus string and top-10 co-varying column pairs (mutual information)
for a protein multiple sequence alignment.

Usage:
    python src/q4_covariation.py <fasta_path>

Alphabet: the 20 standard amino acids plus '-' for gaps (21 symbols). Any
character outside the 20 standard amino acids (including '-') is treated as a
gap, per the assignment's instructions.

Mutual information for columns i, j:
    MI(i,j) = sum_a sum_b f_ij(a,b) * log2( f_ij(a,b) / (f_i(a) * f_j(b)) )
Terms where f_ij(a,b) == 0 are skipped (by convention, 0*log2(0/x) = 0).

Implementation note: rather than looping over all O(L^2) column pairs and all
21x21 symbol pairs in pure Python, we one-hot encode the alignment into a
(N_sequences, L*21) matrix M and compute M.T @ M in one matrix multiplication.
This single (L*21) x (L*21) matrix simultaneously gives the joint symbol
co-occurrence counts for every pair of columns (and, on its diagonal blocks,
each column's own marginal counts), which is far faster than a naive
Python triple-nested loop for alignments with hundreds of columns.
"""
import sys
import os
from math import log2

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from common import read_fasta_records

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"
ALPHABET = AMINO_ACIDS + GAP  # 21 symbols, gap last
SYMBOL_INDEX = {ch: i for i, ch in enumerate(ALPHABET)}
N_SYMBOLS = len(ALPHABET)


def load_alignment(path):
    seqs = [seq.upper() for _header, seq in read_fasta_records(path)]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not all the same length: {lengths}")
    return seqs, lengths.pop()


def encode(seqs, length):
    """Return integer code matrix (N, L) with any non-standard char -> gap index."""
    gap_idx = SYMBOL_INDEX[GAP]
    n = len(seqs)
    codes = np.full((n, length), gap_idx, dtype=np.int16)
    for r, seq in enumerate(seqs):
        for c, ch in enumerate(seq):
            codes[r, c] = SYMBOL_INDEX.get(ch, gap_idx)
    return codes


def consensus_string(codes, length):
    consensus = []
    for i in range(length):
        col = codes[:, i]
        counts = np.bincount(col, minlength=N_SYMBOLS)
        consensus.append(ALPHABET[int(np.argmax(counts))])
    return "".join(consensus)


def one_hot(codes, n, length):
    onehot = np.zeros((n, length, N_SYMBOLS), dtype=np.float64)
    rows = np.arange(n)[:, None]
    cols = np.arange(length)[None, :]
    onehot[rows, cols, codes] = 1.0
    return onehot.reshape(n, length * N_SYMBOLS)


def mutual_information_all_pairs(codes, n, length):
    M = one_hot(codes, n, length)          # (N, L*21)
    joint_counts = M.T @ M                 # (L*21, L*21)
    marginal_counts = joint_counts.diagonal().reshape(length, N_SYMBOLS)  # (L, 21)
    marginal_p = marginal_counts / n

    joint_counts = joint_counts.reshape(length, N_SYMBOLS, length, N_SYMBOLS)
    joint_p = joint_counts / n

    results = []
    for i in range(length):
        for j in range(i + 1, length):
            pij = joint_p[i, :, j, :]          # (21, 21)
            pi = marginal_p[i]                 # (21,)
            pj = marginal_p[j]                 # (21,)
            denom = np.outer(pi, pj)
            mask = pij > 0
            if not np.any(mask):
                continue
            mi = np.sum(pij[mask] * np.log2(pij[mask] / denom[mask]))
            results.append((mi, i, j))
    results.sort(key=lambda t: -t[0])
    return results


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python q4_covariation.py <fasta_path>")
    seqs, length = load_alignment(sys.argv[1])
    n = len(seqs)
    codes = encode(seqs, length)

    consensus = consensus_string(codes, length)
    print("Number of sequences:", n)
    print("Alignment length:", length)
    print("Consensus string:")
    print(consensus)

    mi_results = mutual_information_all_pairs(codes, n, length)
    print("Top 10 co-varying column pairs (1-based positions), by mutual information:")
    for mi, i, j in mi_results[:10]:
        print(f"  columns ({i + 1}, {j + 1}): MI = {mi:.4f} bits")
```

### Output (on `data/Q4/DHFR_aligned.fa`, 667 sequences × 160 columns)

**Consensus string:**
```
MISLIVAMAENGVIGKDNDLPWHLPEDLKYFKRLTLGKPVIMGRKTWESIGRPLPGRRNIVLSRDPDYQAEGVEVVHSLEEALALAG-VEEVFVIGGAEIYAQALP-ADRLYLTEIDAEFEGDTFFPEIDPDEWEEVSREEHPADEKNGYDYTFVTYERK
```

**Top 10 co-varying column pairs (1-based positions), by mutual information:**

| Rank | Columns (i, j) | MI (bits) |
|---|---|---|
| 1 | (149, 150) | 1.2409 |
| 2 | (150, 151) | 1.1677 |
| 3 | (68, 70)   | 1.1111 |
| 4 | (148, 150) | 1.0829 |
| 5 | (145, 146) | 1.0511 |
| 6 | (13, 122)  | 1.0165 |
| 7 | (149, 151) | 1.0131 |
| 8 | (6, 8)     | 1.0072 |
| 9 | (151, 152) | 0.9992 |
| 10 | (58, 74)  | 0.9862 |

Note several of the top pairs are adjacent/near-adjacent columns (e.g.
148–151), which is expected since neighbouring columns in a real protein
alignment often co-vary due to local structural constraints (e.g. a shared
loop or indel boundary) in addition to any long-range structural contacts.

---

## Question 5: A restriction-nuclease case study (chr22)

### Source code (`src/q5_ecorv.py`)

```python
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
```

### (a) Exact EcoRV sites

Count of exact, case-insensitive `GATATC` matches in `data/Q1/chr22.fa`:
**5024**

### (b) Reverse complement and why one strand suffices

The reverse complement of `GATATC` is **`GATATC`** itself (verified by the
script) — the site is palindromic.

Because DNA strands are antiparallel, the recognition site's presence at a
given location is a property of the *base pair*, not of a single strand: if
the top strand reads `GATATC` (5'→3') at some position, the bottom strand,
read in its own 5'→3' direction at that same physical location, reads the
reverse complement of `GATATC` — which, since the site is palindromic, is
again `GATATC`. So a palindromic site is automatically identical on both
strands at every occurrence; scanning the reverse-complement strand would
only ever rediscover the exact same genomic locations already found on the
forward strand (mirrored), never reveal a new cutting site. Hence counting
matches on the forward strand alone already accounts for every place the
double-stranded site (and therefore the enzyme's cut) occurs; counting both
strands would double-count identical sites rather than find additional ones.

### (c) Star-activity relaxed sites

Pattern used: `GA[ACGT]ATC` (mismatch tolerated only at position 3, restricted
to the four standard bases, case-insensitive), counted with the same
overlapping search as parts (a)/(b).

- Check string `TTGATATCAAGAGATCCTTCCGAAATCACGT`: 1 exact + 2 relaxed-only = 3
  total matches — matches the assignment's worked example exactly.
- On `data/Q1/chr22.fa`: **31370** total `GA[N]ATC` matches, of which
  **5024** are the exact site (X = T, matching part (a)) and **26346** are
  relaxed-only sites (X ≠ T) that only the star-activity mutant would cut.

### (d) Fragment lengths after digestion

Cutting chr22 at every exact `GATATC` match (blunt cut after the 3rd base,
i.e. `GAT/ATC`) and excluding all `N` characters when measuring each
fragment's length:

- **Number of fragments:** 5025
- **Median fragment length (N's excluded):** 3561 bp

**Fragment length distribution** (log–log histogram, `results/q5d_fragment_length_hist.png`):

![EcoRV fragment length distribution](results/q5d_fragment_length_hist.png)

The distribution is unimodal on a log scale, peaking around 1 kb and
right-skewed with a long tail out past 10^5 bp, consistent with a Poisson-like
process of randomly spaced 6-bp cut sites across a large chromosome.

---

## Question 6: Open reading frame (ORF) finder

### ORF definition and coordinate convention used

Within each of the 3 reading frames per strand, we scan codon by codon:
every `ATG` opens a candidate ORF; a run of codons between two in-frame stop
codons (or between a stop codon and the end of the frame) may contain several
`ATG`s, and each of them starts its own ORF that all close at the *same* next
in-frame stop codon (a stop always terminates every currently open ORF in
that frame). An `ATG` with no downstream in-frame stop codon before the end
of the genome is an incomplete ORF and is excluded, since it has no defined
end position or terminated protein.

All positions are reported in the **original genome's 1-based forward-strand
numbering**, regardless of which strand the ORF is actually on. For a
`+`-strand ORF, start < end. For a `-`-strand ORF, start > end: "start" is
the genomic position of the first base of the `ATG` as read 5'→3' on the
minus strand (the *higher* forward-strand coordinate), and "end" is the last
base of the stop codon (the *lower* forward-strand coordinate) — this keeps
"start"/"end" meaning "where translation begins/ends" on each ORF's own
strand, while all coordinates stay expressed in one consistent numbering
system.

### Source code (`src/q6_orf_finder.py`)

```python
"""Q6: ORF finder for the E. coli genome.

Usage:
    python src/q6_orf_finder.py <genome_fasta_path> [--outdir results]

Writes results/q6_orfs.csv with columns: strand,id,start,end,protein
(only ORFs whose translated protein is longer than 100 amino acids, on both
strands), and prints the single longest ORF across both strands.

ORF definition used: within a single reading frame, scan codon by codon.
Every ATG start codon begins a candidate ORF; if a run of codons (between two
in-frame stop codons, or between a stop codon and the end of the frame)
contains one or more ATGs, each of those ATGs starts its own ORF that all
share the *same* next in-frame stop codon (a stop always closes every
still-open ORF in that frame). An ATG with no in-frame stop codon before the
end of the genome is an incomplete ORF and is excluded, since it has no
translated stop position.

Coordinate convention: positions are always reported in the original genome's
1-based forward-strand numbering, regardless of which strand the ORF is on.
For a "+"-strand ORF, start < end (5' -> 3' matches increasing coordinates).
For a "-"-strand ORF, start > end: "start" is the genomic position of the
first base of the ATG as read 5'->3' on the minus strand (i.e. the *higher*
forward-strand coordinate), and "end" is the last base of the stop codon
(the *lower* forward-strand coordinate) -- this keeps "start"/"end" meaning
"where translation begins/ends" on each ORF's own strand.
"""
import sys
import os
import csv
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from common import read_fasta_records, reverse_complement, translate_codon, STOP_CODONS

MIN_PROTEIN_LEN = 100


def find_orfs_in_frame(seq, frame):
    """Yield (start0, end0_exclusive, protein) for ORFs in this 0-based frame,
    where start0/end0 are 0-based positions into `seq` (the strand-local
    sequence), end0_exclusive is one past the last base of the stop codon."""
    n = len(seq)
    codon_starts = range(frame, n - 2, 3)
    active_starts = []  # 0-based positions (into seq) of open ATGs
    for pos in codon_starts:
        codon = seq[pos:pos + 3]
        if codon == "ATG":
            active_starts.append(pos)
        elif codon in STOP_CODONS:
            if active_starts:
                stop_end = pos + 3
                for start in active_starts:
                    protein_codons = range(start, pos, 3)
                    protein = "".join(translate_codon(seq[c:c + 3]) for c in protein_codons)
                    yield start, stop_end, protein
                active_starts = []
    # any positions left in active_starts have no downstream in-frame stop -> discarded


def find_orfs_one_strand(seq):
    for frame in range(3):
        yield from find_orfs_in_frame(seq, frame)


def build_orf_table(genome_seq):
    genome_len = len(genome_seq)
    rows = []  # dict: strand, start(1-based genome coord), end, protein, length

    # Forward strand
    for start0, end0, protein in find_orfs_one_strand(genome_seq):
        rows.append({
            "strand": "+",
            "start": start0 + 1,
            "end": end0,
            "protein": protein,
        })

    # Reverse strand: work on the reverse complement, then map coords back
    rc_seq = reverse_complement(genome_seq)
    for start0, end0, protein in find_orfs_one_strand(rc_seq):
        # rc index i (0-based) <-> genome 1-based position (genome_len - i)
        genome_start = genome_len - start0        # position of ATG's first base, minus strand
        genome_end = genome_len - end0 + 1         # position of stop codon's last base
        rows.append({
            "strand": "-",
            "start": genome_start,
            "end": genome_end,
            "protein": protein,
        })

    return rows


def write_csv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["strand", "id", "start", "end", "protein"])
        for i, row in enumerate(rows, start=1):
            writer.writerow([row["strand"], i, row["start"], row["end"], row["protein"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_path")
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    records = list(read_fasta_records(args.fasta_path))
    header, genome_seq = records[0]
    genome_seq = genome_seq.upper()
    print("Genome:", header, "length:", len(genome_seq))

    all_rows = build_orf_table(genome_seq)
    print("Total ORFs found (both strands, complete only):", len(all_rows))

    long_rows = [r for r in all_rows if len(r["protein"]) > MIN_PROTEIN_LEN]
    long_rows.sort(key=lambda r: (0 if r["strand"] == "+" else 1, r["start"]))
    print(f"ORFs with protein length > {MIN_PROTEIN_LEN} aa:", len(long_rows))

    out_csv = os.path.join(args.outdir, "q6_orfs.csv")
    write_csv(long_rows, out_csv)
    print("Wrote:", out_csv)

    longest = max(all_rows, key=lambda r: len(r["protein"]))
    print("\n=== (c) Longest ORF across both strands ===")
    print("Strand:", longest["strand"])
    print("Genome start:", longest["start"], "end:", longest["end"])
    nt_length = len(longest["protein"]) * 3 + 3  # + stop codon
    print("Protein length (aa, excl. stop):", len(longest["protein"]))
    print("Full ORF length (nt, incl. stop codon):", nt_length)
    print("Protein sequence:")
    print(longest["protein"])
```

Check: translating `ATGGCCATGGCGCCCAGAACTGGGCCCTGA` with our codon table gives
`MAMAPRTGP` (stop codon excluded), matching the assignment's worked example.

### (a) Forward-strand ORFs

All forward-strand ORFs longer than 100 amino acids are written to
`results/q6_orfs.csv` (columns: `strand,id,start,end,protein`), submitted as
a separate CSV file alongside this document, per the assignment's
instructions.

### (b) Reverse-strand ORFs

Reverse-strand ORFs are found by running the same frame-scanning routine on
the reverse complement of the genome, then mapping coordinates back to the
original forward-strand numbering (see convention above) and labelling them
`-`. These rows are appended into the *same* `results/q6_orfs.csv` file
(the script builds both strands' ORFs before writing the CSV once).

Combined total across both strands: **153321** complete ORFs found (all
lengths); **30894** of these translate to proteins longer than 100 amino
acids and appear in `results/q6_orfs.csv`.

### (c) The single longest ORF

The longest ORF across both strands is on the **forward (`+`) strand**,
genome position **2044911–2052014** (1-based), with:

- Protein length: **2367 amino acids**
- Full ORF length (including the stop codon): **7104 nt**

Its translated protein sequence is:

```
MLARSGKVSMATKKRSGEEINDRQILCGMGIKLRRLTAGICLITQLAFPMAAAAQGVVNAATQQPVPAQIAIANANTVPYTLGALESAQSVAERFGISVAELRKLNQFRTFARGFDNVRQGDELDVPAQVSEKKLTPPPGNSSDNLEQQIASTSQQIGSLLAEDMNSEQAANMARGWASSQASGAMTDWLSRFGTARITLGVDEDFSLKNSQFDFLHPWYETPDNLFFSQHTLHRTDERTQINNGLGWRHFTPTWMSGINFFFDHDLSRYHSRAGIGAEYWRDYLKLSSNGYLRLTNWRSAPELDNDYEARPANGWDVRAESWLPAWPHLGGKLVYEQYYGDEVALFDKDDRQSNPHAITAGLNYTPFPLMTFSAEQRQGKQGENDTRFAVDFTWQPGSAMQKQLDPNEVAARRSLAGSRYDLVDRNNNIVLEYRKKELVRLTLTDPVTGKSGEVKSLVSSLQTKYALKGYNVEATALEAAGGKVVTTGKDILVTLPAYRFTSTPETDNTWPIEVTAEDVKGNLSNREQSMVVVQAPTLSQKDSSVSLSTQTLNADSHSTATLTFIAHDAAGNPVVGLVLSTRHEGVQDITLSDWKDNGDGSYTQILTTGAMSGTLTLMPQLNGVDAAKAPAVVNIISVSSSRTHSSIKIDKDRYLSGNPIEVTVELRDENDKPVKEQKQQLNNAVSIDNVKPGVTTDWKETADGVYKATYTAYTKGSGLTAKLLMQNWNEDLHTAGFIIDANPQSAKIATLSASNNGVLANENAANTVSVNVADEGSNPINDHTVTFAVLSGSATSFNNQNTAKTDVNGLATFDLKSSKQEDNTVEVTLENGVKQTLIVSFVGDSSTAQVDLQKSKNEVVADGNDSVTMTATVRDAKGNLLNDVMVTFNVNSAEAKLSQTEVNSHDGIATATLTSLKNGDYRVTASVSSGSQANQQVNFIGDQSTAALTLSVPSGDITVTNTAPQYMTATLQDKNGNPLKDKEITFSVPNDVASKFSISNGGKGMTDSNGVAIASLTGTLAGTHMIMARLANSNVSDAQPMTFVADKDRAVVVLQTSKAEIIGNGVDETTLTATVKDPSNHPVAGITVNFTMPQDVAANFTLENNGIAITQANGEAHVTLKGKKAGTHTVTATLGNNNTSDSQPVTFVADKASAQVVLQISKDEITGNGVDSATLTATVKDQFDNEVNNLPVTFSSASSGLTLTPGVSNTNESGIAQATLAGVAFGEKTVTASLANNGASDNKTVHFIGDTAAAKIIELAPVPDSIIAGTPQNSSGSVITATVVDNNGFPVKGVTVNFTSNAATAEMTNGGQAVTNEQGKATVTYTNTRSSIESGARPDTVEASLENGSSTLSTSINVNADASTAHLTLLQALFDTVSAGETTSLYIEVKDNYGNGVPQQEVTLSVSPSEGVTPSNNAIYTTNHDGNFYASFTATKAGVYQLTATLENGDSMQQTVTYVPNVANAEITLAASKDPVIADNNDLTTLTATVADTEGNAIANTEVTFTLPEDVKANFTLSDGGKVITDAEGKAKVTLKGTKAGAHTVTASMTGGKSEQLVVNFIADTLTAQVNLNVTEDNFIANNVGMTRLQATVTDGNGNPLANEAVTFTLPADVSASFTLGQGGSAITDINGKAEVTLSGTKSGTYPVTVSVNNYGVSDTKQVTLIADAGTAKLASLTSVYSFVVSTTEGATMTASVTDANGNPVEGIKVNFRGTSVTLSSTSVETDDRGFAEILVTSTEVGLKTVSASLADKPTEVISRLLNASADVNSATITSLEIPEGQVMVAQDVAVKAHVNDQFGNPVAHQPVTFSAEPSSQMIISQNTVSTNTQGVAEVTMTPERNGSYMVKASLPNGASLEKQLEAIDEKLTLTASSPLIGVYAPTGATLTATLTSANGTPVEGQVINFSVTPEGATLSGGKVRTNSSGQAPVVLTSNKVGTYTVTASFHNGVTIQTQTTVKVTGNSSTAHVASFIADPSTIAATNTDLSTLKATVEDGSGNLIEGLTVYFALKSGSATLTSLTAVTDQNGIATTSVKGAMTGSVTVSAVTTAGGMQTVDITLVAGPADTSQSVLKSNRSSLKGDYTDSAELRLVLHDISGNPIKVSEGMEFVQSGTNVPYIKISAIDYSLNINGDYKATVTGGGEGIATLIPVLNGVHQAGLSTTIQFTRAEDKIMSGTVSVNGTDLPTTTFPSQGFTGAYYQLNNDNFAPGKTAADYEFSSSASWVDVDATGKVTFKNVGSNSERITATPKSGGPSYVYEIRVKSWWVNAGEAFMIYSLAENFCSSNGYTLPRANYLNHCSSRGIGSLYSEWGDMGHYTTDAGFQSNMYWSSSPANSSEQYVVSLATGDQSVFEKLGFAYATCYKNL
```

**Protein identification:** This ORF corresponds to ***yeeJ*** (an inverse
autotransporter adhesin). Several independent lines of evidence support this:

1. *Length and genome position*: `yeeJ` is a well-known outlier in the
   E. coli K-12 MG1655 genome for being unusually long; published work on
   MG1655 reports the YeeJ protein at 2358 amino acids (our translation:
   2367 aa — the small difference is consistent with normal
   strain/annotation-version variation in exact start-codon usage, not a
   different gene).
2. *Domain architecture*: the translated sequence is dominated by dozens of
   short, highly repetitive ~90–100 residue blocks (e.g. repeated
   `...TVTASL...NG...` / `...VTFTLP...` motifs), matching the literature
   description of YeeJ's passenger domain as containing **13 bacterial
   immunoglobulin-like ("Big") domain repeats**, followed by a distinct
   C-terminal region — consistent with the last part of our sequence
   (`...LGFAYATCYKNL`) corresponding to YeeJ's C-type lectin domain.
3. *Function*: YeeJ is described as an inverse autotransporter that binds
   peptidoglycan and promotes biofilm formation in *E. coli* (Meuskens et
   al., *Scientific Reports* 2017, https://www.nature.com/articles/s41598-017-10902-0).

---

## Notes on data and reproducibility

All five input files (`data/Q1/input1.fa`, `data/Q1/chr22.fa`,
`data/Q3/BRCA_aligned.fa`, `data/Q4/DHFR_aligned.fa`,
`data/Q6/ecoli_genome.fa`) and every script (`src/*.py`) referenced above are
included in this submission and can be re-run directly, e.g.:

```bash
python src/q1_count_nucleotides.py data/Q1/chr22.fa
python src/q2_find_motif.py data/Q1/chr22.fa CGTAACC
python src/q3_consensus_profile.py data/Q3/BRCA_aligned.fa
python src/q4_covariation.py data/Q4/DHFR_aligned.fa
python src/q5_ecorv.py data/Q1/chr22.fa --outdir results
python src/q6_orf_finder.py data/Q6/ecoli_genome.fa --outdir results
```
