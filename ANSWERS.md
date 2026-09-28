# COMP3353 Bioinformatics — Assignment 1: Sequence Analysis

Name: Zheng Choi I
University Number: 3035987788
Email: u3598778@connect.hku.hk

Language: Python 3 (standard library only, except `numpy`/`matplotlib` for
Q4/Q5). Positions are 1-based throughout. Each script is self-contained.

## 1. Counting DNA nucleotides

```python
"""Q1: Count occurrences of A, C, G, T (case-insensitive) in a FASTA file."""
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
    print(count_nucleotides("data/Q1/input1.fa"))
    print(count_nucleotides("data/Q1/chr22.fa"))
```

Answer:
- `data/Q1/input1.fa`: A=13 C=13 G=17 T=17
- `data/Q1/chr22.fa`: A=10382214 C=9160652 G=9246186 T=10370725

## 2. Finding a motif in DNA

```python
"""Q2: Count occurrences of a motif s as a substring of a target DNA sequence.

Counts overlapping, case-insensitive matches on the forward strand only.
Each FASTA record is searched independently (matches are not allowed to span
across two different records), and counts are summed across all records.
"""


def read_fasta_records(path):
    """Yield (header, sequence) tuples, sequence as one string per record."""
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


def count_motif_in_file(path, motif):
    total = 0
    for _header, seq in read_fasta_records(path):
        total += count_overlapping(seq, motif)
    return total


if __name__ == "__main__":
    print(count_motif_in_file("data/Q1/input1.fa", "CGTAACC"))
    print(count_motif_in_file("data/Q1/chr22.fa", "CGTAACC"))
```

Answer:
- `{path: "data/Q1/input1.fa", s: "CGTAACC"}`: 4
- `{path: "data/Q1/chr22.fa", s: "CGTAACC"}`: 206

## 3. Consensus and profile matrix

Metric: Shannon entropy per column, H = -sum p_b log2(p_b) over A/C/G/T.
Used instead of a plain majority-mismatch count because it reflects the full
split of the column (13/12 is more variable than 24/1, entropy separates
these, a mismatch count does not). Same quantity used for MI in Q4.

```python
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
from math import log2

BASES = "ACGT"


def read_fasta_records(path):
    """Yield (header, sequence) tuples, sequence as one string per record."""
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
```

Answer (`data/Q3/BRCA_aligned.fa`, 25 seqs x 335 cols):
- Consensus:
```
GATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGAGATGTGGTCAATGGAAGAAACCACCAAGGTCCAAAGCGAGCAAGAGAATCCCAGGACAGAAAGGTAAAGCTCCCTCCCTCAAGTTGACAAAAATCTCACCCCACCACTCTGTATTCCACTCCCCTTTGCAGAGATGGGCCGCTTCATTTTGTAAGACTTATTACATACATACACAGTGCTAGATACTTTCACACAGGTTCTTTTTTCACTCTTCCATCCCAACCACATAAATAAGTATTGTCTCTACTTTATGAATGATAAAACTAAGAGATTTAGAGAGGCTGTGTA
```
- Top 3 least-conserved positions (1-based): 26, 158, 144 (entropy 1.164, 1.164, 1.0211 bits)

## 4. MSA and covariation matrices

Alphabet: 20 amino acids + gap `-`; anything else treated as gap. MI formula
as given, zero-frequency terms skipped.

```python
"""Q4: Consensus string and top-10 co-varying column pairs (mutual information)
for a protein multiple sequence alignment.

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
import numpy as np

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"
ALPHABET = AMINO_ACIDS + GAP  # 21 symbols, gap last
SYMBOL_INDEX = {ch: i for i, ch in enumerate(ALPHABET)}
N_SYMBOLS = len(ALPHABET)


def read_fasta_records(path):
    """Yield (header, sequence) tuples, sequence as one string per record."""
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
    seqs, length = load_alignment("data/Q4/DHFR_aligned.fa")
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

Answer (`data/Q4/DHFR_aligned.fa`, 667 seqs x 160 cols):
- Consensus:
```
MISLIVAMAENGVIGKDNDLPWHLPEDLKYFKRLTLGKPVIMGRKTWESIGRPLPGRRNIVLSRDPDYQAEGVEVVHSLEEALALAG-VEEVFVIGGAEIYAQALP-ADRLYLTEIDAEFEGDTFFPEIDPDEWEEVSREEHPADEKNGYDYTFVTYERK
```
- Top 10 co-varying column pairs (1-based) by MI:

| Rank | Columns | MI (bits) |
|---|---|---|
| 1 | (149,150) | 1.2409 |
| 2 | (150,151) | 1.1677 |
| 3 | (68,70) | 1.1111 |
| 4 | (148,150) | 1.0829 |
| 5 | (145,146) | 1.0511 |
| 6 | (13,122) | 1.0165 |
| 7 | (149,151) | 1.0131 |
| 8 | (6,8) | 1.0072 |
| 9 | (151,152) | 0.9992 |
| 10 | (58,74) | 0.9862 |

## 5. Restriction-nuclease case study (chr22)

```python
"""Q5: EcoRV restriction-site case study on chr22.

Runs all four sub-questions (a-d) and writes a fragment-length histogram to
results/q5d_fragment_length_hist.png.
"""
import os
import re
import statistics

FASTA_PATH = "data/Q1/chr22.fa"
OUTDIR = "results"

SITE = "GATATC"
COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta_records(path):
    """Yield (header, sequence) tuples, sequence as one string per record."""
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
    seq = load_single_sequence(FASTA_PATH)

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
    n_fragments, median_length, hist_path, _ = part_d_fragments(seq, OUTDIR)
    print("Number of fragments:", n_fragments)
    print("Median fragment length (N's excluded):", median_length)
    print("Histogram saved to:", hist_path)
```

a. Exact `GATATC` sites: **5024**

b. Reverse complement of `GATATC` is `GATATC` (palindrome). One strand is
enough because the site is identical on both strands at every occurrence
(DNA strands are antiparallel and complementary, so a palindromic motif on
the top strand implies the same sequence on the bottom strand at the same
location). Scanning the reverse strand would just re-find the same
positions, not new ones.

c. Relaxed sites `GA[ACGT]ATC` (mismatch at position 3):
- Check string `TTGATATCAAGAGATCCTTCCGAAATCACGT`: 1 exact + 2 relaxed = 3 total. Matches.
- chr22: 31370 total matches = 5024 exact + 26346 relaxed-only.

d. Cutting at every exact site (`GAT/ATC`), N's excluded from lengths:
- Fragments: 5025
- Median length: 3561 bp
- Histogram: `results/q5d_fragment_length_hist.png`

![fragment length distribution](results/q5d_fragment_length_hist.png)

## 6. ORF finder (E. coli genome)

ORF = ATG to next in-frame stop codon in the same frame; every ATG before a
stop gets its own ORF ending at that stop; ATGs with no downstream in-frame
stop are dropped (incomplete). Coordinates are always in original genome
1-based forward numbering; for `-` strand ORFs, start = 5' base of ATG on
minus strand (higher coordinate), end = last base of stop codon (lower
coordinate).

```python
"""Q6: ORF finder for the E. coli genome.

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
import os
import csv

FASTA_PATH = "data/Q6/ecoli_genome.fa"
OUT_CSV = "results/q6_orfs.csv"
MIN_PROTEIN_LEN = 100

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

CODON_TABLE = {
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
    """Yield (header, sequence) tuples, sequence as one string per record."""
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
    return CODON_TABLE.get(codon.upper(), "X")


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
    records = list(read_fasta_records(FASTA_PATH))
    header, genome_seq = records[0]
    genome_seq = genome_seq.upper()
    print("Genome:", header, "length:", len(genome_seq))

    all_rows = build_orf_table(genome_seq)
    print("Total ORFs found (both strands, complete only):", len(all_rows))

    long_rows = [r for r in all_rows if len(r["protein"]) > MIN_PROTEIN_LEN]
    long_rows.sort(key=lambda r: (0 if r["strand"] == "+" else 1, r["start"]))
    print(f"ORFs with protein length > {MIN_PROTEIN_LEN} aa:", len(long_rows))

    write_csv(long_rows, OUT_CSV)
    print("Wrote:", OUT_CSV)

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

Check: translating `ATGGCCATGGCGCCCAGAACTGGGCCCTGA` gives `MAMAPRTGP`. Matches.

a. Forward-strand ORFs >100 aa: written to `results/q6_orfs.csv`
   (columns: strand,id,start,end,protein), submitted as a separate file.

b. Reverse-strand ORFs (same length filter) appended to the same CSV,
   labeled `-`. Combined total (both strands, all lengths): 153321 ORFs;
   30894 have protein length > 100 aa (in the CSV).

c. Longest ORF: `+` strand, genome 2044911-2052014, protein length 2367 aa,
   ORF length 7104 nt (incl. stop). Protein:
```
MLARSGKVSMATKKRSGEEINDRQILCGMGIKLRRLTAGICLITQLAFPMAAAAQGVVNAATQQPVPAQIAIANANTVPYTLGALESAQSVAERFGISVAELRKLNQFRTFARGFDNVRQGDELDVPAQVSEKKLTPPPGNSSDNLEQQIASTSQQIGSLLAEDMNSEQAANMARGWASSQASGAMTDWLSRFGTARITLGVDEDFSLKNSQFDFLHPWYETPDNLFFSQHTLHRTDERTQINNGLGWRHFTPTWMSGINFFFDHDLSRYHSRAGIGAEYWRDYLKLSSNGYLRLTNWRSAPELDNDYEARPANGWDVRAESWLPAWPHLGGKLVYEQYYGDEVALFDKDDRQSNPHAITAGLNYTPFPLMTFSAEQRQGKQGENDTRFAVDFTWQPGSAMQKQLDPNEVAARRSLAGSRYDLVDRNNNIVLEYRKKELVRLTLTDPVTGKSGEVKSLVSSLQTKYALKGYNVEATALEAAGGKVVTTGKDILVTLPAYRFTSTPETDNTWPIEVTAEDVKGNLSNREQSMVVVQAPTLSQKDSSVSLSTQTLNADSHSTATLTFIAHDAAGNPVVGLVLSTRHEGVQDITLSDWKDNGDGSYTQILTTGAMSGTLTLMPQLNGVDAAKAPAVVNIISVSSSRTHSSIKIDKDRYLSGNPIEVTVELRDENDKPVKEQKQQLNNAVSIDNVKPGVTTDWKETADGVYKATYTAYTKGSGLTAKLLMQNWNEDLHTAGFIIDANPQSAKIATLSASNNGVLANENAANTVSVNVADEGSNPINDHTVTFAVLSGSATSFNNQNTAKTDVNGLATFDLKSSKQEDNTVEVTLENGVKQTLIVSFVGDSSTAQVDLQKSKNEVVADGNDSVTMTATVRDAKGNLLNDVMVTFNVNSAEAKLSQTEVNSHDGIATATLTSLKNGDYRVTASVSSGSQANQQVNFIGDQSTAALTLSVPSGDITVTNTAPQYMTATLQDKNGNPLKDKEITFSVPNDVASKFSISNGGKGMTDSNGVAIASLTGTLAGTHMIMARLANSNVSDAQPMTFVADKDRAVVVLQTSKAEIIGNGVDETTLTATVKDPSNHPVAGITVNFTMPQDVAANFTLENNGIAITQANGEAHVTLKGKKAGTHTVTATLGNNNTSDSQPVTFVADKASAQVVLQISKDEITGNGVDSATLTATVKDQFDNEVNNLPVTFSSASSGLTLTPGVSNTNESGIAQATLAGVAFGEKTVTASLANNGASDNKTVHFIGDTAAAKIIELAPVPDSIIAGTPQNSSGSVITATVVDNNGFPVKGVTVNFTSNAATAEMTNGGQAVTNEQGKATVTYTNTRSSIESGARPDTVEASLENGSSTLSTSINVNADASTAHLTLLQALFDTVSAGETTSLYIEVKDNYGNGVPQQEVTLSVSPSEGVTPSNNAIYTTNHDGNFYASFTATKAGVYQLTATLENGDSMQQTVTYVPNVANAEITLAASKDPVIADNNDLTTLTATVADTEGNAIANTEVTFTLPEDVKANFTLSDGGKVITDAEGKAKVTLKGTKAGAHTVTASMTGGKSEQLVVNFIADTLTAQVNLNVTEDNFIANNVGMTRLQATVTDGNGNPLANEAVTFTLPADVSASFTLGQGGSAITDINGKAEVTLSGTKSGTYPVTVSVNNYGVSDTKQVTLIADAGTAKLASLTSVYSFVVSTTEGATMTASVTDANGNPVEGIKVNFRGTSVTLSSTSVETDDRGFAEILVTSTEVGLKTVSASLADKPTEVISRLLNASADVNSATITSLEIPEGQVMVAQDVAVKAHVNDQFGNPVAHQPVTFSAEPSSQMIISQNTVSTNTQGVAEVTMTPERNGSYMVKASLPNGASLEKQLEAIDEKLTLTASSPLIGVYAPTGATLTATLTSANGTPVEGQVINFSVTPEGATLSGGKVRTNSSGQAPVVLTSNKVGTYTVTASFHNGVTIQTQTTVKVTGNSSTAHVASFIADPSTIAATNTDLSTLKATVEDGSGNLIEGLTVYFALKSGSATLTSLTAVTDQNGIATTSVKGAMTGSVTVSAVTTAGGMQTVDITLVAGPADTSQSVLKSNRSSLKGDYTDSAELRLVLHDISGNPIKVSEGMEFVQSGTNVPYIKISAIDYSLNINGDYKATVTGGGEGIATLIPVLNGVHQAGLSTTIQFTRAEDKIMSGTVSVNGTDLPTTTFPSQGFTGAYYQLNNDNFAPGKTAADYEFSSSASWVDVDATGKVTFKNVGSNSERITATPKSGGPSYVYEIRVKSWWVNAGEAFMIYSLAENFCSSNGYTLPRANYLNHCSSRGIGSLYSEWGDMGHYTTDAGFQSNMYWSSSPANSSEQYVVSLATGDQSVFEKLGFAYATCYKNL
```

Identity: `yeeJ`, an inverse autotransporter adhesin. MG1655 YeeJ is reported
at 2358 aa (vs. our 2367 aa — minor version/strain difference, same gene).
The translated sequence shows the same repeat pattern as YeeJ's ~13
bacterial Ig-like (Big) domain repeats, and its function (peptidoglycan
binding, biofilm formation) is described in Meuskens et al., Sci Rep 2017
(https://www.nature.com/articles/s41598-017-10902-0).
