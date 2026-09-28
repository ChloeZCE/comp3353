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
