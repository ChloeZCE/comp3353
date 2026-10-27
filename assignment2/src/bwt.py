"""Q4: Burrows-Wheeler Transform encoder and decoder ('$' is the end marker).

Usage:
    python bwt.py encode <fasta_path>      # prints the BWT of the sequence
    python bwt.py decode <bwt_txt_path>    # prints the original sequence
"""
import sys

from fasta import read_fasta

SENTINEL = "$"


def suffix_array(text):
    """Suffix array by prefix doubling: O(n log^2 n), fine for whole genomes."""
    n = len(text)
    rank = [ord(c) for c in text]
    sa = list(range(n))
    k = 1
    while True:
        key = lambda i: (rank[i], rank[i + k] if i + k < n else -1)
        sa.sort(key=key)
        new_rank = [0] * n
        for a, b in zip(sa, sa[1:]):
            new_rank[b] = new_rank[a] + (key(a) != key(b))
        rank = new_rank
        if rank[sa[-1]] == n - 1:
            return sa
        k *= 2


def bwt_encode(seq):
    """BWT = last column of the sorted rotations of seq + '$'.

    Sorting the rotations of a '$'-terminated string is the same as sorting
    its suffixes, and the last-column symbol of rotation i is text[i-1].
    """
    text = seq + SENTINEL
    return "".join(text[i - 1] for i in suffix_array(text))


def bwt_decode(bwt):
    """Invert the BWT with the LF mapping (last-to-first column)."""
    # first[c] = row of the first occurrence of c in the sorted first column F
    counts = {}
    for c in bwt:
        counts[c] = counts.get(c, 0) + 1
    first, total = {}, 0
    for c in sorted(counts):
        first[c] = total
        total += counts[c]
    # LF[i] = row in F holding the same character occurrence as L[i]
    seen, lf = {}, [0] * len(bwt)
    for i, c in enumerate(bwt):
        lf[i] = first[c] + seen.get(c, 0)
        seen[c] = seen.get(c, 0) + 1
    # Row 0 is the rotation starting with '$'; its last char is the final
    # symbol of the sequence.  Walk LF backwards through the original text.
    out, row = [], 0
    for _ in range(len(bwt) - 1):
        out.append(bwt[row])
        row = lf[row]
    return "".join(reversed(out))


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("encode", "decode"):
        sys.exit("Usage: python bwt.py encode <fasta_path> | decode <bwt_txt_path>")
    mode, path = sys.argv[1], sys.argv[2]
    if mode == "encode":
        print(bwt_encode(read_fasta(path)[0][1]))
    else:
        with open(path) as f:
            print(bwt_decode("".join(f.read().split())))
