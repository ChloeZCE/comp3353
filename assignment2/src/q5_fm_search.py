"""Q5: count occurrences of a motif with FM-index backward search.

Usage:
    python q5_fm_search.py <fasta_path> <motif> [--show-occ N]

Prints the BWT, the count array C and the occurrence matrix Occ (the full
matrix for short inputs, otherwise its first N rows), then the backward-search
steps and the number of occurrences.
"""
import sys

from bwt import bwt_encode
from fasta import read_fasta


def build_fm_index(bwt):
    alphabet = sorted(set(bwt))
    # C[c] = number of symbols in the text lexicographically smaller than c
    C, total = {}, 0
    for c in alphabet:
        C[c] = total
        total += bwt.count(c)
    # Occ[c][k] = occurrences of c in bwt[:k]   (k = 0 .. n)
    occ = {c: [0] * (len(bwt) + 1) for c in alphabet}
    for k, ch in enumerate(bwt):
        for c in alphabet:
            occ[c][k + 1] = occ[c][k] + (ch == c)
    return C, occ


def backward_search(motif, C, occ, n, trace=None):
    """Return the number of occurrences; SA interval is [lo, hi)."""
    lo, hi = 0, n
    for c in reversed(motif):
        if c not in C:
            return 0
        lo, hi = C[c] + occ[c][lo], C[c] + occ[c][hi]
        if trace is not None:
            trace.append((c, lo, hi))
        if lo >= hi:
            return 0
    return hi - lo


def main():
    args = sys.argv[1:]
    show = 30
    if "--show-occ" in args:
        i = args.index("--show-occ")
        show = int(args[i + 1])
        del args[i:i + 2]
    if len(args) != 2:
        sys.exit("Usage: python q5_fm_search.py <fasta_path> <motif> [--show-occ N]")
    path, motif = args[0], args[1].upper()

    bwt = bwt_encode(read_fasta(path)[0][1])
    n = len(bwt)
    C, occ = build_fm_index(bwt)
    alphabet = sorted(C)

    print(f"BWT (length {n}):")
    print(bwt if n <= 200 else bwt[:200] + "...")
    print("\nCount array C:")
    print("  " + "  ".join(f"C[{c}]={C[c]}" for c in alphabet))
    rows = n + 1 if n + 1 <= show else show
    print(f"\nOccurrence matrix Occ(c, k), k = 0..{n}"
          + ("" if rows == n + 1 else f" (first {rows} rows shown)") + ":")
    print("   k  L[k-1]  " + "  ".join(f"{c:>6}" for c in alphabet))
    for k in range(rows):
        prev = bwt[k - 1] if k else "-"
        print(f"{k:4d}  {prev:>6}  " + "  ".join(f"{occ[c][k]:6d}" for c in alphabet))
    print(f"  (final row k={n}: " + ", ".join(f"{c}={occ[c][n]}" for c in alphabet) + ")")

    trace = []
    count = backward_search(motif, C, occ, n, trace)
    print(f"\nBackward search for {motif}: start with [lo, hi) = [0, {n})")
    for c, lo, hi in trace:
        print(f"  char {c}: lo = C[{c}] + Occ({c},lo), hi = C[{c}] + Occ({c},hi)"
              f" -> [{lo}, {hi})  size {max(hi - lo, 0)}")
    print(f"\nNumber of occurrences of {motif}: {count}")


if __name__ == "__main__":
    main()
