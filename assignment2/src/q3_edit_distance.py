"""Q3: Edit (Levenshtein) distance between the two sequences in a FASTA file.

Usage:
    python q3_edit_distance.py <fasta_path>

D[i][j] = edit distance between s[:i] and t[:j]
D[0][0] = 0, D[i][0] = i, D[0][j] = j
D[i][j] = min(D[i-1][j] + 1,                     # delete s[i-1]
              D[i][j-1] + 1,                     # insert t[j-1]
              D[i-1][j-1] + (s[i-1] != t[j-1]))  # match / substitute
"""
import sys

from fasta import read_fasta


def edit_distance(s, t):
    n, m = len(s), len(t)
    # Only the previous row is needed to compute the current one, so the
    # (n+1) x (m+1) matrix is filled row by row with O(m) memory.
    prev = list(range(m + 1))                 # row i = 0: D[0][j] = j
    for i in range(1, n + 1):
        cur = [i] + [0] * m                   # D[i][0] = i
        si = s[i - 1]
        for j in range(1, m + 1):
            cur[j] = min(prev[j] + 1,
                         cur[j - 1] + 1,
                         prev[j - 1] + (si != t[j - 1]))
        prev = cur
    return prev[m]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python q3_edit_distance.py <fasta_path>")
    (_, s), (_, t) = read_fasta(sys.argv[1])[:2]
    print(edit_distance(s, t))
