import sys


def read_fasta(path):
    seqs = []
    seq = ""
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if seq:
                seqs.append(seq)
            seq = ""
        else:
            seq += line.upper()
    if seq:
        seqs.append(seq)
    return seqs


def edit_distance(s, t):
    n, m = len(s), len(t)
    # D[i][j] = edit distance between the first i bases of s and the first j bases of t
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        D[i][0] = i
    for j in range(m + 1):
        D[0][j] = j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            mismatch = 0 if s[i - 1] == t[j - 1] else 1
            D[i][j] = min(D[i - 1][j] + 1,             # deletion
                          D[i][j - 1] + 1,             # insertion
                          D[i - 1][j - 1] + mismatch)  # match / substitution
    return D[n][m]


if __name__ == "__main__":
    s, t = read_fasta(sys.argv[1])[:2]
    print(edit_distance(s, t))
