from collections import Counter
from math import log2

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


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


def clean(ch):
    return ch if ch in AMINO_ACIDS else "-"


seqs = read_fasta("data/Q4/DHFR_aligned.fa")
n = len(seqs)
length = len(seqs[0])

columns = [[clean(seq[i]) for seq in seqs] for i in range(length)]
freqs = [Counter(col) for col in columns]
for f in freqs:
    for k in f:
        f[k] /= n

consensus = "".join(c.most_common(1)[0][0] for c in freqs)

results = []
for i in range(length):
    for j in range(i + 1, length):
        joint = Counter(zip(columns[i], columns[j]))
        mi = 0.0
        for (a, b), count in joint.items():
            pij = count / n
            mi += pij * log2(pij / (freqs[i][a] * freqs[j][b]))
        results.append((mi, i + 1, j + 1))

results.sort(reverse=True)

print("Consensus string:")
print(consensus)
print("Top 10 co-varying column pairs by mutual information:")
for mi, i, j in results[:10]:
    print(f"  columns ({i}, {j}): MI = {mi:.4f} bits")
