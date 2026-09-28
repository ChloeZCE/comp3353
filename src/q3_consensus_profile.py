from math import log2

BASES = "ACGT"


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


seqs = read_fasta("data/Q3/BRCA_aligned.fa")
length = len(seqs[0])

# profile matrix: count of each base at each column
profile = {b: [0] * length for b in BASES}
for seq in seqs:
    for i, ch in enumerate(seq):
        if ch in BASES:
            profile[ch][i] += 1

# consensus + entropy per column (entropy = variability metric)
consensus = ""
entropy = []
for i in range(length):
    counts = [profile[b][i] for b in BASES]
    total = sum(counts)
    consensus += BASES[counts.index(max(counts))]
    h = -sum((c / total) * log2(c / total) for c in counts if c > 0)
    entropy.append(h)

top3 = sorted(range(length), key=lambda i: -entropy[i])[:3]
top3_positions = [i + 1 for i in top3]  # 1-based

print("Consensus string:")
print(consensus)
print("Top 3 least-conserved positions (1-based):", top3_positions)
print("Entropy at those positions:", [round(entropy[i], 4) for i in top3])
