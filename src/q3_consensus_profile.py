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
n = len(seqs)
length = len(seqs[0])

# profile matrix: count of each base at each column
profile = {b: [0] * length for b in BASES}
for seq in seqs:
    for i, ch in enumerate(seq):
        if ch in BASES:
            profile[ch][i] += 1

# consensus + variability per column
# variability = mismatches to consensus (N - majority count), tiebroken by
# sum-of-pairs mismatches (how many of the 25*24/2 sequence pairs disagree
# at that column) -- this separates a clean two-way split (e.g. 17/8) from
# a messier three-way split (e.g. 17/6/2) that has the same mismatch count
consensus = ""
mismatches = []
pair_mismatches = []
for i in range(length):
    counts = [profile[b][i] for b in BASES]
    total = sum(counts)
    consensus += BASES[counts.index(max(counts))]
    mismatches.append(total - max(counts))

    nonzero = [c for c in counts if c > 0]
    sp = sum(nonzero[a] * nonzero[b] for a in range(len(nonzero)) for b in range(a + 1, len(nonzero)))
    pair_mismatches.append(sp)

top3 = sorted(range(length), key=lambda i: (-mismatches[i], -pair_mismatches[i]))[:3]
top3_positions = [i + 1 for i in top3]  # 1-based

print("Consensus string:")
print(consensus)
print("Top 3 least-conserved positions (1-based):", top3_positions)
print("Mismatches to consensus at those positions:", [mismatches[i] for i in top3])
print("Sum-of-pairs mismatches at those positions:", [pair_mismatches[i] for i in top3])
