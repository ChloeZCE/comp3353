def read_fasta(path):
    lines = []
    for line in open(path):
        if not line.startswith(">"):
            lines.append(line.strip())
    return "".join(lines)


def count_motif(seq, motif):
    seq = seq.upper()
    motif = motif.upper()
    count = 0
    i = seq.find(motif)
    while i != -1:
        count += 1
        i = seq.find(motif, i + 1)  # overlapping matches
    return count


print(count_motif(read_fasta("data/Q1/input1.fa"), "CGTAACC"))
print(count_motif(read_fasta("data/Q1/chr22.fa"), "CGTAACC"))
