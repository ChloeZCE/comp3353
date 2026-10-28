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


def hamming(s, t):
    if len(s) != len(t):
        raise ValueError("the two sequences must have the same length")
    count = 0
    for a, b in zip(s, t):
        if a != b:
            count += 1
    return count


s, t = read_fasta(sys.argv[1])[:2]
print(hamming(s, t))
