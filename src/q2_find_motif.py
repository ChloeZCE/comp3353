"""Q2: Count occurrences of a motif s as a substring of a target DNA sequence.

Counts overlapping, case-insensitive matches on the forward strand only.
Each FASTA record is searched independently (matches are not allowed to span
across two different records), and counts are summed across all records.
"""


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


def count_overlapping(seq, motif):
    """Count overlapping, case-insensitive occurrences of motif in seq."""
    seq_u = seq.upper()
    motif_u = motif.upper()
    count = 0
    start = 0
    while True:
        idx = seq_u.find(motif_u, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1
    return count


def count_motif_in_file(path, motif):
    total = 0
    for _header, seq in read_fasta_records(path):
        total += count_overlapping(seq, motif)
    return total


if __name__ == "__main__":
    print(count_motif_in_file("data/Q1/input1.fa", "CGTAACC"))
    print(count_motif_in_file("data/Q1/chr22.fa", "CGTAACC"))
