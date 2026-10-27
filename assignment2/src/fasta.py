"""Minimal FASTA reader shared by the Assignment 2 scripts."""


def read_fasta(path):
    """Return a list of (header, sequence) tuples; sequences are upper-cased."""
    records = []
    header, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records
