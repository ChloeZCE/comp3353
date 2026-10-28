import sys


def read_fasta(path):
    seq = ""
    for line in open(path):
        if not line.startswith(">"):
            seq += line.strip().upper()
    return seq


def suffix_array(text):
    # Sort suffix start positions by "prefix doubling": first by their first
    # letter, then by their first 2, 4, 8, ... letters, reusing the ranks from
    # the previous round. This never copies the suffixes themselves, so it
    # also works for a whole genome.
    n = len(text)
    rank = [ord(c) for c in text]
    order = list(range(n))
    k = 1
    while True:
        key = lambda i: (rank[i], rank[i + k] if i + k < n else -1)
        order.sort(key=key)
        new_rank = [0] * n
        for a, b in zip(order, order[1:]):
            new_rank[b] = new_rank[a] + (1 if key(a) != key(b) else 0)
        rank = new_rank
        if rank[order[-1]] == n - 1:  # every suffix has its own rank: done
            return order
        k *= 2


def bwt_encode(seq):
    text = seq + "$"
    # Sorting all rotations of text is the same as sorting its suffixes,
    # because "$" is unique and smaller than A/C/G/T.
    order = suffix_array(text)
    # the last letter of the rotation starting at i is text[i - 1]
    return "".join(text[i - 1] for i in order)


def bwt_decode(bwt):
    # first column F = sorted(bwt); first[c] = row where c first appears in F
    first = {}
    for row, c in enumerate(sorted(bwt)):
        if c not in first:
            first[c] = row
    # LF mapping: the k-th c in the last column is the k-th c in the first column
    seen = {}
    lf = []
    for c in bwt:
        lf.append(first[c] + seen.get(c, 0))
        seen[c] = seen.get(c, 0) + 1
    # row 0 starts with "$", so bwt[0] is the last base of the sequence;
    # following LF walks the sequence backwards
    seq = ""
    row = 0
    for _ in range(len(bwt) - 1):
        seq = bwt[row] + seq
        row = lf[row]
    return seq


if __name__ == "__main__":
    mode, path = sys.argv[1], sys.argv[2]
    if mode == "encode":
        print(bwt_encode(read_fasta(path)))
    else:
        print(bwt_decode(open(path).read().strip()))
