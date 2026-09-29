from collections import Counter


def count_nucleotides(path):
    counts = Counter()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            for ch in line.strip().upper():
                if ch in "ACGT":
                    counts[ch] += 1
    return counts["A"], counts["C"], counts["G"], counts["T"]


if __name__ == "__main__":
    print(count_nucleotides("data/Q1/input1.fa"))
    print(count_nucleotides("data/Q1/chr22.fa"))
