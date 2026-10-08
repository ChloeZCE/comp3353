"""Report positions in a `samtools mpileup -f ref` stream whose major allele
differs from the reference (used for Lab 2 Q9)."""

import re
import sys
from collections import Counter

INDEL = re.compile(r"[+-](\d+)")


def base_counts(ref, bases):
    counts = Counter()
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == "^":            # read start; next char is MAPQ
            i += 2
            continue
        if c in "+-":           # indel: skip its length + sequence
            n = INDEL.match(bases, i).group(1)
            i += 1 + len(n) + int(n)
            continue
        if c in ".,":
            counts[ref.upper()] += 1
        elif c.upper() in "ACGTN":
            counts[c.upper()] += 1
        elif c == "*":
            counts["*"] += 1
        i += 1
    return counts


def main(min_frac=0.5, min_depth=10):
    found = False
    for line in sys.stdin:
        chrom, pos, ref, depth, bases = line.split("\t")[:5]
        counts = base_counts(ref, bases)
        total = sum(counts.values())
        if total < min_depth:
            continue
        allele, n = counts.most_common(1)[0]
        if allele != ref.upper() and n / total >= min_frac:
            found = True
            detail = ", ".join(f"{b}:{k}" for b, k in counts.most_common())
            print(f"{chrom}:{pos}  ref={ref.upper()}  alt={allele}  "
                  f"{n}/{total} = {n / total:.1%}  [{detail}]")
    if not found:
        print("no position with a non-reference major allele")


if __name__ == "__main__":
    main()
