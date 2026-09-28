import re
import statistics
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SITE = "GATATC"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

lines = []
for line in open("data/Q1/chr22.fa"):
    if not line.startswith(">"):
        lines.append(line.strip().upper())
seq = "".join(lines)

# (a) exact sites
count_exact = 0
i = seq.find(SITE)
while i != -1:
    count_exact += 1
    i = seq.find(SITE, i + 1)
print("(a) exact GATATC sites:", count_exact)

# (b) reverse complement
rc = SITE.translate(COMPLEMENT)[::-1]
print("(b) reverse complement of GATATC:", rc, "-> palindrome:", rc == SITE)

# (c) relaxed sites, mismatch allowed at position 3: GA[N]ATC
count_relaxed = len(re.findall(r"(?=(GA[ACGT]ATC))", seq))
print("(c) total GA[N]ATC sites:", count_relaxed,
      "(exact:", count_exact, ", relaxed-only:", count_relaxed - count_exact, ")")

# (d) cut at every exact site (after the 3rd base: GAT/ATC), N's excluded from lengths
cuts = []
i = seq.find(SITE)
while i != -1:
    cuts.append(i + 3)
    i = seq.find(SITE, i + len(SITE))

boundaries = [0] + cuts + [len(seq)]
fragment_lengths = []
for start, end in zip(boundaries, boundaries[1:]):
    fragment = seq[start:end]
    fragment_lengths.append(len(fragment) - fragment.count("N"))

print("(d) number of fragments:", len(fragment_lengths))
print("(d) median fragment length:", statistics.median(fragment_lengths))

positive = [l for l in fragment_lengths if l > 0]
bins = np.logspace(0, np.log10(max(positive)), 50)
plt.hist(positive, bins=bins)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Fragment length (bp, N excluded)")
plt.ylabel("Number of fragments")
plt.title("EcoRV digest of chr22: fragment length distribution")
plt.savefig("results/q5d_fragment_length_hist.png", dpi=150)
