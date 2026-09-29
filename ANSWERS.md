# COMP3353 Bioinformatics — Assignment 1: Sequence Analysis

Name: Zheng Choi I

University Number: 3035987788

Email: u3598778@connect.hku.hk

## Question 1: Counting DNA nucleotides

A DNA string is built from an alphabet of four symbols/nucleotides: 'A', 'C', 'G' and 'T'.
Input: path to a FASTA file. Output: four integers, the count of A, C, G, T (case-insensitive).

```python
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
```

Output on `data/Q1/input1.fa`: A=13 C=13 G=17 T=17
Output on `data/Q1/chr22.fa`: A=10382214 C=9160652 G=9246186 T=10370725

## Question 2: Finding a motif in DNA

Input: path to a FASTA file and a motif string s. Output: count of occurrences of s
as a substring of the sequence (forward strand only).

```python
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
```

Output on `{path: "data/Q1/input1.fa", s: "CGTAACC"}`: 4
Output on `{path: "data/Q1/chr22.fa", s: "CGTAACC"}`: 206

## Question 3: Consensus and profile matrix

Input: a FASTA file of equal-length aligned DNA sequences. Output: the consensus
string, and the 3 positions (1-based) where the sequences vary the most.

I used entropy to score each column: for each of A/C/G/T I work out what
fraction of the 25 sequences have that base, then plug those fractions into
H = -Σ p·log2(p). If everyone at a column agrees, entropy is 0. If all four
bases show up equally often, entropy hits its max of 2 bits. I went with
this instead of just counting how many sequences disagree with the majority
base, because it also picks up on how that disagreement is spread out — a
column where the minority is split across two or three different bases
feels more "confused" to me than one with a single, consistent alternative
base, even when the number of mismatches is the same.

```python
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

# consensus + entropy per column
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
```

Output on `data/Q3/BRCA_aligned.fa`:
```
GATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGAGATGTGGTCAATGGAAGAAACCACCAAGGTCCAAAGCGAGCAAGAGAATCCCAGGACAGAAAGGTAAAGCTCCCTCCCTCAAGTTGACAAAAATCTCACCCCACCACTCTGTATTCCACTCCCCTTTGCAGAGATGGGCCGCTTCATTTTGTAAGACTTATTACATACATACACAGTGCTAGATACTTTCACACAGGTTCTTTTTTCACTCTTCCATCCCAACCACATAAATAAGTATTGTCTCTACTTTATGAATGATAAAACTAAGAGATTTAGAGAGGCTGTGTA
```
Top 3 least-conserved positions (1-based): 26, 158, 144

## Question 4: Multiple sequence alignment and covariation matrices

Input: a FASTA file of equal-length aligned protein sequences (20 amino acids + gap "-").
Output: the consensus string and the top 10 co-varying column pairs by mutual information,
MI(i,j) = sum_a sum_b f_ij(a,b) * log2( f_ij(a,b) / (f_i(a) * f_j(b)) ).

```python
from collections import Counter
from math import log2

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


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


def clean(ch):
    return ch if ch in AMINO_ACIDS else "-"


seqs = read_fasta("data/Q4/DHFR_aligned.fa")
n = len(seqs)
length = len(seqs[0])

columns = [[clean(seq[i]) for seq in seqs] for i in range(length)]
freqs = [Counter(col) for col in columns]
for f in freqs:
    for k in f:
        f[k] /= n

consensus = "".join(c.most_common(1)[0][0] for c in freqs)

results = []
for i in range(length):
    for j in range(i + 1, length):
        joint = Counter(zip(columns[i], columns[j]))
        mi = 0.0
        for (a, b), count in joint.items():
            pij = count / n
            mi += pij * log2(pij / (freqs[i][a] * freqs[j][b]))
        results.append((mi, i + 1, j + 1))

results.sort(reverse=True)

print("Consensus string:")
print(consensus)
print("Top 10 co-varying column pairs by mutual information:")
for mi, i, j in results[:10]:
    print(f"  columns ({i}, {j}): MI = {mi:.4f} bits")
```

Output on `data/Q4/DHFR_aligned.fa`:
```
MISLIVAMAENGVIGKDNDLPWHLPEDLKYFKRLTLGKPVIMGRKTWESIGRPLPGRRNIVLSRDPDYQAEGVEVVHSLEEALALAG-VEEVFVIGGAEIYAQALP-ADRLYLTEIDAEFEGDTFFPEIDPDEWEEVSREEHPADEKNGYDYTFVTYERK
```

Top 10 pairs (columns, MI in bits):
(149,150) 1.2409, (150,151) 1.1677, (68,70) 1.1111, (148,150) 1.0829,
(145,146) 1.0511, (13,122) 1.0165, (149,151) 1.0131, (6,8) 1.0072,
(151,152) 0.9992, (58,74) 0.9862

## Question 5: A restriction-nuclease case study

EcoRV recognizes GATATC and cuts blunt in the middle: GAT/ATC. Target sequence: chr22
(`data/Q1/chr22.fa`).

```python
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
```

a. How many perfectly matching EcoRV sites are present in chr22?
   **5024**

b. GATATC's reverse complement is **GATATC** — same sequence, so it's a palindrome.
   I only need to scan the forward strand because of that: the two DNA strands run
   opposite ways but pair up base-for-base, so whatever sequence sits on the top
   strand at a spot is mirrored on the bottom strand at that same spot. Since the
   site folds back onto itself, every site I find on the forward strand is already
   the same site on the reverse strand too — scanning the reverse strand again
   wouldn't turn up anything new, just the same positions from the other side.

c. Mutant EcoRV cuts GAX/ATC (mismatch tolerated at position 3). Check string
   `TTGATATCAAGAGATCCTTCCGAAATCACGT` gives 1 exact + 2 relaxed sites, matching the
   assignment's example. On chr22: **31370** total sites (5024 exact + 26346 relaxed-only).

d. Cutting chr22 at every exact site: **5025** fragments, median length **3561 bp**
   (N's excluded). Distribution:

![fragment length distribution](results/q5d_fragment_length_hist.png)

## Question 6: ORF finder

ORF: start codon ATG to the first in-frame stop codon (TAA/TAG/TGA). Target: E. coli
genome (`data/Q6/ecoli_genome.fa`).

```python
import csv

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def find_orfs(seq):
    orfs = []
    for frame in range(3):
        open_starts = []
        for pos in range(frame, len(seq) - 2, 3):
            codon = seq[pos:pos + 3]
            if codon == "ATG":
                open_starts.append(pos)
            elif codon in STOP_CODONS:
                for start in open_starts:
                    protein = "".join(CODON_TABLE.get(seq[c:c + 3], "X") for c in range(start, pos, 3))
                    orfs.append((start, pos + 3, protein))
                open_starts = []
    return orfs


lines = []
for line in open("data/Q6/ecoli_genome.fa"):
    if not line.startswith(">"):
        lines.append(line.strip().upper())
genome = "".join(lines)
genome_len = len(genome)

orfs = []
for start, end, protein in find_orfs(genome):
    orfs.append(("+", start + 1, end, protein))

rc = revcomp(genome)
for start, end, protein in find_orfs(rc):
    orfs.append(("-", genome_len - start, genome_len - end + 1, protein))

long_orfs = [o for o in orfs if len(o[3]) > 100]
long_orfs.sort(key=lambda o: (o[0], o[1]))

with open("results/q6_orfs.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["strand", "id", "start", "end", "protein"])
    for i, (strand, start, end, protein) in enumerate(long_orfs, 1):
        writer.writerow([strand, i, start, end, protein])

print("total ORFs (both strands):", len(orfs))
print("ORFs with protein > 100 aa:", len(long_orfs))

strand, start, end, protein = max(orfs, key=lambda o: len(o[3]))
print("longest ORF:", strand, start, end, "-", len(protein), "aa")
print(protein)
```

Check: `ATGGCCATGGCGCCCAGAACTGGGCCCTGA` translates to `MAMAPRTGP` — matches.

a. Forward-strand ORFs longer than 100 aa are written to `results/q6_orfs.csv`
   (columns: strand, id, start, end, protein), submitted alongside this document.

b. Reverse-strand ORFs use the same scan on the reverse complement, with coordinates
   mapped back onto the original genome numbering, and are appended to the same CSV
   (labeled "-"). Combined total: 153321 ORFs found on both strands; 30894 of them
   translate to a protein longer than 100 aa (these are the ones in the CSV).

c. Longest ORF: forward strand, genome position 2044911-2052014, 2367 amino acids
   (7104 nt including the stop codon):
```
MLARSGKVSMATKKRSGEEINDRQILCGMGIKLRRLTAGICLITQLAFPMAAAAQGVVNAATQQPVPAQIAIANANTVPYTLGALESAQSVAERFGISVAELRKLNQFRTFARGFDNVRQGDELDVPAQVSEKKLTPPPGNSSDNLEQQIASTSQQIGSLLAEDMNSEQAANMARGWASSQASGAMTDWLSRFGTARITLGVDEDFSLKNSQFDFLHPWYETPDNLFFSQHTLHRTDERTQINNGLGWRHFTPTWMSGINFFFDHDLSRYHSRAGIGAEYWRDYLKLSSNGYLRLTNWRSAPELDNDYEARPANGWDVRAESWLPAWPHLGGKLVYEQYYGDEVALFDKDDRQSNPHAITAGLNYTPFPLMTFSAEQRQGKQGENDTRFAVDFTWQPGSAMQKQLDPNEVAARRSLAGSRYDLVDRNNNIVLEYRKKELVRLTLTDPVTGKSGEVKSLVSSLQTKYALKGYNVEATALEAAGGKVVTTGKDILVTLPAYRFTSTPETDNTWPIEVTAEDVKGNLSNREQSMVVVQAPTLSQKDSSVSLSTQTLNADSHSTATLTFIAHDAAGNPVVGLVLSTRHEGVQDITLSDWKDNGDGSYTQILTTGAMSGTLTLMPQLNGVDAAKAPAVVNIISVSSSRTHSSIKIDKDRYLSGNPIEVTVELRDENDKPVKEQKQQLNNAVSIDNVKPGVTTDWKETADGVYKATYTAYTKGSGLTAKLLMQNWNEDLHTAGFIIDANPQSAKIATLSASNNGVLANENAANTVSVNVADEGSNPINDHTVTFAVLSGSATSFNNQNTAKTDVNGLATFDLKSSKQEDNTVEVTLENGVKQTLIVSFVGDSSTAQVDLQKSKNEVVADGNDSVTMTATVRDAKGNLLNDVMVTFNVNSAEAKLSQTEVNSHDGIATATLTSLKNGDYRVTASVSSGSQANQQVNFIGDQSTAALTLSVPSGDITVTNTAPQYMTATLQDKNGNPLKDKEITFSVPNDVASKFSISNGGKGMTDSNGVAIASLTGTLAGTHMIMARLANSNVSDAQPMTFVADKDRAVVVLQTSKAEIIGNGVDETTLTATVKDPSNHPVAGITVNFTMPQDVAANFTLENNGIAITQANGEAHVTLKGKKAGTHTVTATLGNNNTSDSQPVTFVADKASAQVVLQISKDEITGNGVDSATLTATVKDQFDNEVNNLPVTFSSASSGLTLTPGVSNTNESGIAQATLAGVAFGEKTVTASLANNGASDNKTVHFIGDTAAAKIIELAPVPDSIIAGTPQNSSGSVITATVVDNNGFPVKGVTVNFTSNAATAEMTNGGQAVTNEQGKATVTYTNTRSSIESGARPDTVEASLENGSSTLSTSINVNADASTAHLTLLQALFDTVSAGETTSLYIEVKDNYGNGVPQQEVTLSVSPSEGVTPSNNAIYTTNHDGNFYASFTATKAGVYQLTATLENGDSMQQTVTYVPNVANAEITLAASKDPVIADNNDLTTLTATVADTEGNAIANTEVTFTLPEDVKANFTLSDGGKVITDAEGKAKVTLKGTKAGAHTVTASMTGGKSEQLVVNFIADTLTAQVNLNVTEDNFIANNVGMTRLQATVTDGNGNPLANEAVTFTLPADVSASFTLGQGGSAITDINGKAEVTLSGTKSGTYPVTVSVNNYGVSDTKQVTLIADAGTAKLASLTSVYSFVVSTTEGATMTASVTDANGNPVEGIKVNFRGTSVTLSSTSVETDDRGFAEILVTSTEVGLKTVSASLADKPTEVISRLLNASADVNSATITSLEIPEGQVMVAQDVAVKAHVNDQFGNPVAHQPVTFSAEPSSQMIISQNTVSTNTQGVAEVTMTPERNGSYMVKASLPNGASLEKQLEAIDEKLTLTASSPLIGVYAPTGATLTATLTSANGTPVEGQVINFSVTPEGATLSGGKVRTNSSGQAPVVLTSNKVGTYTVTASFHNGVTIQTQTTVKVTGNSSTAHVASFIADPSTIAATNTDLSTLKATVEDGSGNLIEGLTVYFALKSGSATLTSLTAVTDQNGIATTSVKGAMTGSVTVSAVTTAGGMQTVDITLVAGPADTSQSVLKSNRSSLKGDYTDSAELRLVLHDISGNPIKVSEGMEFVQSGTNVPYIKISAIDYSLNINGDYKATVTGGGEGIATLIPVLNGVHQAGLSTTIQFTRAEDKIMSGTVSVNGTDLPTTTFPSQGFTGAYYQLNNDNFAPGKTAADYEFSSSASWVDVDATGKVTFKNVGSNSERITATPKSGGPSYVYEIRVKSWWVNAGEAFMIYSLAENFCSSNGYTLPRANYLNHCSSRGIGSLYSEWGDMGHYTTDAGFQSNMYWSSSPANSSEQYVVSLATGDQSVFEKLGFAYATCYKNL
```

I think this is the `yeeJ` gene. When I looked it up, MG1655's YeeJ protein is
listed at 2358 amino acids, which is really close to the 2367 aa I got — the
small gap is probably just a different genome build or annotation call, not a
different gene. My translated sequence also has a lot of repeated ~90-residue
chunks, which lines up with YeeJ being described as having around 13 repeated
Ig-like domains. As for what it does: it's an inverse autotransporter that
sticks to peptidoglycan and helps E. coli form biofilms (Meuskens et al.,
Scientific Reports, 2017, https://www.nature.com/articles/s41598-017-10902-0).
