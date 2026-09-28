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
    # ATG -> next in-frame stop codon, per frame. Every ATG before a stop
    # gets its own ORF ending at that stop. ATGs with no stop are dropped.
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

# positions are reported in 1-based genome coordinates on the forward strand;
# for "-" strand ORFs, start/end are still those forward coordinates, just
# read from high to low (start = 5' end of the ORF on the minus strand)
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
