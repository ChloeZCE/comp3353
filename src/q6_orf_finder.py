"""Q6: ORF finder for the E. coli genome.

Usage:
    python src/q6_orf_finder.py <genome_fasta_path> [--outdir results]

Writes results/q6_orfs.csv with columns: strand,id,start,end,protein
(only ORFs whose translated protein is longer than 100 amino acids, on both
strands), and prints the single longest ORF across both strands.

ORF definition used: within a single reading frame, scan codon by codon.
Every ATG start codon begins a candidate ORF; if a run of codons (between two
in-frame stop codons, or between a stop codon and the end of the frame)
contains one or more ATGs, each of those ATGs starts its own ORF that all
share the *same* next in-frame stop codon (a stop always closes every
still-open ORF in that frame). An ATG with no in-frame stop codon before the
end of the genome is an incomplete ORF and is excluded, since it has no
translated stop position.

Coordinate convention: positions are always reported in the original genome's
1-based forward-strand numbering, regardless of which strand the ORF is on.
For a "+"-strand ORF, start < end (5' -> 3' matches increasing coordinates).
For a "-"-strand ORF, start > end: "start" is the genomic position of the
first base of the ATG as read 5'->3' on the minus strand (i.e. the *higher*
forward-strand coordinate), and "end" is the last base of the stop codon
(the *lower* forward-strand coordinate) -- this keeps "start"/"end" meaning
"where translation begins/ends" on each ORF's own strand.
"""
import sys
import os
import csv
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from common import read_fasta_records, reverse_complement, translate_codon, STOP_CODONS

MIN_PROTEIN_LEN = 100


def find_orfs_in_frame(seq, frame):
    """Yield (start0, end0_exclusive, protein) for ORFs in this 0-based frame,
    where start0/end0 are 0-based positions into `seq` (the strand-local
    sequence), end0_exclusive is one past the last base of the stop codon."""
    n = len(seq)
    codon_starts = range(frame, n - 2, 3)
    active_starts = []  # 0-based positions (into seq) of open ATGs
    for pos in codon_starts:
        codon = seq[pos:pos + 3]
        if codon == "ATG":
            active_starts.append(pos)
        elif codon in STOP_CODONS:
            if active_starts:
                stop_end = pos + 3
                for start in active_starts:
                    protein_codons = range(start, pos, 3)
                    protein = "".join(translate_codon(seq[c:c + 3]) for c in protein_codons)
                    yield start, stop_end, protein
                active_starts = []
    # any positions left in active_starts have no downstream in-frame stop -> discarded


def find_orfs_one_strand(seq):
    for frame in range(3):
        yield from find_orfs_in_frame(seq, frame)


def build_orf_table(genome_seq):
    genome_len = len(genome_seq)
    rows = []  # dict: strand, start(1-based genome coord), end, protein, length

    # Forward strand
    for start0, end0, protein in find_orfs_one_strand(genome_seq):
        rows.append({
            "strand": "+",
            "start": start0 + 1,
            "end": end0,
            "protein": protein,
        })

    # Reverse strand: work on the reverse complement, then map coords back
    rc_seq = reverse_complement(genome_seq)
    for start0, end0, protein in find_orfs_one_strand(rc_seq):
        # rc index i (0-based) <-> genome 1-based position (genome_len - i)
        genome_start = genome_len - start0        # position of ATG's first base, minus strand
        genome_end = genome_len - end0 + 1         # position of stop codon's last base
        rows.append({
            "strand": "-",
            "start": genome_start,
            "end": genome_end,
            "protein": protein,
        })

    return rows


def write_csv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["strand", "id", "start", "end", "protein"])
        for i, row in enumerate(rows, start=1):
            writer.writerow([row["strand"], i, row["start"], row["end"], row["protein"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_path")
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    records = list(read_fasta_records(args.fasta_path))
    header, genome_seq = records[0]
    genome_seq = genome_seq.upper()
    print("Genome:", header, "length:", len(genome_seq))

    all_rows = build_orf_table(genome_seq)
    print("Total ORFs found (both strands, complete only):", len(all_rows))

    long_rows = [r for r in all_rows if len(r["protein"]) > MIN_PROTEIN_LEN]
    long_rows.sort(key=lambda r: (0 if r["strand"] == "+" else 1, r["start"]))
    print(f"ORFs with protein length > {MIN_PROTEIN_LEN} aa:", len(long_rows))

    out_csv = os.path.join(args.outdir, "q6_orfs.csv")
    write_csv(long_rows, out_csv)
    print("Wrote:", out_csv)

    longest = max(all_rows, key=lambda r: len(r["protein"]))
    print("\n=== (c) Longest ORF across both strands ===")
    print("Strand:", longest["strand"])
    print("Genome start:", longest["start"], "end:", longest["end"])
    nt_length = len(longest["protein"]) * 3 + 3  # + stop codon
    print("Protein length (aa, excl. stop):", len(longest["protein"]))
    print("Full ORF length (nt, incl. stop codon):", nt_length)
    print("Protein sequence:")
    print(longest["protein"])
