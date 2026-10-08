# COMP3353 Bioinformatics — Lab Session 2: Short-read alignment with BWA & samtools

**Name:** Zheng Choi I  **University Number:** 3035987788

Everything (download → `bwa index` → `bwa mem | samtools view -bS` → sort → index →
Q1–Q6, plus depth/variant numbers for Q8–Q9) is scripted in
[`src/run_lab02.sh`](src/run_lab02.sh):

```bash
conda activate comp3353
bash lab02/src/run_lab02.sh ~/labsession2 | tee lab02/results/lab02_output.txt
```

> Fill each `____` below with the number printed by the script (the matching
> `== Qn` block in `results/lab02_output.txt`).

Setup notes: with samtools ≥ 1.x the sort command is
`samtools sort SRR12351819.bam -o SRR12351819.sorted.bam` (the old
`samtools sort in.bam prefix` form is only for 0.1.19).

---

## Question 1 — number of alignments in `SRR12351819.bam`

```bash
samtools view -c SRR12351819.bam
```

`-c` counts records instead of printing them; every line of the BAM is one
alignment record. **Answer: ____ alignments.**
(The script also breaks this down into mapped / primary / secondary /
supplementary records — bwa mem can emit supplementary (0x800) records for
chimeric reads, so the total may be slightly larger than 2 × number of read pairs.)

## Question 2 — alignments overlapping NC_045512.2:11075-11102

```bash
samtools view -c SRR12351819.sorted.bam NC_045512.2:11075-11102
```

Region queries need the sorted **and indexed** BAM (`.bai`). samtools returns
every alignment whose span overlaps the region, not just those fully inside it.
**Answer: ____ alignments.**

## Question 3 — TLEN (insert size) statistics

```bash
samtools view SRR12351819.bam | awk '$9 > 0 {print $9}' | sort -n | uniq -c \
  | awk 'NR==1{min=$2} {max=$2; if($1>best){best=$1; mode=$2}}
         END{print "min",min,"max",max,"mode",mode,"("best"x)"}'
```

**Why is TLEN positive for one mate and negative for the other?** TLEN is the
signed distance from the leftmost mapped base of the fragment to the rightmost
one. The mate that maps to the left (smaller coordinate, normally the forward
strand read) gets `+TLEN`; the mate on the right gets the same magnitude with a
minus sign, so that the sign tells you which end of the fragment that read is.
Taking only `TLEN > 0` counts each fragment exactly once (and drops 0 = mate
unmapped / different reference).

| statistic | value |
|---|---|
| min TLEN | ____ |
| max TLEN | ____ |
| most frequent TLEN | ____ |

(The max is usually an outlier from a mis-paired / chimeric read pair, far
above the typical library insert size.)

## Question 4 — fraction properly paired

```bash
samtools flagstat SRR12351819.bam
samtools view -c -f 1 SRR12351819.bam   # paired
samtools view -c -f 2 SRR12351819.bam   # properly paired (FLAG 0x2)
```

FLAG bit 0x2 is set by the aligner when both mates map in the expected
orientation (FR) at a distance consistent with the estimated insert size.
**Answer: ____ / ____ = ____ %** (matches the "properly paired" line of
`flagstat`).

## Question 5 — strand

```bash
samtools view -c -F 20 SRR12351819.bam        # mapped (no 0x4) & forward (no 0x10)
samtools view -c -f 16 -F 4 SRR12351819.bam   # mapped & reverse (0x10 set)
```

FLAG 0x10 = read reverse-complemented. Unmapped reads (0x4) are excluded
because they have no strand. **Positive: ____  Negative: ____** — roughly
equal, as expected: in FR paired-end sequencing each fragment yields one
forward and one reverse alignment.

## Question 6 — MAPQ ≥ 20

```bash
samtools view -c -q 20 SRR12351819.bam
```

MAPQ 20 ⇒ P(wrong placement) ≤ 10^(-20/10) = 1 %.
**Answer: ____ alignments, ____ / (Q1 total) = ____ %.**
Unmapped reads have MAPQ 0, so they are part of the excluded fraction.

## Question 7 — IGV screenshot (NC_045512.2:11010-11168)

1. IGV → Genomes → Download Hosted Genome… → *SARS-CoV-2 (NC_045512.2)*
   (or in igv.org/app: Genome → SARS-CoV-2).
2. File → Load from File… → `SRR12351819.sorted.bam` (the `.bai` must be in
   the same folder).
3. Type `NC_045512.2:11010-11168` in the locus box → Go.
4. Right-click the alignment track → **View as pairs** (optionally *Color
   alignments by → insert size and pair orientation*).
5. Take the screenshot (File → Save PNG Image…) and upload it to Moodle.

## Question 8 — end of ORF8

From the NC_045512.2 annotation, ORF8 (gene *ORF8*, protein YP_009724396.1)
spans **27894–28259**, so the **last position of ORF8 is 28259** (the last
base of its TAA stop codon). Type `ORF8` (or `NC_045512.2:28259`) in IGV's
locus box and click on the coverage track at 28259 — the popup shows
"Total count".

Number of aligned reads at 28259 (same as the script's
`samtools depth`/`samtools view -c -F 0x604` output): **____**.

## Question 9 — mutation in ORF3a

ORF3a spans **NC_045512.2:25393–26220**. Navigate to `ORF3a` in IGV and look
for a column where (almost) every read is coloured — a position where the
sample differs from the reference. The script finds such columns with
`samtools mpileup` over this region (`src/pileup_variants.py`).

**Location: ____  Reference allele: ____  Observed allele: ____ (____ % of reads)**

Expected result for this 2020 sample: the well-known ORF3a mutation
**G25563T** (codon 57 CAG→CAT, **Q57H**), i.e. reference **G**, observed
**T**. Confirm the position/allele against the script output and your IGV view
before submitting.
