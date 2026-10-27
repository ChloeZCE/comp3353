# COMP3353 Bioinformatics — Assignment 2: Sequence alignment and variant calling

**Name:** Zheng Choi I\
**University Number:** 3035987788\
**Email:** u3598778@connect.hku.hk

All code is Python 3 (standard library; `scipy` only for Q6) or bash, and is
also in the `src/` folder of this submission. Every Python script reads its
input path from the command line, as the questions require. `src/fasta.py` is a
small FASTA reader that Q2–Q5 import:

```python
@@src/fasta.py@@
```

`bash src/run_all.sh` (run from `assignment2/`) re-runs Q2–Q6 and writes every
output to `results/`.


---

## Question 1 — Sequence alignment with bwa (ERR2337147, *S. aureus*)

Input files (in `data/Q1/`): `ref.fna.gz`, the *S. aureus* NCTC 8325
chromosome (NC_007795.1, 2,821,361 bp), and `ERR2337147_1.fastq.gz` /
`ERR2337147_2.fastq.gz` (184,803 read pairs of 250 bp, Illumina MiSeq).
Tools used: bwa 0.7.18 and samtools 1.24.

### 1(a) Commands and igv visualisation

Install the tools (same conda environment as tutorial 2):

```bash
conda create -n comp3353 -c bioconda -c conda-forge bwa samtools igv
conda activate comp3353
```

I ran every step with `src/run_q1.sh` (run from `assignment2/`), and the full
output is in `results/q1_output.txt`:

```bash
@@src/run_q1.sh@@
```

The steps:

1. `gunzip -k ref.fna.gz`, then `bwa index ref.fna` builds the FM-index of the
   reference (`.amb .ann .bwt .pac .sa`). `samtools faidx ref.fna` writes the
   `.fai` that igv needs.
2. `bwa mem -t 4 ref.fna R1 R2 | samtools view -bS - > ERR2337147.bam` aligns
   the two FASTQ files as pairs and saves the output as BAM.
3. `samtools sort … -o ERR2337147.sorted.bam` orders the records by coordinate.
   `samtools index` writes the `.bai` used for region lookups and by igv.
4. `igv -g data/Q1/ref.fna results/q1/ERR2337147.sorted.bam` loads the genome
   and the alignments. I then typed the region below into the locus box.

`samtools flagstat` on the sorted BAM:

```
@@results/q1_output.txt:1-16@@
```

94.9% of the reads mapped and 93.5% are properly paired. Mean depth over the
whole genome is 30.2× (`samtools depth -a`). bwa estimated the insert size
at 364 ± 80 bp (FR orientation).

**Region shown: `NC_007795.1:2,041,801-2,042,200`** (400 bp). I wanted a
region with good coverage and no ambiguous mappings. The 2 kb window
`NC_007795.1:2,041,001-2,043,000` has mean depth 65.5× (min 43×, max 81×).
It contains 579 reads, 577 of them with MAPQ 60, and no read with MAPQ < 10.

![igv screenshot of NC_007795.1:2,041,801-2,042,200](results/q1_igv.png)

(This screenshot is from igv.js 3.1.2, the web version of IGV, showing the same
sorted BAM and reference; `tools/igv_snapshot/` has the page used. Desktop IGV
shows the same view after `igv -g data/Q1/ref.fna results/q1/ERR2337147.sorted.bam`
with `NC_007795.1:2,041,801-2,042,200` typed in the locus box.)

The top track is the reference sequence. Below it are the coverage histogram
(about 70×) and the piled-up reads. Each grey bar is one read, and the arrow
end shows its strand. Coloured ticks are bases that differ from the reference.
Several columns are coloured in *every* read covering them. Those are not
sequencing errors. They are real differences (SNPs) between the sequenced
isolate and the NCTC 8325 reference strain. Isolated single ticks in one read
are sequencing errors.

### 1(b) FLAG of the 10th record

```
$ samtools view results/q1/ERR2337147.sorted.bam | sed -n 10p | cut -f1-9
ERR2337147.11	163	NC_007795.1	1	60	97S153M	=	105	352
```

Decoded with `src/q1_decode_flag.py`:

```python
@@src/q1_decode_flag.py@@
```

```
@@results/q1_output.txt:23-36@@
```

**FLAG = 163 = 128 + 32 + 2 + 1** (binary `0000 1010 0011`):

| bit (hex) | decimal | name | set? | meaning for this read |
|---|---|---|---|---|
| 0x1 | 1 | PAIRED | **yes** | the read comes from a paired-end template |
| 0x2 | 2 | PROPER_PAIR | **yes** | both mates aligned in the expected FR orientation at a plausible distance (TLEN 352, close to the 364 bp mean insert) |
| 0x4 | 4 | UNMAP | no | the read itself is mapped |
| 0x8 | 8 | MUNMAP | no | the mate is mapped |
| 0x10 | 16 | REVERSE | no | the read aligned to the forward strand |
| 0x20 | 32 | MREVERSE | **yes** | the mate aligned to the reverse strand |
| 0x40 | 64 | READ1 | no | not the first read |
| 0x80 | 128 | READ2 | **yes** | the second read of the pair (from `ERR2337147_2.fastq.gz`) |
| 0x100 | 256 | SECONDARY | no | primary alignment |
| 0x200 | 512 | QCFAIL | no | passed QC |
| 0x400 | 1024 | DUP | no | not marked as a duplicate |
| 0x800 | 2048 | SUPPLEMENTARY | no | not a supplementary (split) alignment |

**Answer:** this is the **second read** of the pair (0x80 set, 0x40 not set),
and its **mate is mapped** (0x8 not set). The other fields agree: RNEXT `=`,
PNEXT 105 (the mate starts at position 105), and the pair is properly paired
on opposite strands.

Its CIGAR `97S153M` and POS 1 are also worth noting. The read starts at the
very first base of the chromosome, and its first 97 bases are soft-clipped.
The *S. aureus* chromosome is circular, so those 97 bases come from the end of
the linear reference sequence, which the aligner cannot wrap around to.

### 1(c) Reads with MAPQ < 10

```
@@results/q1_output.txt:37-60@@
```

All 10 examples are in `results/q1_output.txt`. More counts:

| MAPQ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| records | 1937 | 118 | 161 | 152 | 117 | 165 | 163 | 167 | 94 | 195 |

Of the 355,435 mapped records, 3,269 (0.9%) have MAPQ < 10. Their CIGAR
strings fall into two groups:

1. **Clipped alignments: 2,602 of 3,269 (80%)** have soft (`S`) or hard (`H`)
   clipping, and many also have indels (`I`/`D`). Only part of the read
   aligns. 1,850 of them carry an `SA:Z` tag, so the read is split across two
   places (chimeric or spanning a structural difference from the reference).
   The `…M…H` CIGARs such as `65M186H` and `47M204H` in the list above are the
   short supplementary pieces of those split reads; bwa hard-clips the rest of
   the read in supplementary records. A short aligned piece scores
   low, so its score is close to the next-best hit, and MAPQ drops.
2. **Full-length matches: 585 records** are plain `M` with no clipping and no
   indels (`251M`, `250M`, `249M` …). The whole read aligns, but 574 of them have an `XA:Z` tag listing
   another equally good or nearly equally good location, and 465 have MAPQ 0.
   These reads come from **repeats**: rRNA operons, insertion
   sequences/transposons and duplicated genes, where the aligner cannot tell
   which copy the read came from. MAPQ is the aligner's confidence in the
   *position* of a read, not a measure of how well the bases match. So a
   perfect `250M` alignment can still get MAPQ 0.

In short: low MAPQ comes either from a read that only partly aligns (clipped
CIGAR) or from a read that aligns equally well in more than one place
(repeat; full-length `M` plus an `XA` tag).

**Example read with MAPQ in [1, 5]:**

```
ERR2337147.5185  FLAG 185  NC_007795.1:36551  MAPQ 5  CIGAR 9S240M
NM:i:6  AS:i:210  XS:i:205  XA:Z:NC_007795.1,-42681,9S240M,7;
```

* `9S` means the first 9 bases of the read are **soft-clipped**. They stay in
  SEQ but are not aligned, most likely low-quality bases or bases that do not
  match the reference here.
* `240M` means the next 240 bases align to the reference starting at POS 36,551.
  `M` covers matches *and* mismatches. `NM:i:6` says 6 of the 240 bases differ
  from the reference, and there are no insertions or deletions.
* The read is 9 + 240 = 249 bp long and covers reference positions
  36,551–36,790.

The read aligns in full, so the CIGAR alone does not explain the low MAPQ.
The tags do. The best alignment scores `AS = 210`, and the second-best hit
scores `XS = 205`. That second hit is listed in `XA`: position 42,681 on the
reverse strand, with the *same* CIGAR `9S240M` and 7 mismatches instead of 6.
This stretch of sequence occurs twice in the genome, about 6 kb apart, and the
read fits both copies almost equally well; the best copy wins by just one
mismatch. bwa therefore gives it the primary position but only MAPQ 5, about
a 30% chance (10^(−5/10) ≈ 0.32) that the position is wrong. FLAG 185 also
shows that its mate is unmapped (0x8), so there is no proper pair to help pick
the right copy.

---

## Question 2 — Hamming distance

### Source code (`src/q2_hamming.py`)

```python
@@src/q2_hamming.py@@
```

The script compares the two sequences position by position and counts the
mismatches. It raises an error if the lengths differ. Check: the Rosalind
example (`GAGCCTACTAACGGGAT` / `CATCGTAATGACGGCCT`) gives 7.

### Output on `data/Q2/input.fa`

```
$ python3 src/q2_hamming.py data/Q2/input.fa
79
```

**Answer: the Hamming distance is 79.** Both sequences (`Sequence_1`, `Sequence_2`) are 500 bp long, and they differ at 79 of the 500 positions.

---

## Question 3 — Edit distance

### Source code (`src/q3_edit_distance.py`)

```python
@@src/q3_edit_distance.py@@
```

### How the distance matrix is populated

`D` has `len(s)+1` rows and `len(t)+1` columns. `D(i, j)` is the edit distance
between the first `i` symbols of `s` and the first `j` symbols of `t`.

1. **Initial values.** `D(0,0) = 0`. `D(i,0) = i`: turning `s[:i]` into the empty
   string takes `i` deletions. `D(0,j) = j`: building `t[:j]` from nothing takes
   `j` insertions.
2. **Recurrence.** I fill row by row, left to right. That way, when I compute
   cell `(i, j)`, the three cells it depends on (up, left and up-left) are
   already filled:

   ```
   D(i,j) = min( D(i-1,j)   + 1,                  # delete s_i
                 D(i,j-1)   + 1,                  # insert t_j
                 D(i-1,j-1) + [s_i != t_j] )      # match (cost 0) or substitute (cost 1)
   ```
   The optimal alignment of `s[:i]` and `t[:j]` must end in one of these three
   ways, so the smallest of the three costs is the answer for that cell.
3. **Result.** The edit distance is `D(len(s), len(t))`, the bottom-right cell.

Row `i` only depends on row `i-1`, so the code keeps just two rows in memory.
Run time is still O(|s|·|t|), but memory drops to O(|t|). That matters for
long input sequences. Check: `PLEASANTLY`/`MEANLY` gives 5, the textbook value.

### Output

```
$ python3 src/q3_edit_distance.py data/Q2/input.fa
7
$ python3 src/q3_edit_distance.py data/Q3/input.fa
1704
```

**Answer for `data/Q3/input.fa`: the edit distance is 1704.** The file holds
two hepatitis B virus genomes: X02763.1 (serotype adw2, 3,221 bp) and
NC_003977.2 (strain ayw, 3,182 bp). Computing it fills a 3,222 × 3,183 matrix,
about 10 million cells, which takes about 1.5 s.

1704 is much higher than the true difference between the two viruses. The HBV
genome is circular, and the two records start at different points of the
circle: the first 30 bases of X02763.1 appear at position 1,815 of
NC_003977.2. Edit distance compares the strings as linear sequences, so it has
to delete a large block at one end and insert it again at the other end. When
I rotate NC_003977.2 so that both start at the same base, the edit distance
drops to **335**, about 10% of the genome. That is roughly what you would
expect between two different HBV strains. The 39 bp length difference
(3,221 vs 3,182) is part of those 335 edits.

**Comparison with Q2:** On `data/Q2/input.fa` the edit distance is **7**, but
the Hamming distance is **79**. The Hamming distance is one valid way to edit
`s` into `t`: substitutions only, no gaps. So the edit distance can never be
larger (`edit ≤ Hamming`). Here it is much smaller because the two sequences
differ by a few insertions and deletions, which Hamming distance cannot model.
One optimal set of 7 edits (positions in `Sequence_1`, 1-based):

| position in s | edit | cost |
|---|---|---|
| 159 | substitute C → G | 1 |
| after 200 | insert `GC` | 2 |
| 301–302 | delete `AA` | 2 |
| 435 | substitute C → A | 1 |
| 440 | substitute A → G | 1 |
| | **total** | **7** |

The 2-base insertion after position 200 shifts every later base of `t` two
places to the right. The 2-base deletion at 301–302 shifts them back. Between
the two, about 100 positions are compared out of register. Hamming distance
counts most of them as mismatches: 79 = the 3 real substitutions plus 76
mismatches caused by the shift. Edit distance aligns around the gap and
charges only 2 + 2 for the indels.

---

## Question 4 — BWT encoder and decoder

### Source code (`src/bwt.py`, encoder and decoder)

```python
@@src/bwt.py@@
```

### 4(a) Encoder

I append `$` (sorts before A, C, G, T) to the sequence. The BWT is the last
column of the sorted list of all rotations. Sorting rotations of a
`$`-terminated string is the same as sorting its suffixes, so the code builds
a suffix array by prefix doubling. That needs O(n) memory, not the n² of a full
rotation matrix, so whole-genome inputs work. The last-column character of
the rotation starting at `i` is `text[i-1]`. Check: `BANANA` → `ANNB$AA`.
Random DNA strings match a brute-force sort of all rotations.

```
$ python3 src/bwt.py encode data/Q4/input1.fa
TC$GAATCCGG
$ python3 src/bwt.py encode data/Q4/input2.fa | cut -c1-50
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTGTCTTCTTTCCTTGG
```

**input1** (`ACGTGCACGT`, 10 bp). The sorted rotations of `ACGTGCACGT$` are:

| row | rotation | last column |
|---|---|---|
| 0 | `$ACGTGCACGT` | T |
| 1 | `ACGT$ACGTGC` | C |
| 2 | `ACGTGCACGT$` | $ |
| 3 | `CACGT$ACGTG` | G |
| 4 | `CGT$ACGTGCA` | A |
| 5 | `CGTGCACGT$A` | A |
| 6 | `GCACGT$ACGT` | T |
| 7 | `GT$ACGTGCAC` | C |
| 8 | `GTGCACGT$AC` | C |
| 9 | `T$ACGTGCACG` | G |
| 10 | `TGCACGT$ACG` | G |

Reading the last column top to bottom gives **BWT = `TC$GAATCCGG`**.

**input2** is the SARS-CoV-2 reference genome (NC_045512.2, 29,903 bp). The
first 50 BWT symbols are `AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTGTCTTCTTTCCTTGG`.

**What I find:** the BWT starts with a run of **33 A's**. The genome ends with
a 33-nt poly-A tail (`…GAATGAC` followed by 33 `A`). In the sorted rotation
matrix, row 0 is the rotation starting with `$`, and its last character is the
final base of the genome, an A. The rows right after it start with `A$`,
`AA$`, `AAA$`, …, i.e. the suffixes inside the poly-A tail. Each of those is
preceded by another A, so all of them put an A in the last column. Row 33
(`AAA…A$`, the whole tail) is preceded by the C before the tail, which is
where the run ends.

This is the general property of the BWT: rotations that start with the same
context sort next to each other, and the same context tends to be preceded by
the same base. So the BWT groups identical symbols into runs. Over the whole
input2 BWT there are 21,518 runs for 29,904 symbols (longest run 33). This is
why the BWT compresses well with run-length encoding. The BWT holds exactly
the same symbols as the input plus one `$`, only reordered, which is what
makes it reversible (4b) and searchable (Q5).

### 4(b) Decoder

The decoder uses the **LF mapping**. The first column `F` is `sorted(L)`. The
k-th occurrence of character `c` in `L` and the k-th occurrence of `c` in `F`
are the same character of the text. So `LF(i) = C[c] + (number of c in L[:i])`.
Row 0 of the sorted rotations starts with `$`, so `L[0]` is the last base of
the sequence. Applying `row = LF(row)` repeatedly walks the text backwards,
one base per step, and reversing the collected bases gives the original. It
runs in O(n).

```
$ python3 src/bwt.py decode data/Q4/input3.txt
AGTTGTTGATCTGTGTGAGTCAGACTGCGACAGTTCGAGTCTGAAGCGAGA
$ python3 src/bwt.py encode data/Q4/input1.fa > results/input1.bwt.txt
$ python3 src/bwt.py decode results/input1.bwt.txt
ACGTGCACGT
```

**input3:** the BWT string `AGGGGGCAGGC$GTAGGTTATATCACTCTTAAATTATAGTGACGTCGCTGGG`
(52 symbols) decodes to the 51 bp sequence
**`AGTTGTTGATCTGTGTGAGTCAGACTGCGACAGTTCGAGTCTGAAGCGAGA`**. Encoding that
sequence again gives back exactly the input3 string.

**Is it the same as the original?** **Yes.** Decoding `TC$GAATCCGG` gives
`ACGTGCACGT`, which is identical to the sequence in `data/Q4/input1.fa`
(`run_all.sh` prints `identical: True`).
The BWT can be inverted exactly because the `$` marker fixes where the text ends.
I also tested encode→decode on 300 random strings and one random 200 kb
string, and every one came back unchanged.

---

## Question 5 — Search with BWT (FM-index)

### Source code (`src/q5_fm_search.py`, reuses `bwt_encode` from Q4)

```python
@@src/q5_fm_search.py@@
```

### Method

* `C[c]` = number of symbols in `t$` that are smaller than `c`
  (`$ < A < C < G < T`). It is the row where the rows starting with `c` begin
  in the sorted rotation matrix.
* `Occ(c, k)` = number of times `c` appears in `t'[0..k-1]`, the first `k`
  symbols of the BWT. It is built as cumulative counts with `n+1` rows.
* **Backward search.** Start with the interval of all rows `[lo, hi) = [0, n)`
  and process the motif from its **last** character to its first:
  `lo ← C[c] + Occ(c, lo)`, `hi ← C[c] + Occ(c, hi)`. After every step,
  `[lo, hi)` holds exactly the rows whose rotations start with the suffix of
  the motif read so far. When the whole motif is processed, the number of
  occurrences is `hi − lo`. If the interval becomes empty, the count is 0.
  Each character costs O(1), so a search takes O(|s|) whatever the size of the text.

Check: the counts match a brute-force overlapping regex count on 300 random
strings. On `ACGTACGTTACGTAACG` the trace is `T:[14,18) → G:[11,14) → C:[7,10) → A:[3,6)`, which gives 3.

### Output on {data/Q4/input1.fa, s = "ACGT"} (BWT, C, Occ and search steps)

```
@@results/q5_output.txt:1-30@@
```

**Answer: `ACGT` occurs 2 times in input1.** By hand, the steps are:

1. Start with all rows `[0, 11)`.
2. `T`: `lo = C[T] + Occ(T,0) = 9 + 0 = 9`, `hi = C[T] + Occ(T,11) = 9 + 2 = 11` → `[9, 11)`, the rows starting with `T`.
3. `G`: `lo = 6 + Occ(G,9) = 6 + 1 = 7`, `hi = 6 + Occ(G,11) = 6 + 3 = 9` → `[7, 9)`, the rows starting with `GT`.
4. `C`: `lo = 3 + Occ(C,7) = 3 + 1 = 4`, `hi = 3 + Occ(C,9) = 3 + 3 = 6` → `[4, 6)`, the rows starting with `CGT`.
5. `A`: `lo = 1 + Occ(A,4) = 1 + 0 = 1`, `hi = 1 + Occ(A,6) = 1 + 2 = 3` → `[1, 3)`, the rows starting with `ACGT`.

`hi − lo = 3 − 1 = 2`. Rows 1 and 2 in the rotation table of Q4(a) are
`ACGT$ACGTGC` and `ACGTGCACGT$`, i.e. the matches at positions 7 and 1 of
`ACGTGCACGT`.

### Output on {data/Q4/input2.fa, s = "ACGT"}

```
@@results/q5_output.txt:32-45@@
   ...   (rows 5-29 are all A as well; full output in results/q5_output.txt)
@@results/q5_output.txt:71-79@@
```

(For the long input2 the script prints only the BWT's first 200 symbols and
the first 30 rows of the `Occ` matrix plus its final row. The full matrix has
29,905 rows.)

**Answer: `ACGT` occurs 64 times in the SARS-CoV-2 genome.** A direct scan
of the sequence for overlapping matches also gives 64.

---

## Question 6 — Variant calling with Bayes' theorem

### Model

The reference base is A, and a read reports T. Each read is independent and
has a 1% misclassification rate `e`. The probability that one read shows T
under each genotype `G`:

| G | P(read = T \| G) |
|---|---|
| AA | e = 0.01 (an A misread as T) |
| AT | ½(1−e) + ½e = 0.5 |
| TT | 1 − e = 0.99 |

The number of T reads `k` out of `n` follows a binomial distribution:
`P(D | G) = C(n,k) p_G^k (1−p_G)^(n−k)`. By Bayes' theorem:

```
P(G | D) = P(D | G) P(G) / Σ_G' P(D | G') P(G'),     P(mutant | D) = P(AT | D) + P(TT | D)
```

Prior for AA = 1 − P(AT) − P(TT).

### Source code (`src/q6_variant_calling.py`)

```python
@@src/q6_variant_calling.py@@
```

### Output

```
@@results/q6_output.txt@@
```

### 6(a) 4 A + 1 T out of 5 reads; priors AT = 2 %, TT = 0.1 %, AA = 97.9 %

| G | likelihood P(D\|G) | prior | likelihood × prior | posterior |
|---|---|---|---|---|
| AA | C(5,1)(0.01)(0.99)⁴ = 0.04803 | 0.979 | 0.04702 | **0.9377** |
| AT | C(5,1)(0.5)⁵ = 0.15625 | 0.02 | 0.003125 | **0.0623** |
| TT | C(5,1)(0.99)(0.01)⁴ = 4.95×10⁻⁸ | 0.001 | 4.95×10⁻¹¹ | ≈ 9.9×10⁻¹⁰ |

Evidence P(D) = 0.05015. **P(mutant) = P(AT) + P(TT) ≈ 0.0623 (6.2 %)**:
almost all of that is heterozygous, and P(TT) is about 10⁻⁹. One T in five reads is
three times more likely under AT than under AA (0.156 vs 0.048). But AA
is 49× more likely a priori, so the single T is best explained as a
sequencing error. The patient is most likely **AA** (93.8 %).

### 6(b) 24 A + 6 T out of 30 reads, same priors

| G | likelihood | likelihood × prior | posterior |
|---|---|---|---|
| AA | C(30,6)(0.01)⁶(0.99)²⁴ = 4.665×10⁻⁷ | 4.567×10⁻⁷ | 0.0397 |
| AT | C(30,6)(0.5)³⁰ = 5.530×10⁻⁴ | 1.106×10⁻⁵ | **0.9603** |
| TT | C(30,6)(0.99)⁶(0.01)²⁴ = 5.59×10⁻⁴³ | 5.59×10⁻⁴⁶ | ≈ 5×10⁻⁴¹ |

**P(mutant) ≈ 0.960 (96.0 %)**, essentially all of it heterozygous (AT).
With 30 reads, six T's are very unlikely to come from 1% errors: we would
expect only 0.3 erroneous T's. The likelihood ratio AT:AA is now about 1,185.
That is enough to outweigh the 49:1 prior in favour of AA.

### 6(c) Same 24 A + 6 T, but a rarely-mutated site: AT = 0.05 %, TT = 0.001 %, AA = 99.949 %

| G | likelihood | prior | likelihood × prior | posterior |
|---|---|---|---|---|
| AA | 4.665×10⁻⁷ | 0.99949 | 4.663×10⁻⁷ | 0.6278 |
| AT | 5.530×10⁻⁴ | 0.0005 | 2.765×10⁻⁷ | **0.3722** |
| TT | 5.59×10⁻⁴³ | 0.00001 | 5.59×10⁻⁴⁸ | ≈ 7.5×10⁻⁴² |

**P(mutant) ≈ 0.372 (37.2 %)**. The data are exactly the same as in (b), but
the prior odds against a heterozygote are now about 2,000:1, not 49:1. That
outweighs the likelihood ratio of about 1,185, so AA (error-only) is
now the more probable explanation. The prior matters a lot for rare variants:
the same read evidence can support a confident variant call at one site and
a likely error at another. To call a variant here with confidence we would
need more depth or higher-quality reads.
