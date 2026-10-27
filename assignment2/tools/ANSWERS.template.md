# COMP3353 Bioinformatics — Assignment 2: Sequence alignment and variant calling

**Name:** Zheng Choi I
**University Number:** 3035987788
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

> **TODO before submitting:** replace every ⟨…⟩ placeholder with the real output
> or screenshot. Run `bash src/run_q1.sh` for Q1 and `bash src/run_all.sh` for Q2–Q5.

---

## Question 1 — Sequence alignment with bwa (ERR2337147, *S. aureus*)

### 1(a) Commands and igv visualisation

Install the tools (same conda environment as tutorial 2):

```bash
conda create -n comp3353 -c bioconda -c conda-forge bwa samtools igv
conda activate comp3353
```

Input files (in `data/Q1/`): `ref.fna.gz` (*S. aureus* NCTC 8325 chromosome, NC_007795.1) and `ERR2337147_1/2.fastq.gz` (184,803 read pairs, 250 bp Illumina MiSeq).

Then index the reference, align the two FASTQ files, sort and index the
results, and open them in igv (`src/run_q1.sh`):

```bash
@@src/run_q1.sh@@
```

The commands in order:

1. `bwa index ref.fa` builds the FM-index (BWT + suffix array) of the reference
   (`.amb .ann .bwt .pac .sa`), and `samtools faidx` writes the `.fai` that igv needs.
2. `bwa mem ref.fa R1 R2 | samtools view -bS -` aligns the two files as a pair
   and saves the output as BAM.
3. `samtools sort` orders the records by coordinate, and `samtools index`
   writes the `.bai` index used for region lookups.
4. `igv -g ref.fa ERR2337147.sorted.bam` loads the genome and the alignments.

**Region shown:** ⟨contig:start-end⟩. I picked it from the five
best-covered 1 kb windows that the `samtools depth` step prints. Mean depth
there is ⟨…⟩×.

![igv screenshot](results/q1_igv.png)
⟨insert igv screenshot: coverage track + piled-up reads⟩

### 1(b) FLAG of the 10th record

```
⟨paste: samtools view results/q1/ERR2337147.sorted.bam | sed -n 10p⟩
```

Decoder (`src/q1_decode_flag.py`):

```python
@@src/q1_decode_flag.py@@
```

FLAG = ⟨value⟩. The FLAG is a sum of powers of two, and each bit is one yes/no property:

| bit (hex) | decimal | name | set? | meaning |
|---|---|---|---|---|
| 0x1 | 1 | PAIRED | ⟨ ⟩ | read comes from a paired-end template |
| 0x2 | 2 | PROPER_PAIR | ⟨ ⟩ | both mates aligned with the expected orientation and insert size |
| 0x4 | 4 | UNMAP | ⟨ ⟩ | this read is unmapped |
| 0x8 | 8 | MUNMAP | ⟨ ⟩ | the mate is unmapped |
| 0x10 | 16 | REVERSE | ⟨ ⟩ | this read aligned to the reverse strand |
| 0x20 | 32 | MREVERSE | ⟨ ⟩ | the mate aligned to the reverse strand |
| 0x40 | 64 | READ1 | ⟨ ⟩ | first read of the pair (from `_1.fastq`) |
| 0x80 | 128 | READ2 | ⟨ ⟩ | second read of the pair (from `_2.fastq`) |
| 0x100 | 256 | SECONDARY | ⟨ ⟩ | secondary alignment |
| 0x200 | 512 | QCFAIL | ⟨ ⟩ | failed vendor QC |
| 0x400 | 1024 | DUP | ⟨ ⟩ | PCR/optical duplicate |
| 0x800 | 2048 | SUPPLEMENTARY | ⟨ ⟩ | supplementary (split/chimeric) alignment |

How to read it: bit 0x40 set means **first** read, bit 0x80 set means
**second** read. Bit 0x8 set means the mate is **unmapped**, otherwise it is
mapped. For example, FLAG 99 = 1 + 2 + 32 + 64: paired, proper pair, mate on
the reverse strand, first read, mate mapped. FLAG 163 = 1 + 2 + 32 + 128 is the
same but for the second read.

**Answer:** ⟨first/second⟩ read in the pair; its mate is ⟨mapped/unmapped⟩.

### 1(c) Reads with MAPQ < 10

```
⟨paste the "(c)" block printed by run_q1.sh: count, CIGAR table, clipping fraction, examples⟩
```

**What I find:** ⟨check this against your output⟩ Low-MAPQ reads come in two kinds:

* **Full-length matches** (e.g. `151M`) with MAPQ 0 or close to 0. The read
  aligns perfectly, but just as well to another place in the genome. Examples
  are the rRNA operons (*S. aureus* has 5–6 copies), IS elements and
  transposons, and duplicated genes. bwa picks one of the equal hits at
  random. These reads usually have an `XA:Z` tag listing the other hits, and
  `XS` (second-best score) close to or equal to `AS`. MAPQ measures how sure the
  aligner is about the *position*, not how well the bases match. So a perfect
  `151M` can still have MAPQ 0.
* **Clipped or gapped alignments** (CIGARs with `S`/`H`, `I`/`D`). Only part
  of the read aligns, for example a read that crosses a repeat boundary, a
  chimera, adapter read-through, or reads at contig ends. Their primary/supplementary
  (`SA:Z`) alignments split the read. A shorter aligned part gives a lower score,
  which brings the score closer to the second-best hit and lowers MAPQ.

**Example read with MAPQ in [1, 5]:** ⟨QNAME, MAPQ = ⟨m⟩, CIGAR = ⟨e.g. `30S121M`⟩⟩.
CIGAR operations are read left to right: `M` = aligned base (match *or*
mismatch), `I` = insertion relative to the reference, `D` = deletion from the
read, `S` = soft-clip (bases stay in SEQ but are not aligned), `H` = hard-clip
(bases removed from SEQ). So `30S121M` means the first 30 bases of the read do
not align and the next 121 bases align from POS onward. Only part of the read
supports this position, and ⟨there is an XA/SA tag or XS ≈ AS⟩, so the
placement is ambiguous and MAPQ is low.

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
⟨…⟩
```

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
⟨…⟩
$ python3 src/q3_edit_distance.py data/Q3/input.fa
⟨…⟩
```

**Comparison with Q2:** On `data/Q2/input.fa` the edit distance is ⟨e⟩ and the
Hamming distance is ⟨h⟩. The Hamming distance is one valid way to edit `s` into
`t`: substitutions only, no gaps. So the edit distance is never larger:
`edit ≤ Hamming`. It is smaller when an insertion plus a deletion fixes a
shifted stretch more cheaply than substituting every shifted position. For
example, `ACGTACGT` vs `CGTACGTA` has Hamming distance 8 but edit distance 2.
⟨State whether the two values are equal or not for your file.⟩

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
⟨…⟩
$ python3 src/bwt.py encode data/Q4/input2.fa | cut -c1-50
⟨…⟩
```

**What I find:** ⟨confirm against your output⟩ The BWT contains the same
symbols as the input plus one `$`, only reordered. The reordering puts
identical characters next to each other in long runs (e.g. `AAAAAA`,
`TTTT`). This happens because rotations that start with the same context
(e.g. every `CGT…`) sort next to each other, and in real DNA they tend to be
preceded by the same base. The longer and more repetitive the sequence (input2),
the longer these runs get. That is why the BWT compresses well with run-length
encoding, and why it is useful for search.

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
⟨…⟩
$ python3 src/bwt.py encode data/Q4/input1.fa > results/input1.bwt.txt
$ python3 src/bwt.py decode results/input1.bwt.txt
⟨…⟩
```

**Is it the same as the original?** ⟨Yes — `run_all.sh` prints `identical: True`⟩.
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
$ python3 src/q5_fm_search.py data/Q4/input1.fa ACGT --show-occ 1000
⟨…⟩
```

### Output on {data/Q4/input2.fa, s = "ACGT"}

```
$ python3 src/q5_fm_search.py data/Q4/input2.fa ACGT
⟨…⟩
```

(For the long input2 the script prints only the first 30 rows of the `Occ`
matrix plus its final row. The full matrix has n+1 rows.)

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
