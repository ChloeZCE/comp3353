#!/usr/bin/env bash
# Q1: index, align, sort, index and inspect ERR2337147 (S. aureus) with bwa/samtools.
# Usage (from assignment2/):  bash src/run_q1.sh [ref.fa] [reads_1.fastq.gz] [reads_2.fastq.gz]
# Needs: conda install -c bioconda bwa samtools igv
set -euo pipefail

REF=${1:-data/Q1/reference.fasta}
R1=${2:-data/Q1/ERR2337147_1.fastq.gz}
R2=${3:-data/Q1/ERR2337147_2.fastq.gz}
OUT=results/q1
mkdir -p "$OUT"

# (a) index the reference, align the pair, sort and index the alignments
bwa index "$REF"
samtools faidx "$REF"                                  # lets igv load the genome
bwa mem -t 4 "$REF" "$R1" "$R2" | samtools view -bS - > "$OUT/ERR2337147.bam"
samtools sort "$OUT/ERR2337147.bam" -o "$OUT/ERR2337147.sorted.bam"
samtools index "$OUT/ERR2337147.sorted.bam"
samtools flagstat "$OUT/ERR2337147.sorted.bam"

# help pick a well-covered region for the igv screenshot: mean depth per 1 kb window
samtools depth -a "$OUT/ERR2337147.sorted.bam" \
    | awk '{w=int(($2-1)/1000); s[$1"\t"w*1000+1"-"w*1000+1000]+=$3}
           END{for(k in s) print k"\t"s[k]/1000}' \
    | sort -k3,3gr | head -5
# then:  igv -g "$REF" "$OUT/ERR2337147.sorted.bam"   and go to one of those regions

# (b) the 10th record of the sorted alignment file and its FLAG
samtools view "$OUT/ERR2337147.sorted.bam" | sed -n 10p
python3 src/q1_decode_flag.py "$(samtools view "$OUT/ERR2337147.sorted.bam" | sed -n 10p | cut -f2)"

# (c) reads with MAPQ < 10: CIGAR summary, then examples with MAPQ in [1, 5]
echo "MAPQ<10 (mapped) records: $(samtools view -F 4 "$OUT/ERR2337147.sorted.bam" | awk '$5<10' | wc -l)"
echo "most common CIGARs among them:"
samtools view -F 4 "$OUT/ERR2337147.sorted.bam" | awk '$5<10 {print $6}' | sort | uniq -c | sort -k1,1nr | head -15
echo "fraction with soft/hard clipping (S/H):"
samtools view -F 4 "$OUT/ERR2337147.sorted.bam" | awk '$5<10 {n++; if ($6 ~ /[SH]/) c++} END {print c "/" n}'
echo "examples with 1 <= MAPQ <= 5 (QNAME FLAG RNAME POS MAPQ CIGAR + XA/SA tags):"
samtools view -F 4 "$OUT/ERR2337147.sorted.bam" \
    | awk '$5>=1 && $5<=5 {tags=""; for(i=12;i<=NF;i++) if($i ~ /^(XA|SA|AS|XS|NM):/) tags=tags" "$i;
                            print $1, $2, $3, $4, $5, $6 tags}' | head -10
