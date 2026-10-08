#!/usr/bin/env bash
# COMP3353 Lab Session 2 — full pipeline: download, align with BWA, then
# answer Q1–Q6 (plus the numeric parts of Q8/Q9) with samtools.
#
# Usage:  bash run_lab02.sh [workdir]        (default workdir: ~/labsession2)
# Needs:  bwa, samtools (conda activate comp3353), curl, unzip, python3
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-$HOME/labsession2}"
mkdir -p "$WORK" && cd "$WORK"

REF=genomic.fna
CHR=NC_045512.2
BAM=SRR12351819.bam
SORTED=SRR12351819.sorted.bam

# ---------- 2.1 reference genome ----------
if [[ ! -s $REF ]]; then
    case "$(uname -s)" in
        Darwin) OS=mac ;;
        *)      OS=linux-amd64 ;;
    esac
    [[ -x datasets ]] || { curl -sSL -o datasets \
        "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/$OS/datasets"; chmod +x datasets; }
    ./datasets download virus genome taxon 2697049 --refseq --filename SARS2-reference.zip
    unzip -o -q SARS2-reference.zip
    cp ./ncbi_dataset/data/genomic.fna ./
fi
head -n 10 "$REF"

# ---------- 2.2 – 2.5 reads, alignment, sort, index ----------
if [[ ! -s $BAM ]]; then
    for i in 1 2; do
        f=SRR12351819_$i.fastq.gz
        [[ -s $f ]] || curl -L \
            "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/019/SRR12351819/$f" > "$f"
    done
    bwa index "$REF"
    bwa mem -t 4 "$REF" SRR12351819_1.fastq.gz SRR12351819_2.fastq.gz \
        | samtools view -bS - > "$BAM"
fi
echo "== 2.3 first two alignments"
samtools view "$BAM" | head -n 2 || true   # head closes the pipe early (SIGPIPE)
samtools sort "$BAM" -o "$SORTED"
samtools index "$SORTED"

frac() { python3 -c "print(f'{$1/$2:.4f} ({100*$1/$2:.2f}%)')"; }

# ---------- Q1 ----------
echo; echo "== Q1: number of alignments in $BAM"
TOTAL=$(samtools view -c "$BAM")
echo "samtools view -c $BAM -> $TOTAL  (all records)"
echo "  of which mapped (-F 4):                 $(samtools view -c -F 4 "$BAM")"
echo "  primary mapped (-F 0x904):              $(samtools view -c -F 0x904 "$BAM")"
echo "  secondary/supplementary (-f 0x100/0x800): $(samtools view -c -f 0x100 "$BAM") / $(samtools view -c -f 0x800 "$BAM")"

# ---------- Q2 ----------
echo; echo "== Q2: alignments overlapping $CHR:11075-11102"
echo "samtools view -c $SORTED $CHR:11075-11102 -> $(samtools view -c "$SORTED" "$CHR:11075-11102")"

# ---------- Q3 ----------
echo; echo "== Q3: positive TLEN statistics"
samtools view "$BAM" | awk '$9 > 0 { print $9 }' | sort -n | uniq -c | awk '
    NR == 1 { min = $2 }
    { max = $2; n += $1; if ($1 > best) { best = $1; mode = $2 } }
    END { printf "entries with TLEN>0: %d\nmin TLEN: %d\nmax TLEN: %d\nmost frequent TLEN: %d (seen %d times)\n", n, min, max, mode, best }'

# ---------- Q4 ----------
echo; echo "== Q4: properly paired"
samtools flagstat "$BAM"
PAIRED=$(samtools view -c -f 1 "$BAM")
PROPER=$(samtools view -c -f 2 "$BAM")
echo "paired (-f 1): $PAIRED   properly paired (-f 2): $PROPER"
echo "fraction properly paired = $PROPER / $PAIRED = $(frac "$PROPER" "$PAIRED")"

# ---------- Q5 ----------
echo; echo "== Q5: strand (FLAG 0x10), mapped alignments only"
POS=$(samtools view -c -F 0x14 "$BAM")       # not unmapped, not reverse
NEG=$(samtools view -c -f 0x10 -F 0x4 "$BAM") # reverse, not unmapped
echo "positive strand (-F 20):       $POS"
echo "negative strand (-f 16 -F 4):  $NEG"
echo "(without excluding unmapped: -F 16 -> $(samtools view -c -F 16 "$BAM"), -f 16 -> $(samtools view -c -f 16 "$BAM"))"

# ---------- Q6 ----------
echo; echo "== Q6: MAPQ >= 20"
Q20=$(samtools view -c -q 20 "$BAM")
echo "samtools view -c -q 20 $BAM -> $Q20"
echo "fraction of all alignments = $Q20 / $TOTAL = $(frac "$Q20" "$TOTAL")"

# ---------- Q8 / Q9 (numeric support for the IGV questions) ----------
echo; echo "== Q8: ORF8 = $CHR:27894-28259, last base 28259"
echo "reads overlapping 28259 (IGV default filters: no dups/QC-fail): $(samtools view -c -F 0x604 "$SORTED" "$CHR:28259-28259")"
samtools depth -a -Q 0 -q 0 -r "$CHR:28259-28259" "$SORTED" | awk '{print "samtools depth at "$2": "$3}'

echo; echo "== Q9: variants inside ORF3a ($CHR:25393-26220)"
samtools mpileup -f "$REF" -r "$CHR:25393-26220" -Q 0 -B -d 0 "$SORTED" 2>/dev/null \
    | python3 "$SRC_DIR/pileup_variants.py"
