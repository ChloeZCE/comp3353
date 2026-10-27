#!/usr/bin/env bash
# Runs Q2-Q6 and stores every output under results/.  Usage (from assignment2/): bash src/run_all.sh
set -uo pipefail
cd "$(dirname "$0")"
R=../results; D=../data
run() { echo "\$ $*"; "$@"; echo; }
{ run python3 q2_hamming.py $D/Q2/input.fa; }                                   > $R/q2_output.txt 2>&1
{ run python3 q3_edit_distance.py $D/Q2/input.fa
  run python3 q3_edit_distance.py $D/Q3/input.fa; }                             > $R/q3_output.txt 2>&1
{ run python3 bwt.py encode $D/Q4/input1.fa
  echo "\$ python3 bwt.py encode $D/Q4/input2.fa | cut -c1-50"
  python3 bwt.py encode $D/Q4/input2.fa | cut -c1-50; echo; }                    > $R/q4a_output.txt 2>&1
python3 bwt.py encode $D/Q4/input1.fa > $R/input1.bwt.txt
{ run python3 bwt.py decode $D/Q4/input3.txt
  run python3 bwt.py decode $R/input1.bwt.txt
  echo "original input1.fa sequence:"
  python3 -c "from fasta import read_fasta; print(read_fasta('$D/Q4/input1.fa')[0][1])"
  echo "identical: $(python3 -c "
from fasta import read_fasta; from bwt import bwt_decode
print(bwt_decode(open('$R/input1.bwt.txt').read().strip()) == read_fasta('$D/Q4/input1.fa')[0][1])")"
}                                                                               > $R/q4b_output.txt 2>&1
{ run python3 q5_fm_search.py $D/Q4/input1.fa ACGT --show-occ 1000
  run python3 q5_fm_search.py $D/Q4/input2.fa ACGT; }                           > $R/q5_output.txt 2>&1
{ run python3 q6_variant_calling.py; }                                          > $R/q6_output.txt 2>&1
tail -n +1 $R/q[2-6]*_output.txt | cut -c1-300
