#!/usr/bin/env bash
# Run from assignment2/: re-creates every output file in results/
cd "$(dirname "$0")"
python3 q2_hamming.py ../data/Q2/input.fa                      > ../results/q2_output.txt
{ python3 q3_edit_distance.py ../data/Q2/input.fa
  python3 q3_edit_distance.py ../data/Q3/input.fa; }          > ../results/q3_output.txt
{ python3 bwt.py encode ../data/Q4/input1.fa
  python3 bwt.py encode ../data/Q4/input2.fa | cut -c1-50; }   > ../results/q4a_output.txt
python3 bwt.py encode ../data/Q4/input1.fa                     > ../results/input1.bwt.txt
{ python3 bwt.py decode ../data/Q4/input3.txt
  python3 bwt.py decode ../results/input1.bwt.txt; }           > ../results/q4b_output.txt
{ python3 q5_fm_search.py ../data/Q4/input1.fa ACGT
  echo
  python3 q5_fm_search.py ../data/Q4/input2.fa ACGT; }         > ../results/q5_output.txt
python3 q6_variant_calling.py                                  > ../results/q6_output.txt
rm ../results/input1.bwt.txt
