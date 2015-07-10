#!/bin/csh
# Sixpack translates a sequence to 6 ORFS, but only takes one sequence in one file at a
# time. Trinity outputs all of the sequences as one multiple fasta file.

rm -r input_Sixpack_gg
mkdir input_Sixpack_gg
cd input_Sixpack_gg
echo "starting split_fasta"
~/work/perl/split_fasta.pl ../trinity_out_m2_1_i18/Trinity-GG.fasta 
echo "finished split fasta"

