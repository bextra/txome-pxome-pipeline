#!/bin/csh
# Sixpack translates a sequence to 6 ORFS, but only takes one sequence in one file at a
# time. Trinity outputs all of the sequences as one multiple fasta file.

#mkdir input_Sixpack
rm input_Sixpack/*
cd input_Sixpack
#~/work/perl/split_fasta.pl ../trinity_out_dir/Trinity.fasta 
~/work/perl/split_fasta.pl ../trinity_out_m2_1_i18.Trinity.fasta
 

