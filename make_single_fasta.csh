#!/bin/csh
# After Sixpack has run, it creates seperate files for each single sequence
# These need to be combined into one single fasta file of all possible ORFs
echo "starting cat"
#cat ./output_Sixpack/*.fasta > denovo_Sixpack_ORFs.fasta 
cat ./output_Sixpack/*.fasta > denovo_Sixpack_ORFs_m2_1_i18.fasta 
echo "completed cat"
