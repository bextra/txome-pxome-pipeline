#!/bin/csh

#foreach file in the input directory
# run Sixpack on it to create the 6 ORFs
# Note that orfminsize = minimum NUCLEOTIDE size

#rm -r ./output_Sixpack 
#mkdir ./output_Sixpack
foreach fname (./input_Sixpack/*.fa)
  set froot = $fname:t:r
  sixpack -sequence $fname -orfminsize 60 -nofirstorf -nolastorf -outfile ./output_Sixpack/{$froot}.sixpack -outseq ./output_Sixpack/{$froot}.fasta
end
