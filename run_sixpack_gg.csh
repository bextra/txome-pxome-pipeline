#!/bin/csh

#foreach file in the input directory
# run Sixpack on it to create the 6 ORFs
# Note that orfminsize = minimum NUCLEOTIDE size

rm -r ./output_Sixpack_gg
mkdir ./output_Sixpack_gg
echo "starting Sixpack"
foreach fname (./input_Sixpack_gg/*.fa)
  set froot = $fname:t:r
  sixpack -sequence $fname -orfminsize 60 -nofirstorf -nolastorf -outfile ./output_Sixpack_gg/{$froot}.sixpack -outseq ./output_Sixpack_gg/{$froot}.fasta
end
echo "finished Sixpack"

