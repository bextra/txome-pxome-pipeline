#!/bin/bash
# After Sixpack has run, it creates seperate files for each single sequence
# These need to be combined into one single fasta file of all possible ORFs

echo "starting fasta file concatenation"

# files must be listed for the full folder
# subsetting for fasta only will overload "ls"
ls ./output_Sixpack_gg/ > list_sixpack_outputs.txt
grep "fasta" list_sixpack_outputs.txt > list_sixpack_outputs2.txt
echo "files listed"

echo "" > genome_guided_merged.fasta # prevent adding a non-empty file

while read file; do
	cat ./output_Sixpack_gg/$file >> genome_guided_merged.fasta
done < list_sixpack_outputs2.txt

# clean up
rm list_sixpack_outputs2.txt
rm list_sixpack_outputs.txt

echo "completed concatentation"
