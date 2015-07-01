#!/bin/bash

# Genome-guided Trinity
# the reads were first aligned to rheMac2
# Trinity will partition according to locus, followed by de novo transcriptome assmebly at each locus

#Read alignments must be supplied as coordinate-sorted bam file.
cd /Chanlab/Scratch/kristen/Macaque_transcriptome/BAMs

#SAMPLE_LIST="$(ls)"

SAMPLE_LIST="monkey_2_1_index18.bam"

for FILE in $SAMPLE_LIST
do
	samtools sort $FILE sorted_$FILE
done

for FILE in $SAMPLE_LIST
do
	/Chanlab/Packages/Trinity/trinityrnaseq-2.0.6/Trinity --genome_guided_bam sorted_$FILE --genome_guided_max_intron 10000 -seqType fq --max_memory 20G --CPU 20 --output ./trinity_gg_$FILE --full_cleanup
done
