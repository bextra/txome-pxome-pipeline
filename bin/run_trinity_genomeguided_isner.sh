#!/bin/bash

# Genome-guided Trinity
# the reads were first aligned to rheMac2
# Trinity will partition according to locus, followed by de novo transcriptome assmebly at each locus



#Read alignments must be supplied as coordinate-sorted bam file.
cd /Chanlab/Scratch/kristen/Macaque_transcriptome/BAMs

#SAMPLE_LIST="$(ls *.bam)"
SAMPLE_LIST="$(ls human*.bam)"

for FILE in $SAMPLE_LIST
do
	echo "Sorting $FILE"
	samtools sort -f $FILE sorted_$FILE
done

export JAVA_HOME=/Chanlab/Packages/Java/jdk1.7.0_80/bin/java
export PATH="/Chanlab/Packages/Java/jdk1.7.0_80/bin/:$PATH"

for FILE in $SAMPLE_LIST
do
	time /Chanlab/Packages/Trinity/trinityrnaseq-2.0.6/Trinity --genome_guided_bam /Chanlab/Scratch/kristen/Macaque_transcriptome/BAMs/sorted_$FILE --genome_guided_max_intron 10000 -seqType fq --max_memory 70G --CPU 20 --output ./trinity_gg_$FILE --full_cleanup > stdout_run_$FILE
done
