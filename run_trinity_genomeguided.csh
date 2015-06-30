#!/bin/csh

# Genome-guided Trinity
# the reads were first aligned to rheMac2
# Trinity will partition according to locus, followed by de novo transcriptome assmebly at each locus

#Read alignments must be supplied as coordinate-sorted bam file.
#samtools sort monkey_2_1_index7_accepted_hits.bam sorted_monkey_2_1_index7
samtools sort monkey_2_1_index18_accepted_hits.bam sorted_monkey_2_1_index18
# run Genome-guided Trinity
#Trinity --genome_guided_bam sorted_monkey_2_1_index7.bam --genome_guided_max_intron 10000 -seqType fq --max_memory 8G --CPU 3
Trinity --genome_guided_bam sorted_monkey_2_1_index18.bam --genome_guided_max_intron 10000 -seqType fq --max_memory 8G --CPU 3 --output ./trinity_out_m2_1_i18_gg --full_cleanup

