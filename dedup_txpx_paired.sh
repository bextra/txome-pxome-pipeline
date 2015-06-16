#!/bin/bash

# load samtools
module load samtools/0.1.7a

# array of samples that have paired RNA-Seq and proteomic data
MM_SAMPLES=(monkey_2_1_index18 monkey_3_1_index12 monkey_3_1_index16 monkey_4_1_index5 monkey_4_1_index15)

# remove duplicate reads
for i in ${MM_SAMPLES[@]}
	do
	samtools rmdup -s /share/milklab/output_TopHat/$i/accepted_hits.bam /share/milklab/proteomics/VariantCalling/$i.rmdup.bam
	samtools index /share/milklab/proteomics/VariantCalling/$i.rmdup.bam
	done
