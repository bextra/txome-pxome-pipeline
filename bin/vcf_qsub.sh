#!/bin/bash

# load samtools
module load samtools/1.2

# add local bcftools to path
PATH=$PATH:/share/milklab/proteomics/Tools/bcftools/Current
export PATH

#make writing file names easier
FILEROOT=/share/milklab/proteomics/VariantCalling
BAMPATH=/share/milklab/proteomics/BAM_files

# generate VCF file
# for macaque:
# samtools mpileup -DVuf $FILEROOT/ReferenceData/rheMac2.fa $FILEROOT/*.bam | bcftools call -vm > $FILEROOT/updated_monkey_pxtx_paired.vcf

# for human:
samtools mpileup -DVuf $FILEROOT/ReferenceData/human_ensembl.GRCh37.fa $BAMPATH/human*.bam | bcftools call -vm > $FILEROOT/human_pxtx_paired_all_tech_reps.vcf

