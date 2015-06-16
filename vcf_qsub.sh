#!/bin/bash

# load samtools
module load samtools/1.2

# add local bcftools to path
PATH=$PATH:/share/milklab/proteomics/Tools/bcftools/Current
export PATH

#make writing file names easier
FILEROOT=/share/milklab/proteomics/VariantCalling

# generate VCF file
samtools mpileup -DVuf $FILEROOT/ReferenceData/rheMac2.fa $FILEROOT/*.bam | bcftools view > $FILEROOT/monkey_pxtx_paired.vcf

