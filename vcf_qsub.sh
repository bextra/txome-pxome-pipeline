#!/bin/bash

# load samtools
module load samtools/1.2

# add local bcftools to path
PATH=$PATH:/share/milklab/proteomics/Tools/bcftools/Current
export PATH

#make writing file names easier
FILEROOT=/share/milklab/proteomics/VariantCalling

# generate VCF file
# line below is potentially wrong function call for new bcftools version
# samtools mpileup -DVuf $FILEROOT/ReferenceData/rheMac2.fa $FILEROOT/*.bam | bcftools view > $FILEROOT/monkey_pxtx_paired.vcf

# test correct bcf function
samtools mpileup -DVuf $FILEROOT/ReferenceData/rheMac2.fa $FILEROOT/*.bam | bcftools call -vm > $FILEROOT/updated_monkey_pxtx_paired.vcf

