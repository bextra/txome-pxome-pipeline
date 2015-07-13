#!/bin/bash

# 1. Prepare BAM files for generating VCF
# Remove duplicates from RNA-Seq reads (alleviates pseudo-replication due to PCR amplification)
# Dependencies: samtools (v0.1.7a, other versions may be acceptable if working on the new cluster)
./bin/dedup_txpx_paired.sh

# 1. Make VCF file
# - reference fasta required
# - BAM files required
# Dependencies: samtools (v1.2 or newer), bcftools (v1.2-22)
./bin/vcf_qsub.sh

# 2. Prepare BAM files for use with customProDB
# Optional: if your BAM files are not indexed, change INDEXBAMS to = TRUE
INDEXBAMS=FALSE

if $INDEXBAMS == TRUE
	module load samtools/1.2
	samtools index /share/milklab/BAM_files/*


# generate index file using
module load samtools/0.1.7a
samtools index [file.bam]

# use samtools and bcftools to generate VCF file
use local bcftools in /share/milklab/proteomics/Tools/bcftools/bcftools-1.2-22/
generate with ...
/share/milklab/proteomics/Code/txome-pxome-pipeline/vcf_qsub.sh

