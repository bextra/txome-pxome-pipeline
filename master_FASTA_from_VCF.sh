#!/bin/bash

# 1. Prepare BAM files for generating VCF
# Remove duplicates from RNA-Seq reads (alleviates pseudo-replication due to PCR amplification)
# Dependencies: samtools (v0.1.7a, other versions may be acceptable if working on the new cluster)
#./bin/dedup_txpx_paired.sh


# 2. Make VCF file
# - reference fasta required
# - BAM files required
# Dependencies: samtools (v1.2 or newer), bcftools (v1.2-22)
#qsub -S /bin/bash -pe threaded 2 ./bin/vcf_qsub.sh
# compute resources may change as desired


# 3. Prepare BAM files for use with customProDB
# Optional: if your BAM files are not indexed, change INDEXBAMS to = TRUE
INDEXBAMS=FALSE

if [ $INDEXBAMS == TRUE ]; then
	module load samtools/1.2
	for i in /share/milklab/proteomics/BAM_files/human*.bam
	   do
		echo "indexing $i"
		samtools index $i
	   done
	exit 1;
fi


# 4. Run customProDB
Rscript ./bin/fasta_from_customProDB.R


