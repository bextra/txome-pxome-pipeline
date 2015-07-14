# fasta_from_customProDB.R
# K. Beck

# Objective:
# This script generates a fasta file from input RNA-Seq reads (BAM format) and 
# a VCF to assist with SNP detection and sequence variants.

# # # # # # # # 
#
# SETUP
#
# # # # # # # # 

# if package is not installed uncomment the lines below
# source("http://bioconductor.org/biocLite.R")
# biocLite("customProDB")

# Load required packages
require("customProDB")
require("stringr")

# quality control: make sure customProDB is up to date
customProDBVersion = unlist(str_split(packageDescription("customProDB")$Version, "\\."))
if(!(customProDBVersion[1] >= 1 & customProDBVersion[2] >= 8 & customProDBVersion[3] >= 2)) {
  stop("Please upgrade customProDB to version 1.8.2 or greater\n")
}

# To explore available reference datasets in Biomart
# listMarts() # lists all types of data

setwd("/share/milklab/proteomics/run_customProDB/")

# Path to input of annotation data for UCSC annotated samples
# Will also be path to subsequent retrieved reference data
annotation_path_mm = "/share/milklab/proteomics/run_customProDB/Macaque_Ref"
annotation_path_hs = "/share/milklab/proteomics/run_customProDB/Human_Ref"
# alternatively, could use tmepdir() to generate a unique path


# # # # # # # # 
#
# Get annotation data
#
# # # # # # # # 
# From Ensembl
# Versions corresponding to our BAM files
  # Macaque: EnsRel67 (MMUL1p0_EnsRel67), May 2012
retrieve_macaque_ensembl = FALSE
retrieve_human_ensembl   = FALSE

# Retrieve reference data from Ensembl
# This version matches what DGL originally used to make the BAM files
if(retrieve_macaque_ensembl == TRUE){
	cat("Retrieving annotation data\n")
	ensembl = useMart(biomart="ensembl") # creates a Mart object
	# listDatasets(ensembl) # view available datasets (optional)
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmulatta_gene_ensembl", host="may2012.archive.ensembl.org")
	PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path_mm, 
                          splice_matrix = FALSE, dbsnp=NULL, COSMIC=TRUE)
}

# Retrieve Ensembl human data (option for other users, not needed for milk project)
if(retrieve_human_ensembl == TRUE){
	ensembl = useMart(biomart="ensembl") # creates a Mart object
	# listDatasets(ensembl) # view available datasets (optional) 
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="oct2014.archive.ensembl.org")
	PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path_hs, 
                          splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)
}


# Retrieve annotation data from UCSC for human when necessary
  # Human: GRCh37
# This version matches what DGL originally used to make the BAM files
# NOTE: Annotation data must first be downloaded from online interface per 
# instructions in the package vignette
pepfasta = "/share/milklab/proteomics/run_customProDB/Human_Ref/hg19_GRCh37_proseq.fasta"
CDSfasta = "/share/milklab/proteomics/run_customProDB/Human_Ref/hg19_GRCh37_codingseq.fasta"
PrepareAnnotationRefseq(genome='hg19', CDSfasta, pepfasta, 
                        annotation_path=annotation_path_hs, 
                        splice_matrix = FALSE, dbsnp = "snp138", COSMIC = TRUE
                        )
# option to filter for specific transcript ids by storing them as a vector and 
# then setting transcrip_ids= to the that vector above


# # # # # # # # 
#
# Generate fasta from BAM files 
# filtering on RPKM
#
# # # # # # # # 
cat("Loading customProDB function...\n")
run_customProDB = function(annotation_path="./", outfile="custom.fasta", singleSample = TRUE, path_to_sample="./", correct_chr_name=FALSE){
  # Load annotation data
  cat("Loading annotation data...\n")
  load(paste(annotation_path, "exon_anno.RData", sep = ""))
  load(paste(annotation_path, "ids.RData",       sep = ""))
  load(paste(annotation_path, "proseq.RData",    sep = ""))
  
  # Re-map chromosome names to match BAM file 
  # *optional:* typically required for Ensembl, not for UCSC)
  if (correct_chr_name == TRUE) {
    cat("Correcting chromosome name...\n")
    exon$chromosome_name = paste("chr", exon$chromosome_name, sep = "") 
  }
  
  outf1 = paste(annotation_path, outfile, sep='')
  
  if (singleSample == TRUE) {
    cat("Filtering single sample on RPKM...\n")
    bamFile = path_to_sample
    RPKM = calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
    Outputproseq(RPKM, 1, proteinseq, outf1, ids) # for single sample
    
  } else if (singleSample == FALSE) {
    cat("Filtering multiple samples on RPKM...\n")
    # Load all bam files within one directory
    bamFiles = paste(path_to_sample, 
                     list.files(path_to_sample, pattern = "*.bam$"), sep="") 
    RPKMs = sapply(bamFiles, function(x) 
      calculateRPKM(x, exon, proteincodingonly = TRUE, ids)
      )
    # For multiple samples
    OutputsharedPro(RPKMs, cutoff=1, share_sample=2, proteinseq, outf1, ids) 
  
  }
}
cat("Completed loading customProDB function...\n")



# Example for Macaque single test sample
# run_customProDB(annotation_path= annotation_path_mm, 
#   outfile="monkey_2_1_index18-customProDB.fasta", singleSample=TRUE, 
#   correct_chr_name=TRUE, 
#   path_to_sample="/share/milklab/proteomics/BAM_files/monkey_2_1_index18.bam")

# Macaque for Cluster - multiple bio reps
# All bams are in a single folder
# run_customProDB(annotation_path=annotation_path_mm, 
#   outfile="monkey_all_customProDB.fasta", singleSample=FALSE, 
#   correct_chr_name=TRUE, path_to_sample="/share/milklab/proteomics/BAM_files/")

# Human for cluster
# All bams are in a single folder
run_customProDB(annotation_path= annotation_path_hs, 
  outfile="hs_multiple_sample_CustomProDB.fasta", singleSample=FALSE, 
  correct_chr_name=FALSE, 
  path_to_sample="/share/milklab/proteomics/BAM_files/Human/")

# # # # # # # # 
#
# Generate variant called fasta
#
# # # # # # # # 

# Note: customProDB was written for variant annotation for only human and mouse
# Developer response here: https://support.bioconductor.org/p/69592/
# dbsnp would not accept any of the parameters suggested in the package vignette
# or in the help page (6/13/15).
# Per the developer's reply, you can construct dbSNP using snp138
# https://support.bioconductor.org/p/69670/



# Load in VCF File
vcffile = "/share/milklab/proteomics/VariantCalling/human_pxtx_paired_all_tech_reps.vcf"
vcf = InputVcf(vcffile)

# quality control vcf
if (table(values(vcf[[1]])[['INDEL']])[2] < 5) {
	stop("Less than 5 INDELs check VCF for quality\n")
}

# subset the object for indels
index    = which(values(vcf[[1]])[['INDEL']]==TRUE) # get rows that are indels
indelvcf = vcf[[1]][index] # subset the object

# subset the object for SNVs
index  = which(values(vcf[[1]])[['INDEL']]==FALSE)
SNVvcf = vcf[[1]][index]

load(paste(annotation_path_hs, "ids.RData", sep=""))
txdb = loadDb(paste(annotation_path_hs, "txdb.sqlite", sep=""))


SNVloc = Varlocation(SNVvcf,txdb,ids)
table(SNVloc$location)

indelloc = Varlocation(indelvcf,txdb,ids)
table(indelloc[,'location'])


load(paste(annotation_path_hs, "exon_anno.RData", sep=""))
load(paste(annotation_path_hs, "dbsnpinCoding.RData", sep=""))
load(paste(annotation_path_hs, "cosmic.RData", sep=""))

postable_snv   = Positionincoding(SNVvcf, exon, dbsnpinCoding, COSMIC=cosmic)
postable_indel = Positionincoding(indelvcf, exon)

load(paste(annotation_path_hs, "procodingseq.RData", sep=""))
txlist = unique(postable_snv[, 'txid'])
codingseq = procodingseq[procodingseq[, 'tx_id'] %in% txlist,]
mtab = aaVariation(postable_snv, codingseq)
table(mtab$vartype) # 1259 non-synonymous changes


outfile = "snv.fasta"
load(paste(annotation_path_hs, "proseq.RData", sep=""))
OutputVarproseq(mtab, proteinseq, outfile, ids)


txlist_indel = unique(postable_indel[, 'txid'])
codingseq_indel = procodingseq[procodingseq[, 'tx_id'] %in% txlist_indel, ]
# there are 30 out of the 35 indels (this number is down from previous ~1500, why?)
outfile = "indel_human.fasta"


Outputaberrant(postable_indel, 
               coding=codingseq_indel, # there are 5 less than postable_indel
               proteinseq=proteinseq, 
               outfile=outfile, ids=ids)



# *Note:* tested the easyRun and easyRunMul functions for human, but they both 
# resulted in an non-descript error