# fasta_from_customProDB.R
# K. Beck

# Objective:
# This script generates a fasta file from input RNA-Seq reads (BAM format) and a VCF to assist with SNP detection and sequence variants.

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

# To explore available reference datasets in Biomart
# listMarts() # lists all types of data

annotation_path_mm = "/share/milklab/proteomics/VariantCalling/output_customProDB/"
annotation_path_hs = "/share/milklab/proteomics/VariantCalling/output_customProDB/"
# alternatively, could use tmepdir() to generate a unique path

# # # # # # # # 
#
# Get annotation data
#
# # # # # # # # 
# From Enseml
# Versions corresponding to our BAM files
  # Macaque: EnsRel67 (MMUL1p0_EnsRel67), May 2012
retrieveReference_mm = FALSE
retrieveReference_hs = FALSE

if(retrieveReference_mm == TRUE){
	cat("Retrieving annotation data\n")
	ensembl = useMart(biomart="ensembl") # creates a Mart object
	# listDatasets(ensembl) # view available datasets (optional)
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmulatta_gene_ensembl", host="may2012.archive.ensembl.org")
	PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path_mm, 
                          splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)
}

if(retrieveReference_hs == TRUE){
	ensembl = useMart(biomart="ensembl") # creates a Mart object
	# listDatasets(ensembl) # view available datasets (optional) 
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="oct2014.archive.ensembl.org")
	PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path_hs, 
                          splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)

# Retrieve annotation data from UCSC for human
  # Human: GRCh37
# This version matches what DGL originally used to make the BAM files
# Note: Annotation data must first be downloaded from online interface per vignette
# pepfasta = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/hg19_GRCh37_proseq.fasta"
# CDSfasta = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/hg19_GRCh37_codingseq.fasta"
#PrepareAnnotationRefseq(genome='hg19', CDSfasta, pepfasta, annotation_path=annotation_path_hs, 
#                         splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)
# option to filter for specific transcript ids by storing them as a vector and then setting transcrip_ids= to the that vector above


# # # # # # # # 
#
# Generate fasta from BAM files 
# filtering on RPKM
#
# # # # # # # # 

run_customProDB = function(annotation_path="./", outfile="custom.fasta", singleSample = TRUE, path_to_sample="./", correct_chr_name=FALSE){
 
  # Load annotation data
  cat("Loading annotation data...\n")
  load(paste(annotation_path, "exon_anno.RData", sep = ""))
  load(paste(annotation_path, "ids.RData",       sep = ""))
  load(paste(annotation_path, "proseq.RData",    sep = ""))
  
  # re-map chromosome names to match BAM file (if needed, typcial for Ensembl file, not for UCSC)
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
    bamFiles = paste(path_to_sample, list.files(path_to_sample, pattern = "*.bam$"), sep="") # Load all bam files within one directory
    RPKMs = sapply(bamFiles, function(x) calculateRPKM(x, exon, proteincodingonly = TRUE, ids))
    pro = OutputsharedPro(RPKMs, cutoff=1, share_sample=2, proteinseq, outf1, ids) # for multiple samples
  
  }
}

# Macaque
run_customProDB(annotation_path= annotation_path_mm, outfile="monkey_2_1_index18-customProDB.fasta", singleSample=TRUE, correct_chr_name=TRUE,
                path_to_sample="~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/monkey_2_1_index18_accepted_hits.bam")

# Human
# For multiple samples this will run on all bam files within one folder
run_customProDB(annotation_path= annotation_path_hs, outfile="multiple_sample_Custo  mProDB.fasta", singleSample=FALSE, correct_chr_name=FALSE, path_to_sample="~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/Human/")

# # # # # # # # 
#
# Generate fasta from BAM files 
# with variant calling
#
# # # # # # # # 
vcffile = "/share/milklab/proteomics/VariantCalling/macaque_var.flt.vcf"
vcf = InputVcf(vcffile)

# vcf[[1]][1:3] # pull an example range of variants
if (table(values(vcf[[1]])[['INDEL']]) < 5) {
	cat("Warning: less than 5 INDELs check VCF for quality\n")
}



# ------------------------ #

# Research and Notes
# https://www.biostars.org/p/63429/



# ------------------------ #
# Check if our BAM file matches the exon data in database
tmp = read.table("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/tmp.txt", header = FALSE)
table(tmp)
table(exon$chromosome_name)

# tmp
# chr1   chr10   chr11   chr12   chr13   chr14   chr15   chr16   chr17   chr18   chr19    chr2   chr20   chr21   chr22 
# 1147846  424346  650613 1888193  242957  307738  347684  332824  747426  132260  538888  779291  227360  116425  216985 
# chr3    chr4    chr5    chr6    chr7    chr8    chr9    chrM    chrX    chrY 
# 895167 2798164  573319  655589  594033  430142  688995  237872  323777   71465 




chrY_tmp = read.table("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/tmp2.txt", sep= "\t", fill = TRUE, header = FALSE)

# regions of homoology between X and Y chromosomes
# http://en.wikipedia.org/wiki/Pseudoautosomal_region

options(scipen = 10)
par(mfrow = c(1,1))
hist(chrY_tmp$V4, breaks = 25,
     xlim = c(0, max(chrY_tmp$V4, na.rm = TRUE)))

plot(density(chrY_tmp$V4, na.rm = TRUE))











