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
session = sessionInfo()
customProDBVersion = unlist(str_split(session$otherPkgs$customProDB$Version, "\\."))
if(!(customProDBVersion[1] >= 1 & customProDBVersion[2] >= 8 & customProDBVersion[3] >= 2)) {
  message("Please upgrade customProDB to version 1.8.2 or greater\n")
}

# To explore available reference datasets in Biomart
# listMarts() # lists all types of data

setwd("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/VariantCalling/")
# KB local only

annotation_path_mm = "/share/milklab/proteomics/run_customProDB/" 

annotation_path_hs = "./"
# Cluster
# annotation_path_hs = "/share/milklab/proteomics/run_customProDB/"
# alternatively, could use tmepdir() to generate a unique path


# # # # # # # # 
#
# Get annotation data
#
# # # # # # # # 
# From Ensembl
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
                          splice_matrix = FALSE, dbsnp=NULL, COSMIC=TRUE)
}



# Ensembl human data as option for other users
if(retrieveReference_hs == TRUE){
	ensembl = useMart(biomart="ensembl") # creates a Mart object
	# listDatasets(ensembl) # view available datasets (optional) 
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="oct2014.archive.ensembl.org")
	PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path_hs, 
                          splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)
}
# Retrieve annotation data from UCSC for human when necessary
  # Human: GRCh37
# This version matches what DGL originally used to make the BAM files
# Note: Annotation data must first be downloaded from online interface per vignette
pepfasta = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/hg19_GRCh37_proseq.fasta"
CDSfasta = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/hg19_GRCh37_codingseq.fasta"
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
  
  # re-map chromosome names to match BAM file 
  # (optional: typically required for Ensembl, not for UCSC)
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
    # for multiple samples
    OutputsharedPro(RPKMs, cutoff=1, share_sample=2, proteinseq, outf1, ids) 
  
  }
}
cat("Completed loading customProDB function...\n")




# Macaque on Cluster - single test sample
# run_customProDB(annotation_path= annotation_path_mm, 
#   outfile="monkey_2_1_index18-customProDB.fasta", singleSample=TRUE, 
#   correct_chr_name=TRUE, 
#   path_to_sample="/share/milklab/proteomics/BAM_files/monkey_2_1_index18.bam")

# Macaque for Cluster - multiple bio reps
# run_customProDB(annotation_path=annotation_path_mm, 
#   outfile="monkey_all_customProDB.fasta", singleSample=FALSE, 
#   correct_chr_name=TRUE, path_to_sample="/share/milklab/proteomics/BAM_files/")

# Human
# TODO update for cluster
#run_customProDB(annotation_path= annotation_path_hs, 
#   outfile="multiple_sample_CustomProDB.fasta", singleSample=FALSE, 
#   correct_chr_name=FALSE, 
#   path_to_sample="~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/Human/")

run_customProDB(annotation_path= annotation_path_hs, 
  outfile="hs_multiple_sample_CustomProDB.fasta", singleSample=FALSE, 
  correct_chr_name=FALSE, 
  path_to_sample="~/Work/1_Milk/RNA-Seq_Guided_Proteomics/BAMs/Human/")

# # # # # # # # 
#
# Generate fasta from BAM files 
# with variant calling
#
# # # # # # # # 

# There were significant issues trying to get the variant annotation portion
# of customProDB to work. These are described below.
# Code that is present represents the portions that work unless otherwise
# noted with an error message

# dbSNP data is only retrievable from customprodb for human and mouse organisms
# https://support.bioconductor.org/p/69592/
# dbsnp would not accept any of the parameters suggested in the package vignette
# or in the help page. This has been posted to a bioconductor help page:
# https://support.bioconductor.org/p/69670/

# per the developer's reply, you can construct it with snp138

# Load in VCF File
# monkey vcf on cluster
# vcffile = "/share/milklab/proteomics/VariantCalling/updated_monkey_pxtx_paired.vcf"
vcffile = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/VariantCalling/human_pxtx_paired.vcf"
vcf = InputVcf(vcffile)

# vcf[[1]][1:3] # pull an example range of variants
if (table(values(vcf[[1]])[['INDEL']])[2] < 5) {
	cat("Warning: less than 5 INDELs check VCF for quality\n")
}

# subset the object for indels
index <- which(values(vcf[[1]])[['INDEL']]==TRUE) # get rows that are indels
indelvcf <- vcf[[1]][index] # subset the object

# subset the object for SNVs
index <- which(values(vcf[[1]])[['INDEL']]==FALSE)
SNVvcf <- vcf[[1]][index]

load(paste(annotation_path_hs, "ids.RData", sep=""))

txdb <- loadDb(paste(annotation_path_hs, "txdb.sqlite", sep=""))
# this chr is named 1, 2, 3, etc.

SNVloc <- Varlocation(SNVvcf,txdb,ids)
table(SNVloc$location)
# TODO fix this bug... does the chr name in the sqlite db need to be prefixed with chr?
# QQ all locations are returned as unknown... why?
# Warning: In .Seqinfo.mergexy(x, y) :
#  The 2 combined objects have no sequence levels in common.

indelloc <- Varlocation(indelvcf,txdb,ids)
table(indelloc[,'location'])
# TODO is this caused by the bug above?
# still listed as all unknown

load("exon_anno.RData")
load("dbsnpinCoding.RData")
load("cosmic.RData")

postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding, COSMIC=cosmic)
postable_indel <- Positionincoding(indelvcf, exon)
# TODO another thing to troubleshoot
# Error in .Call2("solve_user_SEW0", start, end, width, PACKAGE = "IRanges") : 
#  solving row 1: range cannot be determined from the supplied arguments (too many NAs)

# TODO the dbsnp and COSMIC info is only retrieval for human and mouse
# TODO could re-make a dbsnp database
  # need to count variant allele's and make them separated by a comma
  # need to convert chromosome positions to a range
  # need to change strandedness from + to - etc
# TODO what about cosmic? what is that?


# continuation with human data
load("procodingseq.RData")
txlist <- unique(postable_snv[, 'txid'])
codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist,]
mtab <- aaVariation(postable_snv, codingseq)
table(mtab$vartype) # 1259 non-synonymous changes
# Warning messages:
#   1: In .Method(..., deparse.level = deparse.level) :
#   number of columns of result is not a multiple of vector length (arg 3)

outfile <- "snv_human.fasta"
load("proseq.RData")

OutputVarproseq(mtab, proteinseq, outfile, ids)
# results in 704 proteins non-synonymous SNVs into FASTA file

txlist_indel <- unique(postable_indel[, 'txid'])
codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] %in% txlist_indel, ]
# there are 30 out of the 35 indels (this number is down from previous ~1500, why?)
outfile <- "indel_human.fasta"

# TODO work from here

Outputaberrant(postable_indel, 
               coding=codingseq_indel, # there are 5 less than postable_indel
               proteinseq=proteinseq, 
               outfile=outfile, ids=ids)
# Error in translate(DNAStringSet(total[, "coding"])) : 
#   error in evaluating the argument 'x' in selecting a method for function 'translate': Error in width(x) : NAs in 'x' are not supported




# ------------------------ #

# Research and Notes
# https://www.biostars.org/p/63429/



# ------------------------ #
# Check if our BAM file matches the exon data in database
#tmp = read.table("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/tmp.txt", header = FALSE)
#table(tmp)
#table(exon$chromosome_name)

# tmp
# chr1   chr10   chr11   chr12   chr13   chr14   chr15   chr16   chr17   chr18   chr19    chr2   chr20   chr21   chr22 
# 1147846  424346  650613 1888193  242957  307738  347684  332824  747426  132260  538888  779291  227360  116425  216985 
# chr3    chr4    chr5    chr6    chr7    chr8    chr9    chrM    chrX    chrY 
# 895167 2798164  573319  655589  594033  430142  688995  237872  323777   71465 




#chrY_tmp = read.table("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/tmp2.txt", sep= "\t", fill = TRUE, header = FALSE)

# regions of homoology between X and Y chromosomes
# http://en.wikipedia.org/wiki/Pseudoautosomal_region

#options(scipen = 10)
#par(mfrow = c(1,1))
#hist(chrY_tmp$V4, breaks = 25,
#     xlim = c(0, max(chrY_tmp$V4, na.rm = TRUE)))

#plot(density(chrY_tmp$V4, na.rm = TRUE))











# *Note:* tested the easyRun and easyRunMul functions, but they both resulted in
# an non-descript error