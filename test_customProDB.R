# test customProDB package

# if package is not installed uncomment the lines below
# source("http://bioconductor.org/biocLite.R")
# biocLite("customProDB")

# TODO try running this without
  # otherwise you are bound by the Ensembl or RefSeq information provided
# PrepareAnnotationRefseq()

# Load required packages
library("customProDB")
library("biomaRt")

# To explore available datasets in Biomart
# listMarts() # lists all types of data

annotation_path_hs = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/"
annotation_path_mm = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Macaque/"
# alternatively, could use tmepdir()

# # # # # # # # 
#
# Get annotation data
#
# # # # # # # # 
# Versions corresponding to our BAM files
  # Human: GRCh37
  # Macaque: EnsRel67 (MMUL1p0_EnsRel67), May 2012

# Retrieve annotation data from Ensembl for macaque
ensembl = useMart(biomart="ensembl") # creates a Mart object
# listDatasets(ensembl)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmulatta_gene_ensembl", host="may2012.archive.ensembl.org")
PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path_mm, 
                          splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)

# To pull current annotation data from Ensembl, human
# ensembl.hs = useMart(biomart="ensembl") # creates a Mart object
# listDatasets(ensembl) 
# ensembl.hs = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="oct2014.archive.ensembl.org")
# PrepareAnnotationEnsembl(mart=ensembl.hs, annotation_path=annotation_path, 
#                          splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)


# Retrieve annotation data from UCSC for human
# This version matches what DGL originally used to make the BAM files
# Note: Annotation data must first be downloaded from online interface per vignette
pepfasta = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/hg19_GRCh37_proseq.fasta"
CDSfasta = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Make_FASTA_customProDB/Human/hg19_GRCh37_codingseq.fasta"
PrepareAnnotationRefseq(genome='hg19', CDSfasta, pepfasta, annotation_path=annotation_path_hs, 
                         splice_matrix = FALSE, dbsnp = NULL, COSMIC = FALSE)
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
run_customProDB(annotation_path= annotation_path_hs, outfile="multiple_sample_CustomProDB.fasta", singleSample=FALSE, correct_chr_name=FALSE,
                path_to_sample="~/Work/1_Milk/RNA-Seq_Guided_Proteomics/Reads/Human/")

# # # # # # # # 
#
# Generate fasta from BAM files 
# with variant calling
#
# # # # # # # # 
vcffile = "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/VariantCalling/macaque_var.flt.vcf"
vcf = InputVcf(vcffile)

# vcf[[1]][1:3] # pull an example range of variants
table(values(vcf[[1]])[['INDEL']])
table(values(vcf[[1]])[['IDV']])
table(values(vcf[[1]])[['IDV']])

# Headers in VCF file
# INDEL       IDV       IMF        DP       VDB       RPB       MQB       BQ
# MQSB       SGB      MQ0F      I16       QS        PL      PL.1      PL.2
## TODO what is a VCF file supposed to look like?? Is this accurate?


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











