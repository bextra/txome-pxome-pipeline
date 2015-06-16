
# Goals rename files or move them into a folder that makes sense in preparation for Scaffold
  # Move .xml files into respective FASTA folder
  # 4 categories:
  # 1. reference Ensembl
  # 2. de novo transcriptome
  # 3. genome guided transcriptome
  # 4. customProDB filtered on RPKM (from Ensembl too)

library("stringr")

setwd("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/") # destination path
data_path   = "/Volumes/Moose/RNA-Seq-Guided-Proteomics/CorrectAnimal-mzML/" # includes bio-matched reference databases
# data_path = "/Volumes/Moose/RNA-Seq-Guided-Proteomics/All-mzML/"           # includes alignment with bio matched databases and bio unmatched databases
# data_path = "/Volumes/Moose/RNA-Seq-Guided-Proteomics/Scaffold-UnmatchedSpectra-All/mzML/"

# Create a look up table with database and destination information
# TODO update this with the reference database
db_lookup = data.frame("db_name"   = c("macaque_denovo_Sixpack_ORFs_gg+decoy+crap.fasta", "macaque_denovo_Sixpack_ORFs+decoy+cRAP.fasta", "monkey_2_1_index7_customProDB+decoy+crap.fasta", "Blast2GO_Macaque-decoy+cRAP.fasta"),
                       "dest_name" = c("~/Work/1_Milk/RNA-Seq_Guided_Proteomics/FASTA_with_ref_transcriptome/", "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/FASTA_denovo_transcriptome/", "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/FASTA_customProDB_std/", "~/Work/1_Milk/RNA-Seq_Guided_Proteomics/FASTA_BLAST2GO/")
)
# Currently this is not being used

# Get GPM, sample, and reference data base info from files
output_gpm     = list.files(data_path)
in_db_matched  = sapply(output_gpm, function(x) system(paste("grep -m 1 '<file type=\"peptide\" ' ", data_path, x, sep=""), intern=TRUE))
in_out_matched = system(paste("grep '<bioml' ", data_path, "*", sep=""), intern =TRUE)
# This prefix is in the header and exists uniquely in each file

# Process xml file excerpts for file and database info
current_loc = str_extract(in_out_matched, ".*GPM[0-9]+.xml") # file path to data file
fasta_desc  = str_replace(str_extract(in_db_matched, "/fasta/.+.fasta"), "/fasta/", "") # database used in peptide to protein process

# Filter for the macaque sample used for the database
# <file type="peptide" URL="../fasta/monkey_2_1_index7_customProDB+decoy+crap.fasta"/>
bio_sample = str_extract(in_out_matched, "MMU[0-9]+")
# bio_sample_unmatched = str_extract(in_out_matched, "GPM[0-9]+.mgf") # this needs to be completed for the unmatched spectra because scaffold renames the input file

# gpm_moving  = as.data.frame(cbind(current_loc, bio_sample, bio_sample_unmatched, fasta_desc)) # database with both
gpm_moving  = as.data.frame(cbind(current_loc, bio_sample, fasta_desc)) # database with both
# you don't really have to move any file's until you want to load to scaffold - can just process them all from the database

#-- Extract performance parameters from GPM output --#
# TheGPM output description: http://www.thegpm.org/TANDEM/api/operf.html
# We want modelling statistics as those will be the totals prior to refinement

# Report total fasta db size
# <note label="modelling, total proteins used">32080</note> 
# this is the number of protein entries in the database including decoys and cRAP proteins
db_size  = sapply(output_gpm, function(x) system(paste("grep -m 1 '<note label=\"modelling, total proteins used\"' ", data_path, x, sep=""), intern=TRUE))
db_size  = as.numeric(str_extract(db_size, "[0-9]+"))

# Calculate the percent spectra identified
spec_assigned  = sapply(output_gpm, function(x) system(paste("grep -m 1 '<note label=\"modelling, total spectra assigned\"' ", data_path, x, sep=""), intern=TRUE))
spec_assigned  = as.numeric(str_extract(spec_assigned, "[0-9]+"))

spec_used  = sapply(output_gpm, function(x) system(paste("grep -m 1 '<note label=\"modelling, total spectra used\"' ", data_path, x, sep=""), intern=TRUE))
spec_used  = as.numeric(str_extract(spec_used, "[0-9]+")) # total spectra per xml file is same regardless of database

unique_assigned  = sapply(output_gpm, function(x) system(paste("grep -m 1 '<note label=\"modelling, total unique assigned\"' ", data_path, x, sep=""), intern=TRUE))
unique_assigned  = as.numeric(str_extract(unique_assigned, "[0-9]+"))

percent_assigned = spec_assigned/spec_used*100

gpm_moving= cbind(gpm_moving, db_size, spec_assigned, spec_used, unique_assigned, percent_assigned)

# Display some valuable calculations
levels(gpm_moving$fasta_desc)

# DB size with description
# 211392 (gg) 246012 (de novo)  33382 (customProDB)  72820 (blast2go)

## subset full version for just reference database to compare with "wrong animal"
reference_only = subset(gpm_moving, gpm_moving$fasta_desc == "Blast2GO_Macaque-decoy+cRAP.fasta")


# Display summary statistics for each fasta type
sapply(levels(gpm_moving$fasta_desc), function(x) summary(gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == x & gpm_moving$bio_sample=="MMU34863")]))


# Extract vectors based on fasta description alone - all biological samples
# original_run            = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "Blast2GO_Macaque-decoy+cRAP.fasta")]
# genome_guided           = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "denovo_Sixpack_ORFs_m2_1_i18_gg+decoy+cRAP.fasta")]
# de_novo                 = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "denovo_Sixpack_ORFs_m2_1_i18+decoy+cRAP.fasta")]
# customPro_DB            = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "monkey_2_1_index18-customProDB+decoy+cRAP.fasta")]
# unmatched_genome_guided = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "macaque_denovo_Sixpack_ORFs_gg+decoy+crap.fasta")]
# unmatched_de_novo       = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "macaque_denovo_Sixpack_ORFs+decoy+cRAP.fasta")]
# unmatch_customProDB     = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "monkey_2_1_index7_customProDB+decoy+crap.fasta")]



# Extract vectors for MMU34683 that demonstrate matched and unmatched reference data
original_run            = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "Blast2GO_Macaque-decoy+cRAP.fasta" & gpm_moving$bio_sample=="MMU34863")]
genome_guided           = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "denovo_Sixpack_ORFs_m2_1_i18_gg+decoy+cRAP.fasta" & gpm_moving$bio_sample=="MMU34863")]
de_novo                 = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "denovo_Sixpack_ORFs_m2_1_i18+decoy+cRAP.fasta" & gpm_moving$bio_sample=="MMU34863")]
customPro_DB            = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "monkey_2_1_index18-customProDB+decoy+cRAP.fasta" & gpm_moving$bio_sample=="MMU34863")]
unmatched_genome_guided = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "macaque_denovo_Sixpack_ORFs_gg+decoy+crap.fasta" & gpm_moving$bio_sample=="MMU34863")]
unmatched_de_novo       = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "macaque_denovo_Sixpack_ORFs+decoy+cRAP.fasta" & gpm_moving$bio_sample=="MMU34863")]
unmatch_customProDB     = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "monkey_2_1_index7_customProDB+decoy+crap.fasta" & gpm_moving$bio_sample=="MMU34863")]




# Extract vectors for 2-step annotation where unmatched spectra are re-ran against the database
# original_run            = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "Blast2GO_Macaque-decoy+cRAP.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]
# genome_guided           = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "denovo_Sixpack_ORFs_m2_1_i18_gg+decoy+cRAP.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]
# de_novo                 = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "denovo_Sixpack_ORFs_m2_1_i18+decoy+cRAP.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]
# customPro_DB            = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "monkey_2_1_index18-customProDB+decoy+cRAP.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]
# unmatched_genome_guided = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "macaque_denovo_Sixpack_ORFs_gg+decoy+crap.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]
# unmatched_de_novo       = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "macaque_denovo_Sixpack_ORFs+decoy+cRAP.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]
# unmatch_customProDB     = gpm_moving$percent_assigned[which(gpm_moving$fasta_desc == "monkey_2_1_index7_customProDB+decoy+crap.fasta" & gpm_moving$bio_sample_unmatched=="GPM00500038272.mgf")]




# Plot varying databases to determine if there is a difference between biologically matched or unmatched samples
boxplot(original_run, genome_guided, de_novo, customPro_DB, unmatched_genome_guided, unmatched_de_novo, unmatch_customProDB,
        names = c("original run", "genome guided", "de novo", "customProDB", "unmatched genome guided", "unmatched de novo", "unmatched custom pro DB"),
        col = c("grey", "purple", "orange", "green", "purple", "orange", "green"),
        main = "Comparison of database performance for spetral matching",
        ylab = "Percent spectra assigned",
        xlab = "Database description"
        )

# Determine distribution of percent spectra aligned for all database types
par(lwd=3)
plot(density(original_run),
     main = "Comparison of database performance for spectral assignment in databases
     with matched or unmatched transcriptome-proteome data",
     ylab = "Density",
     xlab = "Percent Spectra Assigned",
     ylim = c(0,0.09)
     )
lines(density(genome_guided), col="purple")
lines(density(de_novo), col="orange")
lines(density(customPro_DB), col="green")
lines(density(unmatched_genome_guided), lty="dashed", col="purple")
lines(density(unmatched_de_novo), lty="dashed", col="orange")
lines(density(unmatch_customProDB), lty="dashed", col="green")
legend("topright", fill=c("black", rep(c("purple", "orange", "green"), each= 2)),
       legend=c("BLAST2GO", "genome guided", "unmatched genome guided", "de novo", "unmatched de novo", "custom Pro DB", "unmatched customProDB"), 
       lty=c("solid", rep(c("solid", "dashed"), 3)),
       horiz=FALSE
       )



# Determine distribution differences of transcriptome based methods alone - plot for DGL grant
par(lwd=3)
plot(density(genome_guided),
     main = "Comparison of database performance for spectral assignment in databases
     with matched or unmatched transcriptome-proteome data",
     col = "purple",
     ylab = "Density",
     xlab = "Percent Spectra Assigned",
     xlim = c(5,70),
     ylim = c(0,0.065)
)
lines(density(de_novo), col="orange")
lines(density(unmatched_genome_guided), lty="dashed", col="purple")
lines(density(unmatched_de_novo), lty="dashed", col="orange")
legend("topright", fill=c(rep(c("purple", "orange"), each= 2)),
       legend=c("genome guided", "unmatched genome guided", "de novo", "unmatched de novo"), 
       lty=c(rep(c("solid", "dashed"), 2)),
       horiz=FALSE
)




