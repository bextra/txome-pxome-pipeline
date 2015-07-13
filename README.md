# txome-pxome-pipeline

# Objective
Gernerate fasta files that are derived from RNA-Seq Data for use in peptide-to-protein assignment process

# Generate FASTA Databases from Different Sources
1. De novo transcriptome

2. Genome-guided transcriptome

3. Filter on RPKM with CustomProDB
Note:
SNV proteins are annotated in the fasta headed as where the single amino acid change is indicated:
>NP_000009_P65L
Indels are annotated in the fasta header as where the insertion or deletion is indicated:
>NP_001157750_527:ACC>AC

4. Generate variant called database with CustomProDB

