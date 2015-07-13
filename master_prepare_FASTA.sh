#!/bin/bash

FILE=example.fasta
FROOT="$(echo $FILE | cut -d "." -f 1)"

decoy_fasta.pl $FILE > $FILE+decoy.fasta
cat $FILE+decoy.fasta /share/milklab/proteomics/prepare_FASTAs/cRAP.fasta > $FILE+decoy+cRAP.fasta
