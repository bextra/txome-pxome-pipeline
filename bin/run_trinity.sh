#!/bin/bash

#Trinity --seqType fq --single monkey_2_1_index7_trimmed.fastq.gz --max_memory 8G --CPU 3

export JAVA_HOME=/Chanlab/Packages/Java/jdk1.7.0_80/bin/java
export PATH="/Chanlab/Packages/Java/jdk1.7.0_80/bin/:$PATH"

Trinity --seqType fq --single /share/milklab/proteomics/input_Trinity/monkey_2_1_index18_trimmed.fastq.gz --max_memory 8G --CPU 3 --output ./trinity_out_m2_1_i18 --full_cleanup --verbose

