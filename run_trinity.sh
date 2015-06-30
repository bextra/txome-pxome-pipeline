#!/bin/bash

#Trinity --seqType fq --single monkey_2_1_index7_trimmed.fastq.gz --max_memory 8G --CPU 3
Trinity --seqType fq --single /share/milklab/proteomics/input_Trinity/monkey_2_1_index18_trimmed.fastq.gz --max_memory 8G --CPU 3 --output ./trinity_out_m2_1_i18 --full_cleanup --verbose

