#!/bin/sh

source ~/.bash_profile

mkdir fastQC


for i in *.gz
do

$TRIM --nextseq 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i 

done


mkdir trimmed-fastq TrimQC_Stats
mv *trimmed.fq.gz trimmed-fastq
mv *trimming_report.txt TrimQC_Stats
