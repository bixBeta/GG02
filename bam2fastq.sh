#!/bin/bash

for i in *.bam
do

iSUB=`basename $i .bam`
bedtools bamtofastq -i $i -fq $iSUB.R1.fq -fq2 $iSUB.R2.fq

done



