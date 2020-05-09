#!/bin/sh

for i in *.bam
do
  iSUB=`echo $i | cut -d "_" -f1`
  P=`echo "'name="\"$iSUB.plus\"\'`
  M=`echo "'name="\"$iSUB.minus\"\'`

  bedtools genomecov -bg \
  -trackline -trackopts $P \
  -ibam $i \
  -strand + > $iSUB.plusStrand.bg

  bedtools genomecov -bg \
  -trackline -trackopts $M \
  -ibam $i \
  -strand - > $iSUB.minusStrand.bg

done



