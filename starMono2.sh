#!/bin/sh
if [ "$1" = "help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "<fastqFile> <genomeDirPath>"
    echo ""
    echo " NOTE: this script retains unmapped reads + spits unmapped fastq's"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

else
#--outSAMunmapped Within \


    iSUB=`echo $1 | cut -d "_" -f5`

        STAR \
        --runThreadN 12 \
        --genomeDir $2 \
        --readFilesIn $1 \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outReadsUnmapped Fastx \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${iSUB}. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts


fi
