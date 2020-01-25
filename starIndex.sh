#!/bin/sh


#  starIndex.sh
#  RSC-alpha
#
#  Created by Faraz Ahmed on 3/22/19.
#  

if [ "$1" = "help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "<fastaFile> <gtfFile> <genomeDirPath>"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo " NOTE: add --sjdbGTFtagExonParentTranscript Parent for gff annotations"
    echo ""
    echo ""

    exit 1

else


    #use mkdir to make genomeDir before running the code!!!
    # running star to generate genome indexes

    STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir $3 \
    --genomeFastaFiles $1 \
    --sjdbGTFfile $2 \
    --sjdbOverhang 100 \
    --limitGenomeGenerateRAM 152003700778



fi
