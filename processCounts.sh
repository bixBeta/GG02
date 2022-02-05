#!/bin/bash

if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<1> or <2>"
    echo "    where 1 is for first strand and 2 is for reverse strands, 0 for unstranded counts"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1


elif [ "$1" = "1" ]
    then
        for i in *ReadsPerGene.out.tab
        do
        awk 'NR > 4 {print $1 "\t" $3}' $i > $i.rawCounts
        done
        mkdir rawCounts
        mv *.rawCounts rawCounts

elif [ "$1" = "2" ]
    then
        for i in *ReadsPerGene.out.tab
        do
        awk 'NR > 4 {print $1 "\t" $4}' $i > $i.rawCounts
        done
        mkdir rawCounts
        mv *.rawCounts rawCounts

else
    for i in *ReadsPerGene.out.tab
    do
    awk 'NR > 4 {print $1 "\t" $2}' $i > $i.rawCounts
    done
    mkdir rawCounts
    mv *.rawCounts rawCounts
fi
