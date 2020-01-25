#!/bin/sh

if [ "$1" = "help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "<title>"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    exit 1

else

scp /home/fa286/bin/scripts/qc.atac.Rmd .

Rscript /home/fa286/bin/scripts/knit.atacQC.R $1

rm qc.atac.Rmd

fi

