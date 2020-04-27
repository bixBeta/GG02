#!/bin/sh

#SBATCH -J fq.screen
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000


if [ "$1" = "help" ] || [ -z "$1" ]
then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "  sbatch "$0"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

else
        for i in *.gz
        do
          fastq_screen --outdir fq.screen_out --conf /home/fa286/bin/scripts/my.fastq.conf $i
        done

        cd fq.screen_out
        multiqc -n fq.screen.multiqc.report .

        cd ..
fi

