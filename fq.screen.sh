#!/bin/sh

#SBATCH -J fq.screen
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

for i in *.gz
do
  fastq_screen --outdir fq.screen_out --conf /home/fa286/bin/scripts/my.fastq.conf $i
done

cd fq.screen_out
multiqc -n fq.screen.multiqc.report .

cd ..
