for i in *bam
do

iSUB=`basename $i .bam`

#/home/fa286/bin/HOMER/bin/makeTagDirectory $iSUB.tagDir $i -genome /workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.dna.toplevel.fa -checkGC -keepOne
/home/fa286/bin/HOMER/bin/makeTagDirectory $iSUB.tagDir $i -genome /workdir/genomes/Mus_musculus/mm10/UCSC/genome.fa -checkGC -keepOne


done



