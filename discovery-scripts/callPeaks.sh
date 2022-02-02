for i in *tagDir
do

    iSUB=` basename $i .tagDir`

    /home/fa286/bin/HOMER/bin/findPeaks $i -style tss > $iSUB.tss.peaks.out.txt

done
