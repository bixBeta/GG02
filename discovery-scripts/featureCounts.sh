for i in *bam

do
	iSUB=`echo $i | cut -d . -f1` 	
	featureCounts -a TagDirs/peaks.OUT/SAF3 -F SAF -o $iSUB.readCountinPeaks.txt $i

done

