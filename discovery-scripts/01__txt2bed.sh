for i in *.txt; do
	#statements
	iSUB=`basename $i .mm10.tss.peaks.out.txt`


	grep -v "^#" $i | awk 'OFS ="\t" {print $2, $3, $4 , $5, $1}' > ${iSUB}.bed

done
