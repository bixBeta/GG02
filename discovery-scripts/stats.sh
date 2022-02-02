for i in *.bam
do
    iSUB=`echo $i | cut -d "." -f1,2`
    samtools index $i
    samtools flagstat $i > ${iSUB}.flagstat
    samtools idxstats $i > ${iSUB}.idxstats
    
done
