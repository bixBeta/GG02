for i in *.bam
do
    iSUB=`echo $i | cut -d "." -f1`

    java -jar /programs/bin/picard-tools/picard.jar \
    CollectInsertSizeMetrics \
          I=$i \
          O=${iSUB}_insert_size_metrics.txt \
          H=${iSUB}_insert_size_histogram.pdf \
          M=0.5
done

mkdir insertSizeDists
mv *.txt *.pdf insertSizeDists

cd insertSizeDists
    pdfunite *.pdf insertSizes.pdf
cd ..


