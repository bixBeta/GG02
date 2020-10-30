#!/bin/sh

if [[ $1 = "help" ]] || [[ -z $1 ]]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "<hg38, mm10, cat or chicken> <PIN>"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

elif [[ $1 = "hg38" ]]; then

    BED12='/workdir/genomes/Homo_sapiens/hg38/UCSC/genes.bed12'

    for i in *.gz

    do
        DIR='/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star'
        iSUB=`echo $i | cut -d '_' -f5`

        STAR \
        --runThreadN 12 \
        --genomeDir $DIR \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts

    done

elif [[ $1 = "mm10" ]]; then

    BED12='/workdir/genomes/Mus_musculus/mm10/UCSC/BED12/mm10.ucsc.bed12'

    for i in *_trimmed.fq.gz

    do
        DIR='/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star'
        iSUB=`echo $i | cut -d '_' -f5`

        STAR \
        --runThreadN 12 \
        --genomeDir $DIR \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts

    done

elif [[ $1 = "cat" ]]; then

    BED12='/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/Felis_catus.Felis_catus_9.0.95.bed12'

    for i in *_trimmed.fq.gz

    do
        DIR='/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir'
        iSUB=`echo $i | cut -d '_' -f5`

        STAR \
        --runThreadN 12 \
        --genomeDir $DIR \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts

    done

elif [[ $1 = "chicken" ]]; then

    BED12='/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/Gallus_gallus.Gallus_gallus-5.0.bed12'

    for i in *_trimmed.fq.gz

    do
        DIR='/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star'
        iSUB=`echo $i | cut -d '_' -f5`

        STAR \
        --runThreadN 12 \
        --genomeDir $DIR \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts

    done

else
    echo ""
    echo "Please select appropriate genome Index"
    echo ""
    exit 1

fi


#####################################################################################################
PIN=$2
source activate RSC

    multiqc -f -n ${PIN}.star.multiqc.report .

    mkdir STAR.COUNTS STAR.BAMS STAR.LOGS
    mv *.ReadsPerGene.out.tab STAR.COUNTS
    mv *.bam STAR.BAMS
    mv *.out *.tab *_STARtmp *.list *star.multiqc.report_data STAR.LOGS

    cd STAR.BAMS

    for i in *.bam

        do
        samtools index -b $i
        done

    cd ..

    geneBody_coverage.py -r $BED12 -i STAR.BAMS/ -o $PIN
    mkdir geneBodyCov
    mv *geneBodyCoverage.* log.txt geneBodyCov

