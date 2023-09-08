#!/bin/bash
#SBATCH -J ATACseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000



usage(){

    echo "A T A C - S E Q   W O R K F L O W - @bixBeta"
    echo
    echo

    echo "Usage: bash" $0 "[-h arg] [-p arg] [-d args] [-t arg] [-g arg] [-q arg] [-a arg] [-Q arg] [-F arg] "
    echo
    echo "---------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help"
    echo "[-p] --> Project Identifier Number"
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: _,1 "
    echo "[-t] --> Trimming <nextseq or nova>;"
    echo "[-g] --> Reference Genome <mm10, hg38, dm6>"
    echo "[-a] --> Reference Genome <bwa or bt2 >"
    echo "[-q] --> Execute atacQC.R script <yes>"
    echo ""
    echo "---------------------- MACS2 PARAMS ------------------------"
    echo "    [-Q] --> macs2 q val cutoff, default = 0.05"
    echo "    [-F] --> macs2 fold enrichment cutoff, default = 5"
    echo "---------------------------------------------------------------------------------------------------------------"
}



declare -A genomeDir

genomeDir=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/bwa.index/Homo_sapiens.GRCh38.dna.toplevel.fa" \
)

declare -A genomeDirBT2

genomeDirBT2=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/bowtie2index/mm10" \
["hg38"]="" \
["dm6"]="/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Bowtie2.Index/dm6"
)

declare -A gtfs

gtfs=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.96.gtf" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.96.gtf" \
["dm6"]="/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.32.106.gtf"
)

declare -A gAlias # for compatibility with atacQC.R
gAlias=(
["mm10"]="mouse" \
["hg38"]="human" \
["dm6"]="fly"
)

declare -A gSize # for macs2
gSize=(
["mm10"]="mm" \
["hg38"]="hs" \
["dm6"]="dm"
)

declare -A blkList
blkList=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/no-chr-mm10-blacklist.v2.bed" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/no-chr-hg38-blacklist.v2.bed" \
["dm6"]="/workdir/tools/blacklists/dm6-blacklist.v2.bed"
)


trimPE(){

        cd fastqs
        ls -1 *_R1.fastq* > .R1
        ls -1 *_R2.fastq* > .R2
        paste -d " " .R1 .R2 > Reads.list

        readarray fastqs < Reads.list
        mkdir fastQC

                for i in "${fastqs[@]}"
                do
                        #trim_galore --nextseq 20 --gzip -j 8 --length 50  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
                
                    A=`echo $i | cut -d " " -f1`
                    B=`echo $i | cut -d " " -f2`
                    iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`
                    

                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif fastp -z 4 -w 20 \
                    --length_required 50 --qualified_quality_phred 20 \
                    --trim_poly_g \
                    -i $A \
                    -I $B \
                    -o ${iSUB}_val_1.fq.gz \
                    -O ${iSUB}_val_2.fq.gz \
                    -h ${iSUB}.fastp.html \
                    -j ${iSUB}.fastp.json
                
                done

        mkdir TrimQC_stats trimmed_fastqs
        mv *_trimming_report.txt *json *html TrimQC_stats
        mv *_val* trimmed_fastqs
        mv TrimQC_stats fastQC trimmed_fastqs ../

        cd ..
}


trimNovaPE(){

                cd fastqs
                ls -1 *_1.fq* > .R1
                ls -1 *_2.fq* > .R2
                paste -d " " .R1 .R2 > Reads.list

                readarray fastqs < Reads.list
                mkdir fastQC

                for i in "${fastqs[@]}"
                do
                        #trim_galore --nextseq 20 --gzip -j 8 --length 50  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
                
                    A=`echo $i | cut -d " " -f1`
                    B=`echo $i | cut -d " " -f2`
                    iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`
                    

                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif fastp -z 4 -w 20 \
                    --length_required 50 --qualified_quality_phred 20 \
                    --trim_poly_g \
                    -i $A \
                    -I $B \
                    -o ${iSUB}_val_1.fq.gz \
                    -O ${iSUB}_val_2.fq.gz \
                    -h ${iSUB}.fastp.html \
                    -j ${iSUB}.fastp.json
                
                done

                mkdir TrimQC_stats trimmed_fastqs
                mv *_trimming_report.txt *json *html TrimQC_stats
                mv *_val* trimmed_fastqs
                mv TrimQC_stats fastQC trimmed_fastqs ..

                cd ..
}




alignPE.bwa(){

        cd trimmed_fastqs

        ls -1 *_val_1.fq.gz > .trR1
        ls -1 *_val_2.fq.gz > .trR2
        paste -d " " .trR1 .trR2 > Trimmed.list

        readarray trimmedFastqs < Trimmed.list

        for i in "${trimmedFastqs[@]}"

        do
                # INDEX="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa"
                iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

                bwa mem -t 24 -M -R "@RG\tID:${iSUB}\tSM:${iSUB}\tPL:ILLUMINA\tLB:${iSUB}\tPU:1" ${genomeDir[${DIR}]} $i \
                | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${iSUB}.bam
        done

                mkdir primary-BAMS
                mv *.bam primary-BAMS
                mv primary-BAMS ..
                cd ..
}

alignPE.bt2(){

        cd trimmed_fastqs

        ls -1 *_val_1.fq.gz > .trR1
        ls -1 *_val_2.fq.gz > .trR2
        paste -d " " .trR1 .trR2 > Trimmed.list

        readarray trimmedFastqs < Trimmed.list

          for i in "${trimmedFastqs[@]}"
            do

              iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`
              A=`echo $i | cut -d " " -f1`
              B=`echo $i | cut -d " " -f2`

              (/workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif bowtie2 \
              --no-unal \
              -x ${genomeDirBT2[${DIR}]} \
              -1 $A -2 $B \
              --threads 12 \
              -S - | /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${iSUB}.bam)2>${iSUB}.log

            done

                mkdir primary-BAMS
                mv *.bam primary-BAMS
                mv primary-BAMS ..
                cd ..
}


sort(){
    
            cd primary-BAMS
            for i in *.bam
            do
            /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools sort $i > `echo  $i | cut -d "." -f1`.sorted.bam
            done

            for i in *.sorted.bam
            do
                /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools index $i
            done

                    # alignment stats etc. on raw bams
                    for i in *.sorted.bam
                    do
                        iSUB=`echo $i | cut -d "." -f1`
                        /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools flagstat $i > ${iSUB}.primary.flagstat
                        /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools idxstats $i > ${iSUB}.primary.idxstats
                    done
                

                mkdir primary_stats
                mv *.primary.flagstat *.primary.idxstats primary_stats/
                cd ..
                pwd
}





rmMT(){
                cd primary-BAMS
                for i in *.sorted.bam
                do

                    iSUB=`echo $i | cut -d "." -f1`

                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools view -H `ls -1 *.sorted.bam | head -1` | cut -f2 | grep "SN:" |  cut -d ":" -f2 | grep -v "MT\|_\|\." | xargs /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools view -b $i > ${iSUB}.noMT.bam

                done

                for i in *.noMT.bam
                do
                    iSUB=`echo $i | cut -d "." -f1`
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools index $i
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools flagstat $i > ${iSUB}.noMT.flagstat
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools idxstats $i > ${iSUB}.noMT.idxstats
                done

                mkdir noMT_stats
                mv *noMT.flagstat *noMT.idxstats noMT_stats/

                for i in *.noMT.bam
                do
                    iSUB=`echo $i | cut -d "." -f1`
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif bedtools intersect -v -a $i -b ${blkList[${DIR}]} > ${iSUB}.noBlacklist.noMT.bam
                done

                # for i in *.noBlacklist.noMT.bam
                # do
                #     iSUB=`echo $i | cut -d "." -f1`
                #     /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools index $i
                #     /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools flagstat $i > ${iSUB}.noBlacklist.noMT.flagstat
                #     /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools idxstats $i > ${iSUB}.noBlacklist.noMT.idxstats
                # done

                # mkdir noBLK_stats
                # mv *.noBlacklist.noMT.flagstat *.noBlacklist.noMT.idxstats noBLK_stats/


                cd ..
}

markDups(){
                cd primary-BAMS

                for i in *.noBlacklist.noMT.bam
                do
                    iSUB=`echo $i | cut -d "." -f1`
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif java -jar /myBin/picard.jar \
                    MarkDuplicates \
                    INPUT=$i \
                    OUTPUT=${iSUB}.dupMarked.noBlacklist.noMT.bam \
                    ASSUME_SORTED=true \
                    REMOVE_DUPLICATES=false \
                    METRICS_FILE=${iSUB}.MarkDuplicates.metrics.txt \
                    VALIDATION_STRINGENCY=LENIENT \
                    TMP_DIR=tmp

                done
                cd ..
}

dedupBAM(){
                cd primary-BAMS
                # alignment stats etc. on dupMarked no MT bams
                for i in *.dupMarked.noBlacklist.noMT.bam
                do
                    iSUB=`echo $i | cut -d "." -f1`
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools index $i
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools flagstat $i > ${iSUB}.noBlacklist.noMT.flagstat
                    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools idxstats $i > ${iSUB}.noBlacklist.noMT.idxstats
                done

        for i in *.dupMarked.noBlacklist.noMT.bam
        do
                iSUB=`echo $i | cut -d "." -f1`
                /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools view -b -h -F 0X400 $i > ${iSUB}.DEDUP.bam
        done

                for i in *.DEDUP.bam; do /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools index $i ; /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools idxstats $i > `echo $i | cut -d "." -f1`.DEDUP.idxstats; done
                for i in *.DEDUP.bam; do /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools flagstat $i > `echo $i | cut -d "." -f1`.DEDUP.flagstat; done

                multiqc -n ${PIN}.multiqc.report .

                mkdir dedup-BAMS
                mv *.DEDUP* dedup-BAMS/
                mv dedup-BAMS ..
                cd ..

}

tagDir(){
  cd dedup-BAMS
  for i in *.DEDUP.bam
  do
  iSUB=`echo "$i" | cut -d'.' -f1` # subset to rename
  /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif makeTagDirectory "$iSUB".tag.dir "$i"
  done
  cd ..
}

callPeak(){

  cd dedup-BAMS
    echo "calling peaks on DEDUP bams"
    echo "qval cutoff = $QVAL"
    echo "fe cutoff = $FE"

    mkdir peaks.OUT
    #export PYTHONPATH=/programs/macs2-2.2.7.1-r9/lib64/python3.9/site-packages
    #export PATH=/programs/macs2-2.2.7.1-r9/bin:$PATH

    for i  in *.DEDUP.bam
      do
        iSUB=`echo $i | cut -d "." -f1`
        /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif macs2 callpeak -t $i \
        -f BAMPE \
        -n ${iSUB} \
        -g ${gSize[${DIR}]} \
        -q ${QVAL} \
        --outdir peaks.OUT \
        --nomodel --shift 37 --ext 73 \
        --fe-cutoff ${FE} \
        --keep-dup all
      done

  cd ..

}


mergedPeaks(){

  echo "running module mergedPeaks"
  echo "qval cutoff = $QVAL"
  echo "fe cutoff = $FE"

  cd dedup-BAMS

    allBams=`echo *.DEDUP.bam`

    /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif macs2 callpeak -t ${allBams} \
    -f BAMPE \
    -n allSamplesMergedPeakset \
    -g ${gSize[${DIR}]} \
    -q ${QVAL} \
    --outdir peaks.OUT \
    --nomodel --shift 37 --ext 73 \
    --fe-cutoff ${FE} \
    --keep-dup all

  cd ..
}


saf(){
  # awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${sample}_peaks.narrowPeak > ${sample}_peaks.saf
  cd dedup-BAMS/peaks.OUT
  awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' allSamplesMergedPeakset_peaks.narrowPeak > allSamplesMergedPeakset.saf
  cd ../../
}

#saf

frip(){
  # featureCounts -p -a ${sample}_peaks.saf -F SAF -o readCountInPeaks.txt ${sample}.sorted.marked.filtered.shifted.bam
  cd dedup-BAMS/
    for i  in *.DEDUP.bam
    do
      iSUB=`echo $i | cut -d "." -f1`
      /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif featureCounts -p -a peaks.OUT/allSamplesMergedPeakset.saf -F SAF -o "${iSUB}".readCountInPeaks.txt $i
    done

  cd ..
}

annotatePeaks(){

    cd dedup-BAMS/peaks.OUT
    /workdir/TREx_shared/projects/CHIP_ATAC_V3.sif annotatePeaks.pl allSamplesMergedPeakset.saf ${DIR} -gtf ${gtfs[${DIR}]} > allSamplesMergedPeakset.Annotated.saf
    cd ../..
}

bedGraphs(){
  cd dedup-BAMS
    for i in *.tag.dir
    do
       /workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif makeUCSCfile ${i} -o auto -fsize 1e10 -res 1 -color 106,42,73 -style chipseq
    done

    mkdir tagDirs
        mv *.tag.dir tagDirs
        cd tagDirs
        mkdir bedGraphs
            for i in *.tag.dir
            do
                cd $i
                zcat *.ucsc.bedGraph.gz | awk '{if(NR>1) print "chr"$0; else print $0}' | gzip > `basename *.ucsc.bedGraph.gz .ucsc.bedGraph.gz`.ucsc.bg.gz
                mv *.ucsc.bg.gz ../bedGraphs
                cd ..
            done
        cd ..

    mkdir featureCounts
    mv *.txt featureCounts

    multiqc -n ${PIN}.FRIP.multiqc.report -b "Please note that the featureCounts M Assigned Column refers to Fragments and Not Reads" --ignore tagDirs --ignore peaks.OUT .

  cd ..
}

atacQC(){

    cd dedup-BAMS
    echo "genome alias" = ${gAlias[${DIR}]}
    
    /workdir/TREx_shared/projects/CHIP_ATAC_V2.sif atacQC.R ${gAlias[${DIR}]}
    # ${gAlias[${DIR}]}
    /workdir/TREx_shared/projects/CHIP_ATAC_V2.sif html.atacQC.sh `echo ${PIN}_atacQC`

    cd ..

    /workdir/TREx_shared/projects/CHIP_ATAC_V2.sif tree > folder.structure

}


while getopts "hp:t:g:q:d:a:Q:F:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        echo
        exit 1

    ;;

    p )

        PIN=$OPTARG
        echo "Project Identifier = " $PIN
    ;;

    t )

        T=$OPTARG

    ;;

    g)

        DIR=$OPTARG

    ;;

    q)

        QC=$OPTARG

    ;;

    d)
        DELIM=$OPTARG

    ;;

    a)
        AL=$OPTARG

    ;;

    Q)
        QVAL=$OPTARG

    ;;

    F)
        FE=$OPTARG

    ;;

    \? )
        echo
        echo
        echo
        usage

    ;;

    esac

done
# shift $((OPTIND -1))

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if PIN is provided

if [[ -z "${PIN+x}" ]]; then

    PIN="PIN_Null"
fi

# PARAMETER CHECKS

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if delimiter parameter exists
                    if [[ ! -z "${DELIM+x}" ]]; then
                        #statements
                        if [[ $DELIM == default ]]; then

                        DELIMITER="_"
                        FIELD="1"
                        echo "file naming will be done using the default delimiter settings"
                      else

                        DELIMITER=`echo $DELIM | cut -d , -f1`
                        FIELD=`echo $DELIM | cut -d , -f2-`
                        echo "file naming will be done using the delim = $DELIMITER and field = $FIELD settings"

                      fi

                    fi

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if trimming parameter exists

                    if [[ ! -z "${T+x}" ]]; then

                        if [[ $T == nextseq ]]; then
                            trimPE
                        elif [[ $T == nova ]]; then
                            trimNovaPE
                        else
                        echo "-t only accepts nextseq or nova as arguments"
                        exit 1

                        fi
                    fi

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    # check if macs2 cutoffs are provided
                    if [[ -z "${QVAL+x}" ]]; then

                        QVAL=0.05
                    fi

                    if [[ -z "${FE+x}" ]]; then

                        FE=5
                    fi
                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if genomeDir provided

                    if [[ ! -z "${DIR+x}" ]]; then
                        if [ ${genomeDir[${DIR}]+_} ]; then
                            echo Reference genome selected = $DIR
                            echo

                            if [[ ! -z "${AL+x}" ]]; then

                                if [[ $AL == bwa ]]; then
                                  alignPE.bwa
                                  sort
                                  rmMT
                                  markDups
                                  dedupBAM
                                  callPeak
                                  mergedPeaks
                                  saf
                                  frip
                                  tagDir
                                  annotatePeaks
                                  bedGraphs

                                elif [[ $AL == bt2 ]]; then
                                #alignPE.bt2
                                #   sort
                                   rmMT
                                   markDups
                                   dedupBAM
                                   callPeak
                                   mergedPeaks
                                   saf
                                   frip
                                   tagDir
                                   annotatePeaks
                                   bedGraphs

                                else
                                    echo "Please specify the available genome aligner (bt2 or bwa)"
                                    exit 1
                                fi
                            fi


                        else
                            echo "The reference genome provided '"$DIR"' is not available"
                            exit 1

                        fi
                    fi

                    if [[ ! -z "${QC+x}" ]]; then

                        if [[ $QC == yes ]]; then
                            atacQC
                        else
                            echo "-q option only accepts yes as an argument"
                            exit 1
                        fi
                    fi






#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

if [[ -z $1 ]] || [[  $1 = "help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    echo
    exit 1

else
    echo >> 3.5.atac.log
    echo `date -u` >> 3.5.atac.log
    echo "Project Identifier Specified = " $PIN >> 3.5.atac.log
    echo "Reference Genome Specified   = " $DIR >> 3.5.atac.log
    echo "Trimming                     = " $T >> 3.5.atac.log
    echo "macs2 qval cutoff            = " $QVAL >> 3.5.atac.log
    echo "macs2 fe cutoff              = " $FE >> 3.5.atac.log
    echo >> 3.5.atac.log

    echo "ENV INFO: " >> 3.5.atac.log
    echo >> 3.5.atac.log
    echo "BT2 version:" `/workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif bowtie2 --version | head -1 | cut -d ' ' -f3` >> 3.5.atac.log
    echo "BWA version: 0.7.17-r1188" >> 3.5.atac.log
    echo "multiqc version:" `/workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif multiqc --version | cut -d ' ' -f3` >> 3.5.atac.log
    echo "samtools version:" `/workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif samtools --version | head -1 | cut -d ' ' -f2` >> 3.5.atac.log
    echo "macs2 version:" `/workdir/TREx_shared/projects/CHIP_ATAC_DEV.sif macs2 --version | cut -d ' ' -f2` >> 3.5.atac.log
    echo "HOMER version: 4.11.1" >> 3.5.atac.log
    echo "--------------------------------------------------------------------------------------------------" >> 3.5.atac.log

fi
