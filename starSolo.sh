#!/bin/sh

usage(){

    echo "Usage: bash" $0 "[-h] [-p arg] [-g arg] [-r arg]"
    echo
    echo "-------------------------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help"
    echo "[-p] --> Project Identifier Number"
    echo "[-g] --> Reference Genome < hg38, GRCh38, mm10, GRCm38, rat, cat, chicken, horse, ATCC_13047, grape, ercc, lonchura, goose >"
    echo "[-r] --> <SE>, <SES>, <PE> or <PES> "
    echo "-------------------------------------------------------------------------------------------------------------------------------"
}

se(){
      cd trimmed_fastqs

      for i in *_trimmed.fq.gz

      do

        iSUB=`echo $i | cut -d '_' -f5`

        STAR \
        --runThreadN 12 \
        --genomeDir ${genomeDir[${DIR}]} \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts

      done

            source activate RSC
            multiqc -f -n ${PIN}.starSolo.multiqc.report .
            mkdir STAR.SOLO.COUNTS STAR.SOLO.BAMS STAR.SOLO.LOGS
            mv *.ReadsPerGene.out.tab STAR.SOLO.COUNTS
            mv *.bam STAR.SOLO.BAMS
            mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SOLO.LOGS
            mkdir STAR.SOLO
            mv STAR.* *.html STAR.SOLO
            mv STAR.SOLO ..
            cd ..
}

se_split(){

        cd trimmed_fastqs

        for i in *_trimmed.fq.gz

        do

          iSUB=`echo $i | cut -d '_' -f5`

          STAR \
          --runThreadN 12 \
          --genomeDir ${genomeDir[${DIR}]} \
          --readFilesIn $i \
          --readFilesCommand gunzip -c \
          --outSAMstrandField intronMotif \
          --outFilterIntronMotifs RemoveNoncanonical \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix $iSUB. \
          --limitBAMsortRAM 61675612266 \
          --quantMode GeneCounts

        done

                source activate RSC
                multiqc -f -n ${PIN}.starSolo.multiqc.report .
                mkdir STAR.SOLO.COUNTS STAR.SOLO.BAMS STAR.SOLO.LOGS STAR.SOLO.Unmapped
                mv *Unmapped.out.mate* STAR.SOLO.Unmapped
                cd STAR.SOLO.Unmapped
                    for i in *mate*
                        do
                            mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                        done
                cd ..
                mv *.ReadsPerGene.out.tab STAR.SOLO.COUNTS
                mv *.bam STAR.SOLO.BAMS
                mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SOLO.LOGS
                mkdir STAR.SOLO
                mv STAR.* *.html STAR.SOLO
                mv STAR.SOLO ..
                cd ..

}

pe(){

          cd trimmed_fastqs
          ls -1 *_R1_val_1.fq.gz > .trR1
          ls -1 *_R2_val_2.fq.gz > .trR2
          paste -d " " .trR1 .trR2 > Trimmed.list

          readarray trimmedFastqs < Trimmed.list

          for i in "${trimmedFastqs[@]}"

          do

            iSUB=`echo $i | cut -d '_' -f5`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

                    source activate RSC
                    multiqc -f -n ${PIN}.starSolo.multiqc.report .
                    mkdir STAR.SOLO.COUNTS STAR.SOLO.BAMS STAR.SOLO.LOGS
                    mv *.ReadsPerGene.out.tab STAR.SOLO.COUNTS
                    mv *.bam STAR.SOLO.BAMS
                    mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SOLO.LOGS
                    mkdir STAR.SOLO
                    mv STAR.* *.html STAR.SOLO
                    mv STAR.SOLO ..
                    cd ..
}

pe_split(){

          cd trimmed_fastqs
          ls -1 *_R1_val_1.fq.gz > .trR1
          ls -1 *_R2_val_2.fq.gz > .trR2
          paste -d " " .trR1 .trR2 > Trimmed.list

          readarray trimmedFastqs < Trimmed.list

          for i in "${trimmedFastqs[@]}"

          do

            iSUB=`echo $i | cut -d '_' -f5`

            STAR \
            --runThreadN 12 \
            --genomeDir ${genomeDir[${DIR}]} \
            --readFilesIn $i \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${iSUB}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

          done

                    source activate RSC
                    multiqc -f -n ${PIN}.starSolo.multiqc.report .
                    mkdir STAR.SOLO.COUNTS STAR.SOLO.BAMS STAR.SOLO.LOGS STAR.SOLO.Unmapped
                    mv *Unmapped.out.mate* STAR.SOLO.Unmapped
                    cd STAR.SOLO.Unmapped
                        for i in *mate*
                            do
                                mv $i `echo $i | sed "s/Unmapped/not.$DIR/g"`
                            done
                    cd ..							
                    mv *.ReadsPerGene.out.tab STAR.SOLO.COUNTS
                    mv *.bam STAR.SOLO.BAMS
                    mv *.out *.tab *_STARtmp *.list *.multiqc.report_data STAR.SOLO.LOGS
                    mkdir STAR.SOLO
                    mv STAR.* *.html STAR.SOLO
                    mv STAR.SOLO ..
                    cd ..
}

declare -A genomeDir

genomeDir=( ["hg38"]="/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star" \
["mm10"]="/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star" \
["GRCh38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/GRCh38.star" \
["GRCm38"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/GRCm38.star" \
["cat"]="/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir" \
["chicken"]="/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star" \
["horse"]="/workdir/genomes/Equus_caballus/ENSEMBL/Equus_caballus.star" \
["ATCC_13047"]="/workdir/genomes/Enterobacter_cloacae/ATCC_13047/custom/ATCC_13047.GTF" \
["grape"]="/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.12X.43.bed12" \
["rat"]="/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/rat.star" \
["ercc"]="/workdir/genomes/contaminants/ERCC_spikeIns/ercc.star" \
["lonchura"]="/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/lonchura.star" \
["goose"]="/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/goose.star" )


while getopts "hp:g:r:" opt; do
case ${opt} in

    h )
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

    r )

    RUN=$OPTARG

    ;;

    g )

    DIR=$OPTARG

    ;;


    \? )
    echo
    echo
    echo
    usage

    ;;

esac

done

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if genomeDir provided

if [[ ! -z "${DIR+x}" ]]; then
    if [ ${genomeDir[${DIR}]+_} ]; then
        echo Reference genome selected = $DIR
        echo

        if [[ ! -z "${RUN+x}" ]] && [[ $RUN == "PE" ]]; then
                pe
        elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "PES" ]]; then
            pe_split

    elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "SE" ]]; then
      se

    elif [[ ! -z "${RUN+x}" ]] && [[ $RUN == "SES" ]]; then
      se_split

    else
            echo "missing -r option "
            usage
            exit 1
        fi

    else
        echo "The reference genome provided '"$DIR"' is not available"
        echo " OR missing -r "
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
    echo >> solo.run.log
    echo -------------------------------------------------------------------------------------------------- >> solo.run.log
    echo `date` >> solo.run.log
    echo "Project Identifier Specified = " $PIN >> solo.run.log
    echo "Reference Genome Specified   = " $DIR >> solo.run.log
    echo "SE, SES, PE or PES           = " $RUN >> solo.run.log
    echo >> solo.run.log
    echo -------------------------------------------------------------------------------------------------- >> solo.run.log


fi
