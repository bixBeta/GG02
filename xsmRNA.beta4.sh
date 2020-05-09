#!/bin/sh

# source /programs/bin/util/setup_mirdeep2.sh

usage(){

    echo "sm R N A - S E Q   W O R K F L O W - @bixBeta"
    echo
    echo

    echo "Usage: bash" $0 "[-h arg] [-p arg] [-t arg] [-g arg]"
    echo
    echo "---------------------------------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help "
    echo "[-p] --> Project Identifier Number "
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: -,2 (complex field example: 2 | tail -c 4)"
    echo "[-t] --> NextSeq run < yes, no, na > "
    echo "[-g] --> Mapper Genome < hsa, mmu, cel > "
    echo "[-c] --> CleanUP < yes or no > "
    echo "---------------------------------------------------------------------------------------------------------------------------------------"
}


trimSmall(){
    echo "trimSmall"
        mkdir TrimQC_stats fastQC mirDeep2_results
        for i in fastqs/*.gz
        do
            /home/fa286/bin/TrimGalore-0.6.0/trim_galore --nextseq 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done
        mv *_trimming_report.txt TrimQC_stats
        mv *trimmed.fq.gz mirDeep2_results

}

trimHiSeq(){
    echo "trimHiSeq"
        mkdir TrimQC_stats fastQC mirDeep2_results
        for i in fastqs/*.gz
        do
            /home/fa286/bin/TrimGalore-0.6.0/trim_galore --quality 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done
        mv *_trimming_report.txt TrimQC_stats
        mv *trimmed.fq.gz mirDeep2_results

}

fastq2fasta(){
  echo "fastq2fasta"
  cd mirDeep2_results
  gunzip *.gz

  for i in *.fq
  do
    iSUB=`echo $i | cut -d "." -f1`
    fastq2fasta.pl $i > ${iSUB}.fasta
  done
  cd ..
}

config(){
  echo "config"
  cd mirDeep2_results
  ls -1 *.fasta > f1
  readarray fastas < f1


  for i in "${fastas[@]}"
    do

                if  echo $DELIM | grep -q "|"

                then
                echo $i | cut -d ${DELIMITER} -f${FIELD} | ${CCOUNT} >> f2

                else
                echo $i | cut -d ${DELIMITER} -f${FIELD} >> f2

                fi

    done

  paste f1 f2 > config.txt

  CONFIG="config.txt"
  cd ..

}

mapper(){
  echo "mapper"
  cd mirDeep2_results
  DATE=`date +"%m_%d_%H-%M"`
  mapper.pl $CONFIG -d -c -m -s ${PIN}_${DATE}.collapsed.fa
  cd ..
}

quant(){
  echo "quant"
  cd mirDeep2_results
  quantifier.pl -p /workdir/RSC/referenceFiles/miRBase/v22_1/hairpin.fa \
  -m /workdir/RSC/referenceFiles/miRBase/v22_1/mature.fa \
  -t $G -y ${PIN}_${DATE} -r ${PIN}_${DATE}.collapsed.fa -W -d
  cd ..
}


cleanUp(){
  echo "cleanUp"
  cd mirDeep2_results
  rm -r dir_mapper* f1 f2 *_trimmed.fasta *_trimmed.fq
  cd expression_analyses/expression_analyses_${PIN}_${DATE}
      mv *.mrd *.arf ../../
  cd ../../

  mv miRBase.mrd ${PIN}_${DATE}_miRBase.mrd
  mv mature_mapped.arf ${PIN}_${DATE}_mature_mapped.arf
  rm -r expression_analyses
  mkdir  expression_analyses_${PIN}_${DATE}
  mv *.arf *.fa *.mrd *.html *.csv  expression_analyses_${PIN}_${DATE}

  cd expression_analyses_${PIN}_${DATE}
    gzip *.arf *.fa *.mrd
  cd ..

  echo ""
  echo "DONE =) "

}

# PIN="1058"
# fastq2fasta
# config
# mapper


while getopts "hp:t:g:d:c:" opt; do
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


  g )

    G=$OPTARG

  ;;

    d )

    DELIM=$OPTARG

  ;;

    c )

        CLEAN=$OPTARG

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
## check if PIN is provided

if [[ -z "${PIN+x}" ]]; then

    PIN="PIN_Null"
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if delimiter parameter exists

if [[ ! -z "${DELIM+x}" ]]; then
    #statements
    if [[ $DELIM == default ]]; then

    DELIMITER="-"
    FIELD="2"
    echo "file naming will be done using the default delimiter settings"

    else

        DELIMITER=`echo $DELIM | cut -d , -f1`
        FIELD=`echo $DELIM | cut -d , -f2- | cut -d "|" -f1`
        CCOUNT=`echo $DELIM | cut -d , -f2- | cut -d "|" -f2-`

    echo "file naming will be done using the delim = $DELIMITER and field = $FIELD initial params for $DELIM"

  fi

fi
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if trimming parameter exists and run on nextseq 500 series (or 2 color bias )

if [[ ! -z "${T+x}" ]]; then
    #statements

    if [[ $T == yes ]]; then
        trimSmall
    if [[ $T == yes ]]; then
        trimSmall
        fastq2fasta

  else
    trimHiSeq
        fastq2fasta

    fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if genome parameter is  provided

if [[ ! -z "${G+x}" ]]; then
    #statements
  echo "Genome selected --> $G "
  config
  mapper
  quant

elif [[ -z "$G"  ]]; then
  echo "Genome info not provided "
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if clean up parameter exists

if [[ ! -z "${CLEAN+x}" ]]; then
    #statements

    if [[ $CLEAN == yes ]]; then
        cleanUp

  else
    echo 'cleanUP not required at this time'

    fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

if [[ -z $1 ]] || [[  $1 = "--help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    echo
    exit 1
else
    echo
    echo `date` >> beta4.small.run.log
    echo "Project Identifier Specified = " $PIN >> beta4.small.run.log
    echo "Trimming for NextSeq         = " $T >> beta4.small.run.log
    echo "Selected Genome              = " $G >> beta4.small.run.log
    echo -------------------------------------------------------------------------------------------------- >> beta4.small.run.log
fi
