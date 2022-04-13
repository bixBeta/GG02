#!/bin/bash

#SBATCH -J fq.screen
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

usage(){

  echo "FASTQ__SCREEN   W O R K F L O W - @bixBeta"
  echo ""
  echo ""
  echo "Usage: bash" $0 "[-h arg] [-p arg] [-r arg] "
  echo
  echo "---------------------------------------------------------------------------------------------------------------------------"
  echo "[-h] --> Display Help"
  echo "[-p] --> Project Identifier Number"
  echo "[-r] --> nova , next "
  echo "---------------------------------------------------------------------------------------------------------------------------"
  echo ""
  echo "Configured DBs: "
}



declare -A db

db=( ["Human"]="/workdir/genomes/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38"
["Mouse"]="/workdir/genomes/FastQ_Screen_Genomes/Mouse/Mus_musculus.GRCm38"
["Rat"]="/workdir/genomes/FastQ_Screen_Genomes/Rat/Rnor_6.0"
["Drosophila"]="/workdir/genomes/FastQ_Screen_Genomes/Drosophila/BDGP6"
["Worm"]="/workdir/genomes/FastQ_Screen_Genomes/Worm/Caenorhabditis_elegans.WBcel235"
["Yeast"]="/workdir/genomes/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1"
["Arabidopsis_thaliana"]="/workdir/genomes/FastQ_Screen_Genomes/Arabidopsis/Arabidopsis_thaliana.TAIR10"
["Ecoli"]="/workdir/genomes/FastQ_Screen_Genomes/E_coli/Ecoli"
["MT"]="/workdir/genomes/FastQ_Screen_Genomes/Mitochondria/mitochondria"
["PhiX"]="/workdir/genomes/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs"
["Lambda"]="/workdir/genomes/FastQ_Screen_Genomes/Lambda/Lambda"
["Vectors"]="/workdir/genomes/FastQ_Screen_Genomes/Vectors/Vectors"
["Adapters"]="/workdir/genomes/FastQ_Screen_Genomes/Adapters/Contaminants"
["Chicken"]="/workdir/genomes/FastQ_Screen_Genomes/Chicken/genome"
["Dog"]="/workdir/genomes/FastQ_Screen_Genomes/Dog/genome"
["Horse"]="/workdir/genomes/FastQ_Screen_Genomes/Horse/genome"
["Archae"]="/workdir/genomes/FastQ_Screen_Genomes/Archae/archae"
["Bacteria"]="/workdir/genomes/FastQ_Screen_Genomes/Bacteria/bacteria"
["Virus"]="/workdir/genomes/FastQ_Screen_Genomes/Virus/viral"
["Fungi"]="/workdir/genomes/FastQ_Screen_Genomes/Fungi/fungi"
["Protists"]="/workdir/genomes/FastQ_Screen_Genomes/Protists/protists"
["AllrRNA"]="/workdir/genomes/contaminants/SILVA_rRNA/Bowtie2Index/AllrRNA"
["EHV8"]="/workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8" )

printDB() {
  for i in "${!db[@]}"; do echo "[${i}]=${db[$i]}"; done
}


screen(){


          for i in $F
          do
            fastq_screen --outdir ${PIN}_fq.screen_out --conf /home/fa286/bin/scripts/my.fastq.conf $i
          done

          cd $PIN_fq.screen_out

          multiqc -n ${PIN}_fq.screen_multiqc.report .

          cd ..



}




while getopts "hp:r:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        printGenomes
        echo
        exit 1

    ;;

    p )

        PIN=$OPTARG
        echo "Project Identifier = " $PIN
    ;;

    r )

        PLATFORM=$OPTARG
        echo "Sequencing Platform = " $PLATFORM
    ;;


    \?)
        echo
        echo
        echo
        usage

    ;;

    esac

done




if [[ ! -z "${PLATFORM+x}" ]]; then
    #statements

    if   [[ $PLATFORM == next ]]; then

        F="*_R1.fastq.gz"
        screen

    elif [[ $PLATFORM == nova ]]; then

        F="*_val_1.fq.gz"
        screen

    else
        echo "-t only accepts next or nova as arguments"
        exit 1
    fi
fi



if [[ -z $1 ]] || [[  $1 = "help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    printDB
    echo
    exit 1

else
    echo
    echo `date`
    echo "Project Identifier Specified = " $PIN
    echo "Sequencing Platform Specified   = " $DIR
    echo

    echo "ENV INFO: "
    echo

    echo `fastq_screen --version`


    echo --------------------------------------------------------------------------------------------------

fi
