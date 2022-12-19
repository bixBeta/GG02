#!/bin/bash

#SBATCH -J gbcov
#SBATCH -o %x.out
#SBATCH -n 6
#SBATCH --mem-per-cpu=18000

usage(){

    echo "geneBodyCoverage - @bixBeta"
    echo
    echo

    echo "Usage: bash" $0 "[-h arg] [-p arg] [-g arg]"
    echo
    echo "-------------------------------------------------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help "
    echo "[-p] --> Project Identifier Number "
    echo "[-g] --> genome < GRCh38, GRCm38, horse > "

    echo "-------------------------------------------------------------------------------------------------------------------------------------------------------"
}


declare -A bed12

bed12=(	["hg38"]="/workdir/genomes/Homo_sapiens/hg38/UCSC/genes.bed12" \
["mm10"]="/workdir/genomes/Mus_musculus/mm10/UCSC/BED12/mm10.ucsc.bed12" \
["GRCh38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.bed12" \
["GRCm38"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.bed12" \
["cat"]="/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/Felis_catus.Felis_catus_9.0.95.bed12" \
["chicken"]="/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/Gallus_gallus.Gallus_gallus-5.0.bed12" \
["chicken6"]="/workdir/genomes/Gallus_gallus/Galgal6/ENSEMBL/gtf.102/Gallus_gallus.GRCg6a.102.bed12" \
["horse"]="/workdir/genomes/Equus_caballus/EquCab3/ENSEMBL/genebuild-102/Equus_caballus.EquCab3.0.102.bed12" \
["horse2"]="/workdir/genomes/Equus_caballus/EquCab2/ENSEMBL/Equus_caballus.EquCab2.94.BED12" \
["ATCC_13047"]="/workdir/genomes/Enterobacter_cloacae/ATCC_13047/GCF_000025565.1_ASM2556v1_genomic.bed12" \
["grape"]="/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.star" \
["rat"]="/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/Rattus_norvegicus.Rnor_6.0.bed12" \
["lonchura"]="/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/Lonchura_striata_domestica.LonStrDom1.bed12" \
["goose"]="/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/Anser_brachyrhynchus.ASM259213v1.bed12" \
["maize"]="/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/Zea_mays.B73_RefGen_v4.bed12" \
["finch"]="/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/UPDATED.ANNOTS/Taeniopygia_guttata.bTaeGut1_v1.p.bed12" \
["finch2"]="/workdir/genomes/Geospiza_fortis_ground_finch/GeoFor_1.0/NCBI/GCF_000277835.1_GeoFor_1.0_genomic.bed12" \
["dog"]="/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/CanFam3.1.101/Canis_lupus_familiaris.CanFam3.1.101.bed12" \
["aedes"]="/workdir/genomes/Aedes_aegypti/AaegL5.0/NCBI/GCF_002204515.2_AaegL5.0_genomic.bed12" \
["aedesVB"]="/workdir/genomes/Aedes_aegypti/AaegL5.0/Vectorbase/AaegyptiLVP_AGWG_release49/VectorBase-49_AaegyptiLVP_AGWG.BED12" \
["cow"]="/workdir/genomes/Bos_taurus/ARS-UCD1.2_GCA_002263795.2/ENSEMBL/Bos_taurus.ARS-UCD1.2.BED12" \
["BG8"]="/workdir/genomes/Methylomicrobium_album_BG8/ASM21427v3/NCBI/test/ncbi-genomes-2020-08-27/GCF_000214275.2_ASM21427v3_genomic.bed12" \
["ddSmed"]="/workdir/genomes/Schmidtea_mediterranea/dd_Smed_v6/NCBI/dd_Smed_v6.pcf.bed12" \
["yeast"]="/workdir/genomes/Saccharomyces_cerevisiae/R64-1-1_GCA_000146045.2/ENSEMBL/Saccharomyces_cerevisiae.R64-1-1.bed12" )

printBED() {
  for i in "${!bed12[@]}"; do echo "[${i}]=${bed12[$i]}"; done
}


geneBodyCov(){
        # cd STAR*/*.BAMS

        # for i in *.bam
        # do
        #   BASE=`basename $(echo $i) .Aligned.sortedByCoord.out.bam `
        #   mv $i ${BASE}.bam
        # done
        #
        # for i in *.bam
        # do
        #   /programs/samtools-1.15.1-r/bin/samtools index -b $i
        # done
        # cd ..
        # echo
        # echo
        # pwd
        # echo
        # echo
        source activate RSeQC
        geneBody_coverage.py -r ${bed12[${BED}]} -i *.BAMS/ -o ${PIN}
        mkdir geneBodyCov
        mv *geneBodyCoverage.* log.txt geneBodyCov
        cd ..
}


while getopts "hp:t:g:d:c:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo
        printBED
        echo
        echo
        exit 1

    ;;

    p )

        PIN=$OPTARG
        echo "Project Identifier = " $PIN
    ;;


    g )

        BED=$OPTARG
        echo "Bed file Selected = " $BED
        geneBodyCov
    ;;

    \? )
        echo
        echo
        echo
        usage

    ;;


    esac

done


if [[ -z $1 ]] || [[  $1 = "--help"  ]] ; then
    #statements
    echo
    echo
    echo
    usage
    echo
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo
    printBED
    echo
    echo
    exit 1

fi
