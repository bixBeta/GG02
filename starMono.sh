#!/bin/bash

#SBATCH -J RNAseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000
if [ "$1" = "help" ] || [  -z $1  ]; then
  echo ""
  echo "--------------------------------------------------------------------------------------------------------"
  echo "  To run this script, use the following syntax:"
  echo "     bash" $0 "<delim> <hg38, GRCh38, mm10 etc.>"
  echo "--------------------------------------------------------------------------------------------------------"
  echo ""
  echo "******************************************** "
  echo "**********/ EXTENDED GENOME LIST /********** "
  echo "[hg38]=/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star "
  echo "[mm10]=/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star "
  echo "[GRCh38]=/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/GRCh38.star "
  echo "[GRCm38]=/workdir/genomes/Mus_musculus/mm10/ENSEMBL/GRCm38.star "
  echo "[cat]=/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir "
  echo "[chicken]=/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star "
  echo "[horse]=/workdir/genomes/Equus_caballus/ENSEMBL/Equus_caballus.star "
  echo "[rat]=/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/rat.star "
  echo "[ercc]=/workdir/genomes/contaminants/ERCC_spikeIns/ercc.star "
  echo "[lonchura]=/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/lonchura.star "
  echo "[goose]=/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/goose.star "
  echo "[ehv8]=/workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8.star "
  echo "[erdman]=/workdir/genomes/Mycobacterium_tuberculosis/Ensembl_GCA_000668235/GCA_000668235.star "
  echo "[TB]=/workdir/genomes/Mycobacterium_tuberculosis/CDC1551_Ensembl/cdc1551.star "
  echo "[maize]=/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/star.maize "
  echo "[finch]=/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/star.index "
  echo "[dog]=/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/star.index "
  echo ""
  echo ""
  exit 1

else


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
["goose"]="/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/goose.star" \
["ehv8"]="/workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8.star" \
["erdman"]="/workdir/genomes/Mycobacterium_tuberculosis/Ensembl_GCA_000668235/GCA_000668235.star" \
["TB"]="/workdir/genomes/Mycobacterium_tuberculosis/CDC1551_Ensembl/cdc1551.star" \
["maize"]="/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/star.maize" )


  ls -1 *mate1 > .R1
  ls -1 *mate2 > .R2

  paste -d " " .R1 .R2 > reads.list

  DELIMITER=`echo $1 | cut -d , -f1`
  FIELD=`echo $1 | cut -d , -f2-`


  readarray fastqs < reads.list

  for i in "${fastqs[@]}"

  do

    iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

    STAR \
    --runThreadN 12 \
    --genomeDir ${genomeDir[${2}]} \
    --readFilesIn $i \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${iSUB}.${2}.hits. \
    --limitBAMsortRAM 61675612266 \
    --quantMode GeneCounts

  done

fi
