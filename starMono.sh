#!/bin/sh
if [ "$1" = "help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "<fastqFile> <hg38, GRCh38, mm10, GRCm38, rat, cat etc. >"
    echo "--------------------------------------------------------------------------------------------------------"
    echo ""
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

iSUB=`echo $1 | cut -d "_" -f5`

STAR \
--runThreadN 12 \
--genomeDir ${genomeDir[$2]} \
--readFilesIn $1 \
--readFilesCommand gunzip -c \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $iSUB. \
--limitBAMsortRAM 61675612266 \
--quantMode GeneCounts


fi
