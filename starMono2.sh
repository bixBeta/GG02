#!/bin/sh

#SBATCH -J RNAseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

if [ "$1" = "help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "<delim> <fastqFile> <genomeDirPath>"
    echo ""
    echo " NOTE: this script retains unmapped reads + spits unmapped fastq's"
    echo "--------------------------------------------------------------------------------------"
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
    exit 1

else
#--outSAMunmapped Within \

    DELIMITER=`echo $1 | cut -d , -f1`
    FIELD=`echo $1 | cut -d , -f2-`
    iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`


        STAR \
        --runThreadN 12 \
        --genomeDir $3 \
        --readFilesIn $2 \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outReadsUnmapped Fastx \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${iSUB}. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts


fi
