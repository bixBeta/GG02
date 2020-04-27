#!/bin/sh

#SBATCH -J fq.screen
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000


if [ "$1" = "help" ] || [ -z "$1" ]
then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "  sbatch" $0 
    echo "--------------------------------------------------------------------------------------"
    echo " CONFIGURATION FILE "
    echo "--------------------------------------------------------------------------------------"
    echo "THREADS         16"
    echo "DATABASE        Human /workdir/genomes/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38"
    echo "DATABASE        Mouse /workdir/genomes/FastQ_Screen_Genomes/Mouse/Mus_musculus.GRCm38"
    echo "DATABASE        Rat   /workdir/genomes/FastQ_Screen_Genomes/Rat/Rnor_6.0"
    echo "DATABASE        Drosophila    /workdir/genomes/FastQ_Screen_Genomes/Drosophila/BDGP6"
    echo "DATABASE        Worm  /workdir/genomes/FastQ_Screen_Genomes/Worm/Caenorhabditis_elegans.WBcel235"
    echo "DATABASE        Yeast /workdir/genomes/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1"
    echo "DATABASE	Arabidopsis_thaliana  /workdir/genomes/FastQ_Screen_Genomes/Arabidopsis/Arabidopsis_thaliana.TAIR10"
    echo "DATABASS	Ecoli /workdir/genomes/FastQ_Screen_Genomes/E_coli/Ecoli"
    echo "DATABASE  	MT    /workdir/genomes/FastQ_Screen_Genomes/Mitochondria/mitochondria"
    echo "DATABASE        PhiX  /workdir/genomes/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs"
    echo "DATABASE	Lambda    /workdir/genomes/FastQ_Screen_Genomes/Lambda/Lambda"
    echo "DATABASE        Vectors   /workdir/genomes/FastQ_Screen_Genomes/Vectors/Vectors"
    echo "DATABASE        Adapters  /workdir/genomes/FastQ_Screen_Genomes/Adapters/Contaminants"
    echo "DATABASE        Chicken   /workdir/genomes/FastQ_Screen_Genomes/Chicken/genome"
    echo "DATABASE        Dog   /workdir/genomes/FastQ_Screen_Genomes/Dog/genome"
    echo "DATABASE        Horse /workdir/genomes/FastQ_Screen_Genomes/Horse/genome"
    echo "DATABASE        Archae    /workdir/genomes/FastQ_Screen_Genomes/Archae/archae"
    echo "DATABASE        Bacteria  /workdir/genomes/FastQ_Screen_Genomes/Bacteria/bacteria"
    echo "DATABASE        Virus /workdir/genomes/FastQ_Screen_Genomes/Virus/viral"
    echo "DATABASE        Fungi /workdir/genomes/FastQ_Screen_Genomes/Fungi/fungi"
    echo "DATABASE        Protists  /workdir/genomes/FastQ_Screen_Genomes/Protists/protists"
    echo "DATABASE        AllrRNA   /workdir/genomes/contaminants/SILVA_rRNA/Bowtie2Index/AllrRNA"
    echo "DATABASE        EHV8  /workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8"
    echo ""

    exit 1

else
        for i in *.gz
        do
          fastq_screen --outdir fq.screen_out --conf /home/fa286/bin/scripts/my.fastq.conf $i
        done

        cd fq.screen_out
        multiqc -n fq.screen.multiqc.report .

        cd ..
fi

