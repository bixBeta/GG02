DIR="/workdir/genomes/Mus_musculus/rRNA/Bowtie2Index/rRNA"


ls -1 *_1.fq.gz > .trR1
ls -1 *_2.fq.gz > .trR2
paste -d " " .trR1 .trR2 > Trimmed.list

readarray trimmedFastqs < Trimmed.list
for i in "${trimmedFastqs[@]}"
do

        iSUB=`echo $i | cut -d "_" -f1,2`
        A=`echo $i | cut -d " " -f1`
        B=`echo $i | cut -d " " -f2`

        (bowtie2 \
        --no-unal \
        --un-conc-gz ${iSUB}_rRNA_FILTERED_%.fq.gz \
        -x ${DIR} \
        -1 $A -2 $B \
        -S - | /programs/samtools-1.9-r9/bin/samtools view -@ 24 -b -O BAM -o ${iSUB}.bam)2>${iSUB}.log
done


mkdir rRNA_FILTERED_FASTQS
mv *rRNA_FILTERED* rRNA_FILTERED_FASTQ
