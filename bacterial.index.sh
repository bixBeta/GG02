    STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir $3 \
    --genomeFastaFiles $1 \
    --genomeSAindexNbases 5 \
    --sjdbGTFfile $2 \
    --sjdbOverhang 100
