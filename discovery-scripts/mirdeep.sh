fastq2fasta(){
  echo "fastq2fasta"
  gunzip *.gz

  for i in *.fq
  do
    iSUB=`echo $i | cut -d "." -f1`
    fastq2fasta.pl $i > ${iSUB}.fasta
  done
  
}





fastq2fasta

ls -1 *.fasta > f1
COUNTER=`wc -l f1 | cut -d " " -f1`
COUNTERC=`expr $COUNTER + 100 `
seq 101 1 $COUNTERC > f2

paste f1 f2 > config.txt
CONFIG="config.txt"

