#----------------------------------------------------------------------------------------
# get sample names and sample ids

for i in */
do
	cd $i/outs
	DIR=`pwd`
	cat input_samplesheet.csv | cut -d "," -f3 | cut -d "-" -f2 | cut -d "_" -f3 | sed 1,2d >> ../../.ids
	cat input_samplesheet.csv | cut -d "," -f3 | sed 1,2d >> ../../.names
	cd ../../

done

pwd


#----------------------------------------------------------------------------------------
# samples with  multiple run folders; un-comment this for autodetecting samples

# sort .ids | uniq -d > .idoi  


#----------------------------------------------------------------------------------------
# get fastqs path info 

for i in  */
do
	cd $i/outs
	cat input_samplesheet.csv | cut -d "," -f3 | sed 1,2d | awk '{print "'$i'" "outs/fastq_path/ :" $0}' >> ../../.paths.info
	cd ../../

done

pwd

wc -l .idoi


#----------------------------------------------------------------------------------------
# cellranger count

readarray sampleIDs < .idoi

for i in "${sampleIDs[@]}"
do

	ID=$i
	#echo "SAMPLE_ID = $ID"
	
	GREP_PATH=$(grep `echo $i`$ .paths.info | cut -d ":" -f1  | xargs | sed -e 's/ /,/g')
	#echo "FASTQ_PATH = $GREP_PATH" 

	GREP_NAME=$(grep `echo $i`$ .names | xargs | sed -e 's/ /,/g')
	#echo  "SAMPLE_NAME = $GREP_NAME"


	/programs/cellranger-3.0.2/cellranger count --id=$ID \
  	--transcriptome=/workdir/singleCellData/10x_reference_files/refdata-cellranger-GRCh38-3.0.0/ \
  	--fastqs=$GREP_PATH \
  	--sample=$GREP_NAME \
  	--localcores 20 --localmem 250


done




rm .ids .names .paths.info 	

