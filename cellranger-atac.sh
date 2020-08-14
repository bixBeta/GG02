#!/bin/sh

# initial commit FA200323
# version 1

if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<auto> , <no-dups> or <idoi> "
    echo " 	   	auto - autodetects duplicate sampleIDs "
    echo "		no-dups - autodetects all sampleIDs"
    echo "		idoi - expects user to provide a sample ID list named idoi"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

elif [[ "$1" = "auto" ]]; then


    #statements
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
        #----------------------------------------------------------------------------------------
        # samples with  multiple run folders;
        sort .ids | uniq -d > idoi
        #----------------------------------------------------------------------------------------
        # get fastqs path info
        for i in  */
        do
            cd $i/outs
            cat input_samplesheet.csv | cut -d "," -f3 | sed 1,2d | awk '{print "'$i'" "outs/fastq_path/ :" $0}' >> ../../.paths.info
            cd ../../
        done
        #----------------------------------------------------------------------------------------
        # cellranger count
        readarray sampleIDs < idoi
        for i in "${sampleIDs[@]}"
        do
            ID=$i
            GREP_PATH=$(grep `echo $i`$ .paths.info | cut -d ":" -f1  | xargs | sed -e 's/ /,/g')
            GREP_NAME=$(grep `echo $i`$ .names | xargs | sed -e 's/ /,/g')


            echo "SAMPLE_ID = $ID"
            echo "FASTQ_PATH = $GREP_PATH"
            echo "SAMPLE_NAME = $GREP_NAME"
            echo "-----------------------------------"
            echo ""
            /programs/cellranger-atac-1.2.0/cellranger-atac count --id=$ID \
                   --reference=/workdir/singleCellData/10x_reference_files/refdata-cellranger-atac-GRCh38-1.2.0 \
                   --fastqs=$GREP_PATH \
                   --sample=$GREP_NAME \
                   --localcores 20 --localmem 250

        done

    DATE=`date +"%m_%d_%H-%M"`
    mkdir CellRangerATAC-Output_${DATE}
    readarray sampleIDs < idoi
      for i in "${sampleIDs[@]}"
      do
        mv $i CellRangerATAC-Output_${DATE}
      done

elif [[ "$1" = "no-dups" ]]; then


    #statements
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


        # get sample ids
        cat .ids > idoi

        #----------------------------------------------------------------------------------------
        # get fastqs path info
        for i in  */
        do
            cd $i/outs
            cat input_samplesheet.csv | cut -d "," -f3 | sed 1,2d | awk '{print "'$i'" "outs/fastq_path/ :" $0}' >> ../../.paths.info
            cd ../../
        done
        #----------------------------------------------------------------------------------------
        # cellranger count
        readarray sampleIDs < idoi
        for i in "${sampleIDs[@]}"
        do
            ID=$i
            GREP_PATH=$(grep `echo $i`$ .paths.info | cut -d ":" -f1  | xargs | sed -e 's/ /,/g')
            GREP_NAME=$(grep `echo $i`$ .names | xargs | sed -e 's/ /,/g')


            echo "SAMPLE_ID = $ID"
            echo "FASTQ_PATH = $GREP_PATH"
            echo "SAMPLE_NAME = $GREP_NAME"
            echo "-----------------------------------"
            echo ""

            /programs/cellranger-atac-1.2.0/cellranger-atac count --id=$ID \
                   --reference=/workdir/singleCellData/10x_reference_files/refdata-cellranger-atac-GRCh38-1.2.0 \
                   --fastqs=$GREP_PATH \
                   --sample=$GREP_NAME \
                   --localcores 20 --localmem 250

        done

        DATE=`date +"%m_%d_%H-%M"`
        mkdir CellRangerATAC-Output_${DATE}
    readarray sampleIDs < idoi
      for i in "${sampleIDs[@]}"
      do
        mv $i CellRangerATAC-Output_${DATE}
      done


else [[ "$1" = "idoi" ]]


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
        #----------------------------------------------------------------------------------------
        # get fastqs path info
        for i in  */
        do
            cd $i/outs
            cat input_samplesheet.csv | cut -d "," -f3 | sed 1,2d | awk '{print "'$i'" "outs/fastq_path/ :" $0}' >> ../../.paths.info
            cd ../../
        done
        #----------------------------------------------------------------------------------------
        # cellranger count
        readarray sampleIDs < idoi
        for i in "${sampleIDs[@]}"
        do
            ID=$i
            GREP_PATH=$(grep `echo $i`$ .paths.info | cut -d ":" -f1  | xargs | sed -e 's/ /,/g')
            GREP_NAME=$(grep `echo $i`$ .names | xargs | sed -e 's/ /,/g')


            echo "SAMPLE_ID = $ID"
            echo "FASTQ_PATH = $GREP_PATH"
            echo "SAMPLE_NAME = $GREP_NAME"
            echo "-----------------------------------"
            echo ""

            /programs/cellranger-atac-1.2.0/cellranger-atac count --id=$ID \
                   --reference=/workdir/singleCellData/10x_reference_files/refdata-cellranger-atac-GRCh38-1.2.0 \
                   --fastqs=$GREP_PATH \
                   --sample=$GREP_NAME \
                   --localcores 20 --localmem 250


        done

        DATE=`date +"%m_%d_%H-%M"`
        mkdir CellRangerATAC-Output_${DATE}
    readarray sampleIDs < idoi
      for i in "${sampleIDs[@]}"
      do
        mv $i CellRangerATAC-Output_${DATE}
      done

fi
