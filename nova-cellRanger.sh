#!/bin/sh

# initial commit FA200323
# version 1

if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<idoi> "
    # echo " 	   	auto - autodetects duplicate sampleIDs "
    # echo "		no-dups - autodetects all sampleIDs"
    echo "		idoi - expects user to provide a sample ID list named idoi"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

else [[ "$1" = "idoi" ]]

#   for i in */

#   do
#     cd $i
#     echo `pwd` >> ../paths
#     cd ..
#   done

  find . -type d | sed 1d | cut -d / -f2 | sort >> names

  # cellranger count
  cat idoi | while read i 
    do

    ID=$i
    G="$ID$|^R$ID|^$ID"

    GREP_NAME=`egrep -h $G names | xargs | sed -e 's/ /,/g'`

    echo "SAMPLE_ID = $ID"
    echo "FASTQ_PATH = $GREP_NAME"
    #    echo "SAMPLE_NAME = $GREP_NAME"
    echo "-----------------------------------"
    echo ""
    done


  
       /programs/cellranger-3.0.2/cellranger count --id=$ID \
            --transcriptome=/workdir/singleCellData/10x_reference_files/refdata-cellranger-GRCh38-3.0.0/ \
            --fastqs=$GREP_NAME \
            --sample=$GREP_NAME \
            --localcores 20 --localmem 250

    done
  
  


  DATE=`date +"%m_%d_%H-%M"`
  mkdir CellRanger-scRNAseq-Output_${DATE}
  readarray sampleIDs < idoi
  for i in "${sampleIDs[@]}"
  do
  mv $i CellRanger-scRNAseq-Output_${DATE}
  done

  mv paths names idoi CellRanger-scRNAseq-Output_${DATE}
  chmod -R 777 CellRanger-scRNAseq-Output_${DATE}

fi

