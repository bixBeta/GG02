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

  for i in */

  do
    cd $i
    echo `pwd` >> ../paths
    cd ..
  done

  find . -type d | sed 1,2d | cut -d / -f2 | sort >> names

  # cellranger count
  readarray sampleIDs < idoi
  for i in "${sampleIDs[@]}"
  do
      ID=$i
      GREP_PATH=$(grep `echo $i` paths | xargs | sed -e 's/ /,/g')
      GREP_NAME=$(grep `echo $i` names | xargs | sed -e 's/ /,/g')


      echo "SAMPLE_ID = $ID"
      echo "FASTQ_PATH = $GREP_PATH"
      echo "SAMPLE_NAME = $GREP_NAME"
      echo "-----------------------------------"
      echo ""

      # /programs/cellranger-atac-1.2.0/cellranger-atac count --id=$ID \
      #        --reference=/workdir/singleCellData/10x_reference_files/refdata-cellranger-atac-GRCh38-1.2.0 \
      #        --fastqs=$GREP_PATH \
      #        --sample=$GREP_NAME \
      #        --localcores 20 --localmem 250


  done

  # DATE=`date +"%m_%d_%H-%M"`
  # mkdir CellRangerATAC-Output_${DATE}
  # readarray sampleIDs < idoi
  # for i in "${sampleIDs[@]}"
  # do
  # mv $i CellRangerATAC-Output_${DATE}
  # done
  #
  # mv .paths.info .ids .names idoi CellRangerATAC-Output_${DATE}
  # chmod -R 777 CellRangerATAC-Output_${DATE}


fi
