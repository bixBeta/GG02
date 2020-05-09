#!/bin/sh

if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash $0 <delim,field> "
    echo "   comma spearated values for delimiter and field"
    echo "    simple example: -,2 "
    echo "    complex example: \"-,2 | tail -c 4\" or \"-,2 | grep -o ...$\""
    echo "--------------------------------------------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

  else

  DELIM=$1
  for i in *
  do

    DELIMITER=`echo $DELIM | cut -d , -f1`
    FIELD=`echo $DELIM | cut -d , -f2- | cut -d "|" -f1`
    CCOUNT=`echo $DELIM | cut -d , -f2- | cut -d "|" -f2-`

        if  echo $DELIM | grep -q "|"

        then
        echo $i | cut -d ${DELIMITER} -f${FIELD} | ${CCOUNT}

        else
        echo $i | cut -d ${DELIMITER} -f${FIELD}

        fi

  done

fi
