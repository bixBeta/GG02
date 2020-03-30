#!/bin/sh

if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<delimiter> <field>"
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

  else

  for i in *.gz
  do

    echo $i | `echo cut -d ${1} -f${2}`
  done

fi
