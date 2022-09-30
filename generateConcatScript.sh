#!/bin/bash

if [ "$1" = "help" ] || [ -z "$1" ] || [  $1 = "-h" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0
    echo "    Generates cat.sh script "
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1


else

    ls -1 *gz > allSamps

    cut -d _ -f4 allSamps | sort -u | sed 'H;1h;$!d;x;y/\n/,/' > flowID

    A=`cut -d , -f1 flowID`
    B=`cut -d , -f2 flowID`

    ls -1 *$A*gz > .As
    ls -1 *$B*gz > .Bs

    paste .As .Bs > .temp

    cat .temp |while read i; do echo cat $i '>' c`echo $i | cut -d " " -f2`; done >> cat.sh

fi
