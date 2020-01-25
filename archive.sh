#!/bin/sh

if [ "$1" = "help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0  "</path/to/fastqs/> <PIN>  "
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

else

PIN=$2

mkdir /home/RSCshare/RSC/Projects/ARCHIVE/$PIN

rsync -av $1 /home/RSCshare/RSC/Projects/ARCHIVE/$PIN

fi

