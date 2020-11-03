if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<PIN>"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1


else

    echo *.ReadsPerGene.out.tab.rawCounts | sed -e 's/ /    /g' > colnames
    paste *.ReadsPerGene.out.tab.rawCounts > $1.rawCounts.txt.temp

<<<<<<< HEAD
    cat colnames $1.rawCounts.txt > $1.rawCounts.txt
=======
    cat colnames $1.rawCounts.txt.temp > $1.rawCounts.txt
>>>>>>> 74060fac230d814614cbe5dfa4ef795acf9da549

fi
