if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<PIN>"
    echo "    where 1 is for first strand and 2 is for reverse strands, 0 for unstranded counts"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1


else 
    
    echo *.ReadsPerGene.out.tab.rawCounts | sed -e 's/ /    /g' > colnames
    paste *.ReadsPerGene.out.tab.rawCounts > $1.rawCounts.txt

    cat colnames $1.rawCounts.txt

fi

