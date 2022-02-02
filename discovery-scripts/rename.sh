for i in *tagDir

do

  iSUB=` basename $i .tagDir`

  cd $i 

    for j in *.txt

    do

      mv $j $iSUB.$j

    done

  cd ..



done
