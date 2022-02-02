for i in *tagDir

do

  

  cd $i 

    for j in *.txt

    do

      jSUB=` echo $j | cut -d . -f3- `
      
      mv $j $jSUB

    done

  cd ..



done
