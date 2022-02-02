for i in *tagDir

do
	makeUCSCfile $i -strand separate -o auto -fsize 1e10 -res 1 -color 106,42,73 -style tss
done

