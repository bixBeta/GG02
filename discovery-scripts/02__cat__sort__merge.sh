
echo 'concating ----------------------'

cat *bed > concatted__bed

echo 'sorting   ----------------------'

cat concatted__bed | sort -k1,1 -k2,2n > concatted__sorted__bed

echo 'merging  ----------------------'

#bedtools merge -i concatted__sorted__bed > MERGED.bed

#L=`wc -l MERGED.bed`

#seq 1 `echo $L | cut -d ' ' -f1` > .L

#paste .L MERGED.bed > saf.input.bed

bedtools merge -i concatted__sorted__bed -c 4 -o distinct > SAF1
