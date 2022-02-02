$(awk 'BEGIN {OFS="\t"; m=0} \
FNR==NR {d[$1]=1; next} { if(FNR%2==1) {s=substr($1,2,3); \
c=substr($1,match($1,"_x")+2,15); \
if(substr($1,2,25) in d) m=1} \
else { l=length($1); f=substr($1,1,1); \
a[s,l,f,m]++; b[s,l,f,m]=b[s,l,f,m]+c;s=0;c=0;l=0;f=0; m=0}} \
END {print "library\treadlength\tbase1\tmiRBaseMatch\t#distinctReads\t#reads"; \
for (var in a) {split(var,q,SUBSEP); print q[1], q[2], q[3], q[4], a[var], b[var]} }' <(zcat *miRBase.mrd.gz) <(zcat *collapsed.fa.gz) > mirmap_firstbase_readlengthcounts.txt)
