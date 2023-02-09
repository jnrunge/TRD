# get loci

bcftools view -V indels,mnps,other -s $1,$2 ~/data/TRD/1011Matrix.gvcf.gz | bcftools view -i 'F_MISSING==0 && GT[0] == "hom" && GT[1] == "hom" && AN != AC && AC != 0' | grep -v ^# | cut -f 1,2 > ~/data/TRD/Oppo-Homo-Pos/$3.ohloci

# get genotype owners

bcftools view -T ~/data/TRD/Oppo-Homo-Pos/$3.ohloci -s $1,$2 ~/data/TRD/1011Matrix.gvcf.gz | grep -v ^# | cut -f 1,2,4,5,10,11 | gzip > ~/data/TRD/Oppo-Homo-Pos/$3.ohloci.GT.gz

gzip -f ~/data/TRD/Oppo-Homo-Pos/$3.ohloci

touch ~/data/TRD/Oppo-Homo-Pos/$3.done