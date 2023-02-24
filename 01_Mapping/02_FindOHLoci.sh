# get loci

touch ~/data/TRD/Oppo-Homo-Pos/$3.running

bcftools view -V indels,mnps,other -s $1 ~/data/trd/full2489Matrix.vcf.gz | bcftools view -i 'F_MISSING==0 & GT[0] == "hom" & GQ >=30 && FMT/DP >= 20' | bcftools view -e 'type!="snp"' | grep -v ^# | cut -f 1,2 > ~/data/TRD/Oppo-Homo-Pos/$1.loci

bcftools view -V indels,mnps,other -s $2 ~/data/trd/full2489Matrix.vcf.gz | bcftools view -i 'F_MISSING==0 & GT[0] == "hom" & GQ >=30 && FMT/DP >= 20' | bcftools view -e 'type!="snp"' | grep -v ^# | cut -f 1,2 > ~/data/TRD/Oppo-Homo-Pos/$2.loci

. ~/activate.sh fixR

Rscript ~/TRD/01_Mapping/02_combineLoci.r $3 $1 $2

. ~/activate.sh bwaetc

gzip -f ~/data/TRD/Oppo-Homo-Pos/$3.ohloci
gzip -f ~/data/TRD/Oppo-Homo-Pos/$1.loci
gzip -f ~/data/TRD/Oppo-Homo-Pos/$2.loci


# get genotype owners

bcftools view -T ~/data/TRD/Oppo-Homo-Pos/$3.ohloci.gz -s $1,$2 ~/data/trd/full2489Matrix.vcf.gz | grep -v ^# | cut -f 1,2,4,5,9,10,11 | gzip > ~/data/TRD/Oppo-Homo-Pos/$3.ohloci.GT.gz


touch ~/data/TRD/Oppo-Homo-Pos/$3.done