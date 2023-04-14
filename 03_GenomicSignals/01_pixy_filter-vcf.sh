cd ~/data/trd/mapped_reads

# arg 1 = chromosome

bcftools view -r $1 ALL.DP5-95.bcf | bcftools +setGT  -- -t q -n . -e 'FMT/DP>=10' | bcftools +setGT -- -t q -n . -e 'FMT/GQ>=20 | FMT/RGQ>=20' | bcftools +fill-tags | bcftools view -e 'ALT="*" || (type!="snp" && type!="ref")' | bgzip > ALL.DP5-95.$1.DP10.GQ20RGQ20.SNPsRef.vcf.gz

tabix ALL.DP5-95.$1.DP10.GQ20RGQ20.SNPsRef.vcf.gz