source("~/BrusselSprouts/scripts/functions.R")

args=getArgs()
cross=args[1]
sample1=args[2]
sample2=args[3]


loci1=fread(paste0("~/data/TRD/Oppo-Homo-Pos/",sample1,".loci"))
loci2=fread(paste0("~/data/TRD/Oppo-Homo-Pos/",sample2,".loci"))

loci=subset(loci1, paste(V1,V2) %in% paste(loci2$V1,loci2$V2))

fwrite(loci, paste0("~/data/TRD/Oppo-Homo-Pos/",cross,".ohloci"), sep="\t", col.names=FALSE)