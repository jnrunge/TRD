#!~/anaconda3/envs/JupyteR4/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

source("~/BrusselSprouts/scripts/functions.R")


example=fread(cmd=paste("zcat /home/jnrunge/data/TRD/AD/header.AD.csv.gz ",args[1],sep=""), data.table=FALSE, fill=TRUE)
for(i in 1:ncol(example)){
    example[,i]=as.numeric(example[,i])
}
getAF=function(x)
    {
    fields=sum(!is.na(t(example[x,])))
    depth=sum(t(example[x,][1:fields]))
    AF_df=data.frame(AF1=example[x,1]/depth,AF2=example[x,2]/depth,AF3=example[x,3]/depth,AF4=example[x,4]/depth)
    return(AF_df)
}
AF_df=bind_rows(lapply(1:nrow(example),getAF))

polyploidy_signal_fraction=sum(((AF_df$AF1> 0.1 & AF_df$AF1 < 0.4) | (AF_df$AF1< 0.9 & AF_df$AF1>0.6)) | 
    ((AF_df$AF2> 0.1 & AF_df$AF2 < 0.4) | (AF_df$AF2< 0.9 & AF_df$AF2>0.6)) | 
    ((AF_df$AF3> 0.1 & AF_df$AF3 < 0.4) | (AF_df$AF3< 0.9 & AF_df$AF3>0.6)) |
    ((AF_df$AF4> 0.1 & AF_df$AF4 < 0.4) | (AF_df$AF4< 0.9 & AF_df$AF4>0.6)),na.rm = TRUE)/nrow(AF_df)

writeLines(as.character(polyploidy_signal_fraction),args[2])