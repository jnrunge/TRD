source("/home/jnrunge/BrusselSprouts/scripts/functions.R")

args=getArgs()
f=args[1]

MQ=fread(f, data.table=FALSE, header=FALSE)

MQ=subset(MQ, V3!=0)

for(c in unique(MQ$V1)){
    tmp=data.frame(chr=c,
                   pos=seq(from=min(MQ$V2[MQ$V1==c]), to=max(MQ$V2[MQ$V1==c]), by=10000))
    
    if(c==MQ$V1[1]){
        ranges=tmp
    }else{
        ranges=bind_rows(ranges, tmp)
    }
}

getMeanMQ=function(x){
    tmp=subset(MQ, V1 == ranges$chr[x] & V2 <= ranges$pos[x] & V2>ranges$pos[x-1])
    n=nrow(tmp)
    return(data.frame(meanMQ=mean(tmp$V3), n=n))
}

df_MQ=bind_cols(ranges[-1,],bind_rows(lapply(2:nrow(ranges), getMeanMQ)))

fwrite(df_MQ, paste(f,"-windows.tsv.gz",sep=""))