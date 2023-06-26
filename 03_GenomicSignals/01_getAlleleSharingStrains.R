source("~/BrusselSprouts/scripts/functions.R")
library(tidytable)

args=getArgs()
f=args[1]
df=fread(f,data.table=TRUE)
vcf=fread(cmd=paste0("zcat ",paste0(f,".loci.full.vcf.gz")," | grep -v ^##"), data.table = FALSE)
setDF(vcf)

vcf$chrpos=paste(vcf$`#CHROM`,vcf$POS)
df$chrpos=paste(df$chr,df$pos)

vcf=subset(vcf, chrpos %in% df$chrpos)

getA1A2sharersPerLocus=function(x){
    if(vcf$chrpos[x]%in%df$chrpos){
        which_df_row=which(df$chrpos==vcf$chrpos[x])
        vcf_alleles=strsplit(paste(vcf$REF[x],vcf$`ALT...5`[x],sep=","),",",fixed=TRUE)[[1]]
        df_alleles=strsplit(df$alleles[which_df_row],",",fixed=TRUE)[[1]]
        A1=which(vcf_alleles==df_alleles[as.numeric(df$Allele1[which_df_row])+1])-1 # 0 = REF
        A2=which(vcf_alleles==df_alleles[as.numeric(df$Allele2[which_df_row])+1])-1
        t_vcf=t(vcf[x,10:ncol(vcf)])
        A1_homs=which(substr(t_vcf, 1,1)==A1 & substr(t_vcf, 3,3) == A1)
        A2_homs=which(substr(t_vcf, 1,1)==A2 & substr(t_vcf, 3,3) == A2)
        vcf_return=vcf[x,c(1,2,10:ncol(vcf))]
        vcf_return[1,3:ncol(vcf_return)]="Other"
        vcf_return[1,2+A1_homs]="A1_hom"
        vcf_return[1,2+A2_homs]="A2_hom"
        return(vcf_return)
        #return(tibble(Strain=c(colnames(vcf)[A1_homs+9],colnames(vcf)[A2_homs+9]),
        #              Type=c(rep("A1_hom",length(A1_homs)),
        #                   rep("A2_hom",length(A2_homs)))))
    }else{
        return(NULL)
    }
}

summarise_strains=function(x){
    # x = strain
    return(summarise(group_by(vcf_translated, all_of(x)), p=n()/nrow(vcf_translated))%>%mutate(Strain=x)%>%rename(Type=all_of(x))%>%select(Strain,Type,p))
}

wrapperDfLocus=function(x){
    y=which(vcf$chrpos==df$chrpos[x])
    if(length(y)==0){
        return(NULL)
    }
    return(getA1A2sharersPerLocus(y))
}

vcf_translated=setDT(data.table::rbindlist(lapply(1:nrow(df),wrapperDfLocus)))
colnames(vcf_translated)[colnames(vcf_translated)=="ALT...204"]="ALT"


fwrite(vcf_translated, paste0(f,".allelesharing.csv.gz"))