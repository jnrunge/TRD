# here I should just input from,to , and whether that is a square or from here until end
library(data.table)
args=c("~/OneDrive - unistra.fr/TRD/PhenotypingDiploidy/20221205/gitter/Cross_18.JPG.dat",
       "17,6", # col, row, from
       "24,16", # col, row, to
       "range") # range or region or pos (just 'from' field)

file=args[1]

from=args[2]

to=args[3]

type=args[4]

df=fread(file, sep="\t")

if(!file.exists(paste(file, ".bak",sep=""))){
  file.copy(file, paste(file, ".bak",sep=""))
}

if(type=="range"){
  df$size[((df$col==as.numeric(strsplit(from, ",")[[1]][1]) & 
            df$`# row`>=as.numeric(strsplit(from, ",")[[1]][2]))) | 
            ((df$col>as.numeric(strsplit(from, ",")[[1]][1])) & df$col<as.numeric(strsplit(to, ",")[[1]][1])) | 
            (df$col==as.numeric(strsplit(to, ",")[[1]][1]) & 
               df$`# row`<=as.numeric(strsplit(to, ",")[[1]][2]))]=9999
}

if(type=="region"){
  df$size[(df$col>=as.numeric(strsplit(from, ",")[[1]][1]) & 
             df$`# row`>=as.numeric(strsplit(from, ",")[[1]][2])) &  
            (df$col<=as.numeric(strsplit(to, ",")[[1]][1]) & 
               df$`# row`<=as.numeric(strsplit(to, ",")[[1]][2]))]=9999
}

if(type=="pos"){
  df$size[(df$col==as.numeric(strsplit(from, ",")[[1]][1]) & 
             df$`# row`==as.numeric(strsplit(from, ",")[[1]][2]))]=9999
}

writeLines(text = readLines(paste(file, ".bak",sep=""))[which(startsWith(readLines(paste(file, ".bak",sep="")), "#"))], con=file)
fwrite(df, file, sep = "\t",append=TRUE,quote = FALSE)


