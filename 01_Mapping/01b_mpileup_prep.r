library(tidytable)
library(lubridate)
library(stringr)


samples=readLines("~/data/trd/mapped_reads/bamlist")

get_most_common_char <- function(string) {
  # Keep only A, a, C, c, T, t, G, g in the string
 # additional info is removed, because it overlaps with base letters...
  string <- gsub("[^AaCcTtGg]", "", string)
  if(nchar(string)==0){
      return(".")
  }
  # Convert the string to uppercase
  string <- toupper(string)
  
  # Get the most common character
  mode <- names(which.max(table(strsplit(string, "")[[1]])))
  
  # Return the most common character
  return(mode)
}

getAF=function(x, alt){
        # is content of one sample's cell
        ref=str_count(x, "[,.]")
        alt=str_count(x, paste0("[",tolower(alt),alt,"]"))
        return(alt/(alt+ref))
    }
    library(tidyverse)


source("/home/jnrunge/BrusselSprouts/scripts/functions.R")

args=getArgs()

batch_size <- as.numeric(args[3])
i_min=as.numeric(args[1]) # start
i_max=as.numeric(args[2]) # stop

# Open the gzipped file connection
file_conn <- gzfile("~/data/trd/mapped_reads/all_mpileup.txt.gz",open="r")

dt <- data.table()

i=0


if(i_min>0){
    scan(file_conn, nlines = ((i_min-1)*batch_size), what = character())
    i=i_min
}


# Read in the file in batches of lines
while(length(batch <- readLines(file_conn, n = batch_size)) > 0) {
    i=i+1
    if(i>i_max){
        break
    }
    all_mpileup = fread(text=batch, header=FALSE,sep="\t", drop=sort(c(4+(3*0:(length(samples)-1)), 6+(3*0:(length(samples)-1)))))
    
    # do all the modifications and filtering
    
    
    colnames(all_mpileup)[4:ncol(all_mpileup)]=samples
    
    all_mpileup=unite(all_mpileup, allpaste, all_of(samples), sep = "", remove=FALSE)

    all_mpileup=mutate(all_mpileup, allpaste = gsub("\\^\\w[,|.]","", allpaste))
    all_mpileup=filter(all_mpileup, grepl("[AaCcTtGg]", allpaste))

    all_mpileup <- all_mpileup %>% 
      rowwise() %>% 
      mutate(ALT = get_most_common_char(allpaste))


    all_mpileup=filter(all_mpileup, ALT!=".")

    #all_mpileup_bak=all_mpileup

    #all_mpileup=all_mpileup_bak
    
    all_mpileup <- all_mpileup %>% 
      mutate_at(vars(ends_with(".bam")), ~as.numeric(getAF(., ALT)))    
    
    all_mpileup=select(all_mpileup, -allpaste)
    
  # Combine the batch data.table with the previous data.table
  dt <- bind_rows(list(dt, all_mpileup))

}
dt <- dt %>%
  mutate_at(4:(ncol(.)-1), round, digits = 3)
fwrite(dt, paste0("~/data/trd/mapped_reads/all_mpileup.txt.gz", "-", i_min, "-", i_max, ".gz"), sep=",", col.names=TRUE)

# Close the file connection
close(file_conn)
