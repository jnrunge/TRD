# inputs:
# location of data frame that gives the chr, pos_from, pos_to
# i; i.e. the row that we are computing in this instance

source("~/BrusselSprouts/scripts/functions.R")
library(tidytable)
library(SNPRelate)
library(ape)

euclidean_distance <- function(x, y) {
  sqrt(x^2 + y^2)
}

local_phylogeny_output_file<-"/home/jnrunge/data/trd/local_phylogenies.RDS"

# Get paths of input files
my.vcf <- "/home/jnrunge/data/trd/full2489Matrix.vcf.gz"
subset_samples = "none"

if(subset_samples == "none"){
    subset_samples=NULL
    prefix <- basename(my.vcf)
    }else{
    prefix <- basename(subset_samples)
    subset_samples=readLines(subset_samples)
    }

setwd(dirname(my.vcf))

genofile <- snpgdsOpen(paste0(prefix, ".gds"))


args=getArgs()

rng<-fread(args[1])
i<-as.numeric(args[2])
rng_row<-slice(rng,i)

snpset=which(read.gdsn(index.gdsn(genofile, "snp.position"))>=pull(rng_row,from)&
                 read.gdsn(index.gdsn(genofile, "snp.position"))<=pull(rng_row,to)&
                 read.gdsn(index.gdsn(genofile, "snp.chromosome"))==pull(rng_row,chr))

print(head(snpset))

pca_rng <- snpgdsPCA(genofile, num.thread=1, snp.id = snpset)

    x_pca<-pca_rng$eigenvect[,1]
    y_pca<-pca_rng$eigenvect[,2]

ibs <- snpgdsIBS(genofile, num.thread=1, snp.id = snpset)
    loc <- cmdscale(1 - ibs$ibs, k = 2)

    x_ibs <- loc[, 1]
    y_ibs <- loc[, 2]

dissMatrix_rng <- snpgdsDiss(genofile, snp.id=snpset, num.thread = 1)
    colnames(dissMatrix_rng$diss) <- dissMatrix_rng$sample.id
    tr <- bionjs(dissMatrix_rng$diss)

    tr_dist<-ape::cophenetic.phylo(tr)

print("Locking output file.")
local_phylogeny_output_file_lock<-flock::lock(local_phylogeny_output_file)

if(file.size(local_phylogeny_output_file) == 0){
    PCA_distances_from_0=list()
    IBS_MDS_distances_from_0=list()
    rng_trees<-list()
    rng_tree_distances<-list()
    
    
    
    local_phylogeny_data<-list()
    local_phylogeny_data[["PCA_distances_from_0"]]<-PCA_distances_from_0
    local_phylogeny_data[["IBS_MDS_distances_from_0"]]<-IBS_MDS_distances_from_0
    local_phylogeny_data[["rng_trees"]]<-rng_trees
    local_phylogeny_data[["rng_tree_distances"]]<-rng_tree_distances
    

}else{
    local_phylogeny_data<-readRDS(local_phylogeny_output_file)
    
}

local_phylogeny_data[["PCA_distances_from_0"]][[as.character(i)]]<-euclidean_distance(x_pca,y_pca)
local_phylogeny_data[["IBS_MDS_distances_from_0"]][[as.character(i)]]<-euclidean_distance(x_ibs,y_ibs)
local_phylogeny_data[["rng_trees"]][[as.character(i)]]<-tr
local_phylogeny_data[["rng_tree_distances"]][[as.character(i)]]<-tr_dist

saveRDS(local_phylogeny_data, file=local_phylogeny_output_file)

flock::unlock(local_phylogeny_output_file_lock)
print("Unlocking output file.")
