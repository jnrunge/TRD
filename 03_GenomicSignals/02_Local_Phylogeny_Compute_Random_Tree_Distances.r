library(ape)
library(stringr)
library(parallel)

num_cores <- 10

local_phylogeny_output_files_prefix<-"/home/jnrunge/data/trd/local_phylogenies_random_data/local_phylogenies."
local_phylogeny_output_files_postfix<-".RDS"

tr=read.tree("/home/jnrunge/data/trd/diss_global.newick")

outfile<-paste0(dirname(local_phylogeny_output_files_prefix),
                              "/tree_distances.RDS")

# check if some tree distances are up to date in their computation

if(file.exists(outfile)){
    if(file.size(outfile)>0){
        out_lock<-flock::lock(outfile)
        tree_distances_prev<-readRDS(outfile)
        flock::unlock(out_lock)
    }
}

done_files<-list.files(path=dirname(local_phylogeny_output_files_prefix),pattern=paste0(basename(local_phylogeny_output_files_prefix),
                                                                           "[0-9]*",
                                                                           local_phylogeny_output_files_postfix),full.names = TRUE)
done_files<-done_files[file.size(done_files)>0]

if(exists("tree_distances_prev")){
    done_files_IDs<-str_replace_all(done_files,
                                   local_phylogeny_output_files_prefix,
                                   "")
    
    done_files_IDs<-str_replace_all(done_files_IDs,
                                   local_phylogeny_output_files_postfix,
                                   "")
    
    done_files<-done_files[!(done_files_IDs %in% names(tree_distances_prev)) | 
                          file.mtime(done_files)>file.mtime(outfile)]
}

print(done_files)

local_phylogeny_data<-readRDS(done_files[1])
names(local_phylogeny_data)
for(f in done_files[-1]){
    tmp<-readRDS(f)
    
    for(n in names(tmp)){
        local_phylogeny_data[[n]][names(tmp[[n]])]<-tmp[[n]][names(tmp[[n]])]
    }
}

#tree_distances <- sapply(local_phylogeny_data[["rng_trees"]], function(tree) dist.topo(tr, tree, method = "score"))
tree_distances <- mclapply(local_phylogeny_data[["rng_trees"]], 
                           function(tree) dist.topo(tr, tree, method = "score"), 
                           mc.cores = num_cores)
                         
tree_distances<-c(tree_distances_prev,
                  tree_distances)
                         
out_lock<-flock::lock(outfile)
                         
saveRDS(tree_distances, outfile)
                         
flock::unlock(out_lock)