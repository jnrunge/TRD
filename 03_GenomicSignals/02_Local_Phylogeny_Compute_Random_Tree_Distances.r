library(ape)

local_phylogeny_output_files_prefix<-"/home/jnrunge/data/trd/local_phylogenies_random_data/local_phylogenies."
local_phylogeny_output_files_postfix<-".RDS"

tr=read.tree("/home/jnrunge/data/trd/diss_global.newick")

done_files<-list.files(path=dirname(local_phylogeny_output_files_prefix),pattern=paste0(basename(local_phylogeny_output_files_prefix),
                                                                           "[0-9]*",
                                                                           local_phylogeny_output_files_postfix),full.names = TRUE)
done_files<-done_files[file.size(done_files)>0]

local_phylogeny_data<-readRDS(done_files[1])
names(local_phylogeny_data)
for(f in done_files[-1]){
    tmp<-readRDS(f)
    
    for(n in names(tmp)){
        local_phylogeny_data[[n]][names(tmp[[n]])]<-tmp[[n]][names(tmp[[n]])]
    }
}

tree_distances <- sapply(local_phylogeny_data[["rng_trees"]], function(tree) dist.topo(tr, tree, method = "score"))
                         
outfile<-paste0(dirname(local_phylogeny_output_files_prefix),
                              "/tree_distances.RDS")
                         
out_lock<-flock::lock(outfile)
                         
saveRDS(tree_distances, outfile)
                         
flock::unlock(out_lock)