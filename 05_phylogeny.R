#' Script to obtain the phylogenies for fish species from fishtree database
#' Authors: Oscar P. Zapata and Matheus T. Baumgartner

rm(list = ls())

# Load trait data
traits <- read.csv("portfolio.csv", sep = ",")

# ---- Phylogeny ----

library(fishtree)
library(phytools)

# Get phylogenies
phylo <- fishtree_phylogeny(species = traits$Species)

# Missing species
spp <- gsub('[ ]', '_', traits$Species)
miss_spp <- setdiff(spp, phylo$tip.label)
miss_spp[miss_spp == 'Ericymba_buccata'] <- 'Notropis_buccatus' # update species name
spp[spp == 'Ericymba_buccata'] = 'Notropis_buccatus'
rownames(traits) <- spp

# Insert missing species as polytomies in the phylogeny
phylo_tpl <- chronopl(phylo, lambda = .6, 
                      age.max = max(phylo$edge.length)) # generate an ultra-metric tree

for(i in 1:length(miss_spp)){
  phylo_tpl <- add.species.to.genus(phylo_tpl, miss_spp[i], where = 'root')
}

# Export phylogeny
writeNexus(phylo_tpl, "phylogeny.nexus")

