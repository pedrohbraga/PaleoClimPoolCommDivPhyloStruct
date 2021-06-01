#######################################################################################################

#### Load phylogenetic tree ####

Mammal.FaurSven.tree.1 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_1.nex")
Mammal.FaurSven.tree.2 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_2.nex")
Mammal.FaurSven.tree.3 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_3.nex")

Mammal.FaurSven.trees <- c(Mammal.FaurSven.tree.1, Mammal.FaurSven.tree.2, Mammal.FaurSven.tree.3)

Mammal.FaurSven.MCC.tree <- phangorn::maxCladeCred(Mammal.FaurSven.trees)

#Chiroptera.tree

#Mammal.FaurSven.tree

# Filter species from the phylogenetic tree
#ChiropteraPhylo <- drop.tip(Chiroptera.tree, 
#                               Chiroptera.tree$tip.label[-match(colnames(Chiroptera.Comm), 
#                                                                Chiroptera.tree$tip.label)]) #Prunning the tree
#length(ChiropteraPhylo$tip.label); ncol(Chiroptera.Comm)

Chiroptera.FaurSven.tree <- match.phylo.comm(Mammal.FaurSven.MCC.tree, Chiroptera.Comm)$phy
Chiroptera.FaurSven.comm <- match.phylo.comm(Mammal.FaurSven.MCC.tree, Chiroptera.Comm)$comm

# ChiropteraPhylo
is.ultrametric(Chiroptera.FaurSven.tree)

# Force tree to be ultrametric
# Chiroptera.FaurSven.Phylo.ultra <- force.ultrametric(Chiroptera.FaurSven.tree)

# Check if number of species matches between both trait and phylogenetic relationship datasets
length(Chiroptera.FaurSven.tree$tip.label); ncol(Chiroptera.FaurSven.comm)

write.tree(Chiroptera.FaurSven.tree, file = "data/phylogenies/Chiroptera.FaurSven.tree.phy")

###################

library(tidyverse)
library(ggtree)

# or add the entire scale to the x axis with theme_tree2()
ggtree(Chiroptera.FaurSven.tree) + theme_tree2()

# Label the tips
ggtree(Chiroptera.FaurSven.tree) + theme_tree2() + geom_tiplab()


setdiff(colnames(Chiroptera.Comm), Mammal.FaurSven.MCC.tree$tip.label)
