######################################################################################
### Code to prepare the data on bat communities and bat phylogenetic relationships   #
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2022-01-12"                                                          #
#                                                                                    #
######################################################################################

#### Preparing the dataset of bat communities  ---------------------------------------

# Load community data ###

Chiroptera.Comm <- read.csv("data/matrices/world_chiroptera_PropCover_50KM.csv",
                            header = TRUE)

nrow(Chiroptera.Comm); nrow(world_grid_50km_cat_df) # Verify lengths

# Replace periods in species names by underscores

colnames(Chiroptera.Comm) <- gsub('\\.', 
                                  "_", 
                                  colnames(Chiroptera.Comm))

head(Chiroptera.Comm); nrow(Chiroptera.Comm)

# Remove species that do not occur in any community
# To identify them: colnames(Chiroptera.Comm[ , colSums(Chiroptera.Comm) < 1])

# Remove species (i.e. columns) that have sum zero
Chiroptera.Comm <- Chiroptera.Comm[ , colSums(Chiroptera.Comm) >= 1]

Chiroptera.Comm[rowSums(Chiroptera.Comm) < 1, ]

# Drop empty-case levels

world_grid_50km_cat_df$ID_Realm <- droplevels(world_grid_50km_cat_df$ID_Realm)
world_grid_50km_cat_df$ID_Biome <- droplevels(world_grid_50km_cat_df$ID_Biome)

### Preparing the phylogenetic data  -----------------------------------------

# Load phylogenetic trees from Faurby and Svenning, 2015 (doi:10.1016/j.ympev.2014.11.00)
# They are available in Appendix D, within the Supplementary Material

Mammal.FaurSven.tree.1 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_1.nex")
Mammal.FaurSven.tree.2 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_2.nex")
Mammal.FaurSven.tree.3 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_3.nex")

Mammal.FaurSven.trees <- c(Mammal.FaurSven.tree.1, 
                           Mammal.FaurSven.tree.2, 
                           Mammal.FaurSven.tree.3)

# Obtain the maximum credibility tree

# This tree will be used for subsequent analyses, but also see the Supplementary
# Material X, where we repeated our analyses across subsets of trees randomly obtained
# from Mammal.FaurSven.trees and show that they are very similar.

Mammal.FaurSven.MCC.tree <- phangorn::maxCladeCred(Mammal.FaurSven.trees)

# Calculate the log clade credibility scores, which will be used to sample other
# trees to assess whether the results from the maximum credibility clade tree
# differs across different trees

Mammal.FaurSven.MCC.tree.scores <- phangorn::maxCladeCred(Mammal.FaurSven.trees, 
                                                          tree = FALSE)

# Represent log credibility scores in a histogram 

Mammal.FaurSven.MCC.tree.scores.histogram <- hist(Mammal.FaurSven.MCC.tree.scores,
                                                  breaks = 25,
                                                  col = "gray",
                                                  main = "",
                                                  xlab = "log(Clade Credibility for Faurby and Svenning's (2015) tree)")

arrows(x0 = max(Mammal.FaurSven.MCC.tree.scores), 
       y0 = 0.5*max(Mammal.FaurSven.MCC.tree.scores.histogram$counts),
       x1 = max(Mammal.FaurSven.MCC.tree.scores),
       y1 = 0,
       col = "red",
       lwd = 2,
       length = 0.15, 
       angle = 20)

text(x = max(Mammal.FaurSven.MCC.tree.scores),
     y = 0.5*max(Mammal.FaurSven.MCC.tree.scores.histogram$counts)+ 1,
     "Clade credibility of the tree we used",
     adj = c(0, 0.3),
     srt = 90)

arrows(x0 = median(Mammal.FaurSven.MCC.tree.scores), 
       y0 = 0.5*median(Mammal.FaurSven.MCC.tree.scores.histogram$counts),
       x1 = median(Mammal.FaurSven.MCC.tree.scores),
       y1 = 0,
       col = "blue",
       lwd = 2,
       length = 0.15,
       angle = 20)

Mammal.FaurSven.MCC.tree.scores

# If you would like to choose the trees belonging to the 5% quantile of log-cred scores 

Mammal.FaurSven.MCC.tree.scores[Mammal.FaurSven.MCC.tree.scores < quantile(Mammal.FaurSven.MCC.tree.scores, 0.05)]

# Remove species that are unmatched between the phylogeny and the community data

# ChiropteraPhylo <- drop.tip(Chiroptera.tree, 
#                               Chiroptera.tree$tip.label[-match(colnames(Chiroptera.Comm), 
#                                                                Chiroptera.tree$tip.label)]) #Prunning the tree
# length(ChiropteraPhylo$tip.label); ncol(Chiroptera.Comm)

Chiroptera.FaurSven.tree <- picante::match.phylo.comm(Mammal.FaurSven.MCC.tree, 
                                                      Chiroptera.Comm)$phy

Chiroptera.FaurSven.comm <- picante::match.phylo.comm(Mammal.FaurSven.MCC.tree, 
                                                      Chiroptera.Comm)$comm

# Save trees

# MCC tree
write.tree(Mammal.FaurSven.MCC.tree, 
           file = "data/phylogenies/Mammal.FaurSven.MCC.tree.phy")

# Chiroptera MCC tree, filtered species

write.tree(Chiroptera.FaurSven.tree, 
           file = "data/phylogenies/Chiroptera.FaurSven.tree.phy")

# For future analyses, load them from here

# Mammal.FaurSven.MCC.tree <- read.tree("data/phylogenies/Mammal.FaurSven.MCC.tree.phy")
# Chiroptera.FaurSven.tree <- read.tree("data/phylogenies/Chiroptera.FaurSven.tree.phy")

# ChiropteraPhylo
is.ultrametric(Chiroptera.FaurSven.tree)

# Force tree to be ultrametric
# Chiroptera.FaurSven.Phylo.ultra <- force.ultrametric(Chiroptera.FaurSven.tree)

# Check if number of species matches between both trait and phylogenetic relationship
# datasets
length(Chiroptera.FaurSven.tree$tip.label); ncol(Chiroptera.FaurSven.comm)

######################################################################################

# Representation of the resulting phylogenetic tree, just for visual purposes

# library(tidyverse)
# library(ggtree)

# or add the entire scale to the x axis with theme_tree2()
# ggtree(Chiroptera.FaurSven.tree) + theme_tree2()

# Label the tips
# ggtree(Chiroptera.FaurSven.tree) + theme_tree2() + geom_tiplab()

######################################################################################

# Remove objects that will not be used later

rm(
  Mammal.FaurSven.tree.1,
  Mammal.FaurSven.tree.2,
  Mammal.FaurSven.tree.3,
  Mammal.FaurSven.MCC.tree,
  Mammal.FaurSven.trees,
  Mammal.FaurSven.MCC.tree.scores,

)
