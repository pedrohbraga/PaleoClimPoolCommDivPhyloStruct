#################################################################################
# Testing the implementation of expanded grids to calculate the standardized    #
# effect size for mean phylogenetic pairwise distances and mean nearest taxon   #
# distances using PhyloMeasures and picante. This script also includes the      #
# mpd.new() (from Pedro Peres-Neto), which relies on matrix multiplication to   #
# calculate MPD.                                                                #   
#                                                                               #
# Other functions are borrowed from other scripts and will be detailed later.   # 
#
# Note: Future modifications should include the other null models, as well as   #
# improvements in speed supported by benchmarking.                              #
#                                                                               #
# Last Update: "2021-12-27"                                                     #
#                                                                               # 
#################################################################################

# Load required packages

library(phytools)
library(picante)
library(geiger)
library(microbenchmark)
library(matrixStats)
library(PhyloMeasures)
library(snowfall)

# Required functions

generate_community <- function(E,T,preset_nspecies,preset_ncommunities,param){
  repeat {
    # one trait, one environmental variable
    L <- matrix(data = 0, nrow=preset_ncommunities,ncol=preset_nspecies)
    for(j in 1:preset_nspecies){
      L[,j]<-rbinom(preset_ncommunities,1,(1/(1+exp(-param[1]*E+param[2]*T[j]+param[3]*E*T[j]))))
    }
    n_species_c<-sum(colSums(L)!=0) #
    n_communities_c<-sum(rowSums(L)!=0)
    if ((n_species_c==preset_nspecies) & (n_communities_c==preset_ncommunities)){break}
  }
  return(L)
}


mpd.new <- function(species.data,phylo.dist,n.perm){
  n.communities <- nrow(species.data)
  n.species <- ncol(species.data)
  mpd.rnd <- matrix(0,n.perm,n.communities)
  sum.communities <- diag(species.data %*% t(species.data))
  sum.communities <- sum.communities^2-sum.communities
  
  mpd.obs <- rowSums((species.data %*% phylo.dist)*species.data)/sum.communities
  # this way of multiplying is much faster than different options, e.g.:
  # diag((species.data %*% phylo.dist %*% t(species.data)))/sum.communities    
  # diag(tcrossprod(species.data %*% phylo.dist,(species.data)))/sum.communities
  perms <- replicate(n.perm, sample(n.species)) # note that a loop would be faster than replicate here
  # the use of "=" speeds up over "<-"
  for (i in 1:n.communities){
    tmp = t(matrix(species.data[i,perms],n.species,n.perm)) # this step speeds up compared to having it twice calculated below
    mpd.rnd[,i] = rowSums((tmp %*% phylo.dist)*tmp)/sum.communities[i]
  }
  mpd.rand.mean <- apply(mpd.rnd,2,mean)
  mpd.rand.sd <- apply(mpd.rnd,2,sd)
  mpd.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
  array.obs <- t(do.call(cbind, lapply(seq_len(n.perm), function(X) mpd.obs)))
  p.values <- (apply(mpd.rnd < array.obs,2, sum, na.rm = TRUE)+1)/(n.perm+1)
  matrix.result <- data.frame(ntaxa=apply(species.data,1,sum),mpd.obs,mpd.rand.mean,mpd.rand.sd,
                              mpd.z,p.values,n.perm)
  return(matrix.result)
} 

### Main functions

# Implementation using PhyloMeasures

mod.ses.mpd.query <- function(phylo.tree = phylo.tree,
                              perms = 999,
                              species.data = species.data){
  
  taxaShuffle.map <- purrr::map(1:perms, 
                                ~shuffle_tiplabels(phylo.tree))
  
  random.mpd.query <- sapply(X = taxaShuffle.map,
                             FUN = mpd.query,
                             matrix = species.data,
                             USE.NAMES = TRUE,
                             simplify = TRUE,
                             standardize = FALSE)
  
  random.mpd.query.mean <- apply(X = random.mpd.query, 
                                 MARGIN = 1, # by rows
                                 FUN = mean)
  
  random.mpd.query.sd <- apply(random.mpd.query, 
                               1, 
                               sd)
  
  mpd.obs.query <- mpd.query(phylo.tree,
                             species.data)
  
  names(mpd.obs.query) <- row.names(species.data)
  
  ses.mpd.z.query <- (mpd.obs.query - random.mpd.query.mean)/random.mpd.query.sd
  
  mpd.obs.rank <- apply(X = rbind(mpd.obs.query, random.mpd.query.mean), 
                        MARGIN = 2, 
                        FUN = rank)[1, ]
  
  mpd.obs.rank <- ifelse(is.na(random.mpd.query.mean), 
                         NA, 
                         mpd.obs.rank)
  
  return(data.frame(ntaxa = apply(species.data, 1, sum),
                    mpd.obs.query,
                    random.mpd.query.mean,
                    random.mpd.query.sd,
                    ses.mpd.z.query,
                    mpd.obs.p = mpd.obs.rank/(perms + 1),
                    perms))
}

## Implementation using picante

mod.ses.mpd.picante <- function(phylo.tree = phylo.tree,
                                perms = 999,
                                species.data = species.data){
  
  ### picante::mpd() requires cophenetic distances
  
  taxaShuffle.map <- purrr::map(1:perms, 
                                ~shuffle_tiplabels(phylo.tree))
  
  taxaShuffle.map.l <- lapply(taxaShuffle.map, cophenetic)
  
  random.mpd.picante <- sapply(X = taxaShuffle.map.l,
                               FUN = mpd,
                               samp = species.data,
                               USE.NAMES = TRUE,
                               simplify = TRUE)
  
  random.mpd.picante.mean <- apply(X = random.mpd.picante, 
                                   MARGIN = 1, # by rows
                                   FUN = mean)
  
  random.mpd.picante.sd <- apply(random.mpd.picante, 
                                 MARGIN = 1, 
                                 FUN = sd)
  
  mpd.obs.picante <- mpd(dis = cophenetic(phylo.tree),
                         samp = species.data)
  
  names(mpd.obs.picante) <- row.names(species.data)
  
  ses.mpd.z.picante <- (mpd.obs.picante - random.mpd.picante.mean)/random.mpd.picante.sd
  
  mpd.obs.rank <- apply(X = rbind(mpd.obs.picante, 
                                  random.mpd.picante.mean), 
                        MARGIN = 2, 
                        FUN = rank)[1, ]
  
  mpd.obs.rank <- ifelse(is.na(random.mpd.picante.mean), 
                         NA, 
                         mpd.obs.rank)
  
  return(data.frame(ntaxa = specnumber(species.data),
                    mpd.obs.picante,
                    random.mpd.picante.mean,
                    random.mpd.picante.sd,
                    ses.mpd.z.picante,
                    mpd.obs.p = mpd.obs.rank/(perms + 1),
                    perms)
  )
}

# Create species occurence matrices and a phyogenetic tree

# Set seed
set.seed(100)

n.species <- 1000
n.communities <- 60000

simulated.tree <- pbtree(n = n.species, 
                         scale = 1)

phylo.dist <- cophenetic(simulated.tree)

T.phylo <- scale(
  rTraitCont(
    geiger::rescale(simulated.tree, 
                    "delta", 
                    0.01)
  )
) # very conserved phylogenetically


species.data <- generate_community(E = rnorm(n.communities),
                                   T.phylo,
                                   preset_nspecies = n.species,
                                   preset_ncommunities = n.communities,
                                   c(2, 2, 2))

colnames(species.data) <- simulated.tree$tip.label

row.names(species.data) <- sprintf("comm_%02d", 
                                   seq(1, nrow(species.data))
)

# forces one (the first) community to have only one species for testing the
# function

species.data[1, ] <- c(1,
                       rep(0,
                           n.species - 1))


# Setting general parameters 

phylo.tree <- simulated.tree
perms <- 99

## Expanding permutations in a list() grid format with purrr::map()

### PhyloMeasures::mpd.query() requires a tree

shuffle_tiplabels <- function(tree){
  tree$tip.label <- sample(tree$tip.label)
  return(tree)
}

# Dummy runs


system.time({test.mod.ses.mpd.picante <- mod.ses.mpd.picante(phylo.tree = phylo.tree,
                                                             perms = 999,
                                                             species.data = species.data)})

system.time({test.mod.ses.mpd.query <-  mod.ses.mpd.query(phylo.tree = phylo.tree,
                                                          perms = 999,
                                                          species.data = species.data)})

mpd.new.ppn <- mpd.new(species.data, 
                       cophenetic(phylo.tree), 
                       n.perm = perms)

ses.mpd.picante <- ses.mpd(samp = species.data,
                           dis = cophenetic(phylo.tree),
                           null.model = "taxa.labels",
                           runs = 99)

sf.ses.mpd.picante <- sf.ses.mpd(samp = species.data,
                                 dis = cophenetic(phylo.tree),
                                 runs = 99)


#### Some cross-checking

plot(test.mod.ses.mpd.picante$ses.mpd.z.picante,
     test.mod.ses.mpd.query$ses.mpd.z.query)

plot(ses.mpd.picante$mpd.obs.z, 
     test.mod.ses.mpd.query$ses.mpd.z.query)

plot(ses.mpd.picante$mpd.obs.z, 
     ses.mpd.mpd.picante)

plot(ses.mpd.picante$mpd.obs.z, 
     mpd.new.ppn$mpd.z)

plot(ses.mpd.query, 
     mpd.new.ppn$mpd.z)

plot(sf.ses.mpd.phpb$mpd.obs.z, 
     ses.mpd.query)

# Bench press and mark comparisons
### Note: this is slower than using microbenchmark and similar packages

changing.communities.benchpress <- press(
  p = seq(100, 
          1000, 
          length.out = 10),
  {
    n.species <- 1000
    n.communities <- p
    
    simulated.tree <- pbtree(n = n.species, 
                             scale = 1)
    
    phylo.dist <- cophenetic(simulated.tree)
    
    T.phylo <- scale(rTraitCont(geiger::rescale(simulated.tree, 
                                                "delta", 
                                                0.01)
    )
    ) # very conserved phylogenetically
    
    
    species.data <- generate_community(E = rnorm(n.communities),
                                       T.phylo,
                                       preset_nspecies = n.species,
                                       preset_ncommunities = n.communities,
                                       c(2, 2, 2))
    
    colnames(species.data) <- simulated.tree$tip.label
    
    row.names(species.data) <- sprintf("comm_%02d", seq(1, nrow(species.data)))
    
    # forces one (the first) community to have only one species for testing the
    # function
    
    species.data[1, ] <- c(1,
                           rep(0,
                               n.species - 1))
    
    
    # Setting general parameters 
    
    phylo.tree <- simulated.tree
    
    perms <- 999
    n.cores <- 8
    
    bench::mark(
      min_time = Inf,
      
      test.mod.ses.mpd.picante <- mod.ses.mpd.picante(phylo.tree = phylo.tree,
                                                      perms = perms,
                                                      species.data = species.data),
      test.mod.ses.mpd.query <-  mod.ses.mpd.query(phylo.tree = phylo.tree,
                                                   perms = perms,
                                                   species.data = species.data),
      mpd.new.ppn <- mpd.new(species.data, 
                             cophenetic(phylo.tree), 
                             n.perm = perms),
      ses.mpd.picante <- ses.mpd(samp = species.data,
                                 dis = cophenetic(phylo.tree),
                                 null.model = "taxa.labels",
                                 runs = perms),
      mpd.new.parallel.ppn <- mpd.new.parallel(species.data, 
                                               cophenetic(phylo.tree),
                                               perms, 
                                               n.cores = n.cores),
      test.sf.ses.mpd <- sf.ses.mpd(species.data, 
                                    cophenetic(phylo.tree), 
                                    runs = perms, 
                                    cores = n.cores),
      max_iterations = 1,
      min_iterations = 1,
      iterations = 1,
      check = FALSE,
      memory = FALSE
    )
  }
)

# Representing the results from the comparsions

library(ggplot2)
library(tidyr)

changing.communities.benchpress %>%
  summary() %>%
  dplyr::mutate(Approach = as.character(expression)) %>%
  ggplot(
    aes(p, as.numeric(median), color = Approach, group = Approach)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(n.breaks = 25) +
  labs(x = "Number of Communities",
       y = "Median (s)") +
  theme_bw() +
  theme(legend.position = "right")


