######################################################################################
## Utility functions to compute the phylogenetic relatedness of communities
## using the traditional and a new rarefaction approach based on the
## standardized effect sizes of mean phylogenetic pairwise and mean nearest
## taxon distances (based on PhyloMeasures::mpd.query() and
## PhyloMeasures::mntd.query()).
######################################################################################

# mpd

opt.rarefaction.mpd <- function(phylo.tree,
                                species.data,
                                n.species.rarefac,
                                n.rep = 100, 
                                n.cores = 4) {
  
  # species.data represent set of communities in the region
  
  species.names <- colnames(species.data)
  n.species <- ncol(species.data)
  
  # we compute the square-root of the cophenetic distances of phylo.tree, and
  # then convert the distances back to a phylogenetic tree
  
  ses.mpd.z.query <- mpd.query(as.phylo(hclust(as.dist(sqrt(cophenetic(phylo.tree))))),
                               species.data,
                               standardize = TRUE)
  
  names(ses.mpd.z.query) <- row.names(species.data)
  
  # PhyloMeasures::mpd.query() computes mpd = 0 for communities that had less or
  # equal than one species occurrence.
  # The line below replaces the values for communities with less or equal than
  # one species occurrence to NA.
  
  ses.mpd.z.query <- replace(x = ses.mpd.z.query,
                             rowSums(species.data > 0) <= 1,
                             values = NA)
  
  # library(doParallel)
  
  library(foreach)
  library(doSNOW)
  
  # library(tcltk)
  
  nw <- n.cores  # number of workers
  cl <- makeSOCKcluster(nw)
  registerDoSNOW(cl)
  
  # library(doFuture)
  # registerDoParallel(cores = n.cores)
  
  # registerDoFuture()
  # plan(multicore, workers = n.cores)
  
  values.rarefac <- foreach(i = 1:n.rep, 
                            #.noexport = ls()
                            .packages = c("ape", "PhyloMeasures")
                            ) %dopar% {
    choose.species.rarefac <- species.names[sample(n.species)[1:n.species.rarefac]]
    species.data.rarefac <- species.data[, choose.species.rarefac]
    phylo.tree.rarefac <- keep.tip(phylo.tree, 
                                   choose.species.rarefac)
    mpd <- mpd.query(as.phylo(hclust(as.dist(sqrt(cophenetic(phylo.tree.rarefac))))), 
                     species.data.rarefac, 
                     standardize = TRUE)
    
    names(mpd) <- row.names(species.data)
    
    mpd <- replace(
      x = mpd,
      rowSums(species.data > 0) <= 1,
      values = NA
    )
    
    return(list(mpd))
  }
  
  # Stop implicit clusters
  
  stopCluster(cl)
  
  ses.mpd.z.query.rarefac.matrix <- sapply(values.rarefac, function(x) {
    x[[1]]
  })
  
  # Calculate the mean ses mpd from the rarefacted subsets
  
  ses.mpd.z.query.rarefac.mean <- ses.mpd.z.query.rarefac.matrix %>%
    as.data.frame() %>%
    mutate(ses.mpd.z.query.rarefac.mean = rowMeans(., na.rm = TRUE)) %>%
    dplyr::select(ses.mpd.z.query.rarefac.mean)
  
  # 
  
  ses.mpd.z.query.rarefac <- ses.mpd.z.query.rarefac.mean %>%
    mutate(
      ses.mpd.z.query.sec = ses.mpd.z.query,
      nri.sec = -1 * ses.mpd.z.query,
      nri.rarefac.mean = -1 * ses.mpd.z.query.rarefac.mean,
      rarefac.pool.size = n.species.rarefac
    )
  
  output <- list(
    ses.mpd.z.query.rarefac.matrix = ses.mpd.z.query.rarefac.matrix,
    ses.mpd.z.query.rarefac = ses.mpd.z.query.rarefac
  )
  return(output)
}

# mntd

opt.rarefaction.mntd <- function(phylo.tree,
                                 species.data,
                                 n.species.rarefac,
                                 n.rep = 100, 
                                 n.cores = 4) {
  
  # species.data represent set of communities in the region
  
  species.names <- colnames(species.data)
  n.species <- ncol(species.data)
  
  # we compute the square-root of the cophenetic distances of phylo.tree, and
  # then convert the distances back to a phylogenetic tree
  
  ses.mntd.z.query <- mntd.query(as.phylo(hclust(as.dist(sqrt(cophenetic(phylo.tree))))),
                                 species.data,
                                 standardize = TRUE)
  
  names(ses.mntd.z.query) <- row.names(species.data)
  
  # PhyloMeasures::mntd.query() computes mntd = 0 for communities that had less or
  # equal than one species occurrence.
  # The line below replaces the values for communities with less or equal than
  # one species occurrence to NA.
  
  ses.mntd.z.query <- replace(x = ses.mntd.z.query,
                              rowSums(species.data > 0) <= 1,
                              values = NA)
  
  # library(doParallel)
  
  library(foreach)
  library(doSNOW)
  
  # library(tcltk)
  
  nw <- n.cores  # number of workers
  cl <- makeSOCKcluster(nw)
  registerDoSNOW(cl)
  
  # library(doFuture)
  # registerDoParallel(cores = n.cores)
  
  # registerDoFuture()
  # plan(multicore, workers = n.cores)
  
  values.rarefac <- foreach(i = 1:n.rep, 
                            #.noexport = ls()
                            .packages = c("ape", "PhyloMeasures")
  ) %dopar% {
    choose.species.rarefac <- species.names[sample(n.species)[1:n.species.rarefac]]
    species.data.rarefac <- species.data[, choose.species.rarefac]
    phylo.tree.rarefac <- keep.tip(phylo.tree, 
                                   choose.species.rarefac)
    mntd <- mntd.query(as.phylo(hclust(as.dist(sqrt(cophenetic(phylo.tree.rarefac))))), 
                       species.data.rarefac, 
                       standardize = TRUE)
    
    names(mntd) <- row.names(species.data)
    
    mntd <- replace(
      x = mntd,
      rowSums(species.data > 0) <= 1,
      values = NA
    )
    
    return(list(mntd))
  }
  
  # Stop implicit clusters
  
  stopCluster(cl)
  
  ses.mntd.z.query.rarefac.matrix <- sapply(values.rarefac, function(x) {
    x[[1]]
  })
  
  # Calculate the mean ses mntd from the rarefacted subsets
  
  ses.mntd.z.query.rarefac.mean <- ses.mntd.z.query.rarefac.matrix %>%
    as.data.frame() %>%
    mutate(ses.mntd.z.query.rarefac.mean = rowMeans(., na.rm = TRUE)) %>%
    dplyr::select(ses.mntd.z.query.rarefac.mean)
  
  # 
  
  ses.mntd.z.query.rarefac <- ses.mntd.z.query.rarefac.mean %>%
    mutate(
      ses.mntd.z.query.sec = ses.mntd.z.query,
      nti.sec = -1 * ses.mntd.z.query,
      nti.rarefac.mean = -1 * ses.mntd.z.query.rarefac.mean,
      rarefac.pool.size = n.species.rarefac
    )
  
  output <- list(
    ses.mntd.z.query.rarefac.matrix = ses.mntd.z.query.rarefac.matrix,
    ses.mntd.z.query.rarefac = ses.mntd.z.query.rarefac
  )
  return(output)
}

