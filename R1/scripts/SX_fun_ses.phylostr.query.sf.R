#################################################################################
### Modification of the functions picante::ses.mpd() and picante::ses.mntd()  ###
### to allow for parallel computation using SNOW/snowfall                     ###
#                                                                               #
# Description: This code presents the sf.ses.mpd() and sf.ses.mntd() functions  #
# that are heavily-based in picante's ses.mpd() and ses.mntd() functions. It    #
# implements parallel computation for the calculation of the mean pairwise      #
# distances and mean nearest taxon distances using SNOW via the snowfall        # 
# package. The current implementation only addresses the taxa.labels null model.#
#                                                                               #
# Note: Future modifications should include the other null models, as well as   #
# improvements in speed supported by benchmarking.                              #
#                                                                               #
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2021-07-02"                                                     #
#                                                                               # 
#################################################################################

mod.ses.mpd.query.sf <- function(phylo.tree = phylo.tree,
                                 perms = 999,
                                 species.data = species.data,
                                 cores = 20){
  
  mpd.obs.query <- mpd.query(tree = phylo.tree,
                             matrix = species.data)
  
  # Replace the mpd of all communities that have 1 or less species for NA
  
  mpd.obs.query <- replace(x = mpd.obs.query,  
                           rowSums(species.data > 0) <= 1, 
                           values = NA)
  
  # system.time({mpd.obs.query.ppn <- mpd.ppn(phylo.dist = cophenetic(phylo.tree),
  #                                          species.data = species.data) })
  # 
  
  taxaShuffle.map <- purrr::map(1:perms, 
                                ~shuffle_tiplabels(phylo.tree))
  
  # format(object.size(taxaShuffle.map.l), 
  #        units = "Gb",
  #        standard = "auto",
  #        digits = 5L)
  
  # build null models in parallel
  
  sfInit(parallel = TRUE, 
         cpus = cores, 
         slaveOutfile = "logAnalysis.txt")
  
  # load libraries in parallel
  
  sfLibrary(PhyloMeasures)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  
  # export objects 
  sfExport('mpd.obs.query', 
           'perms',
           'taxaShuffle.map')
  
  random.mpd.query <- sfLapply(x = taxaShuffle.map,
                               fun = mpd.query,
                               matrix = species.data)
  
  random.mpd.query <- rlist::list.rbind(random.mpd.query)  
  
  random.mpd.query <-  t(apply(random.mpd.query,
                               1,
                               FUN = replace,
                               rowSums(species.data > 0) <= 1, 
                               values = NA))
  
  random.mpd.query.mean <- sfApply(x = random.mpd.query, 
                                   margin = 2,
                                   fun = mean,
                                   na.rm = TRUE)
  
  random.mpd.query.sd <- sfApply(random.mpd.query, 
                                 margin = 2, # by rows
                                 fun = sd,
                                 na.rm = TRUE)
  
  mpd.obs.rank <- sfApply(x = rbind(mpd.obs.query, 
                                    random.mpd.query.mean), 
                          margin = 2, 
                          fun = rank)[1, ]
  
  mpd.obs.rank <- ifelse(is.na(random.mpd.query.mean), 
                         NA, 
                         mpd.obs.rank)
  
  sfStop()
  
  names(mpd.obs.query) <- row.names(species.data)
  
  ses.mpd.z.query <- (mpd.obs.query - random.mpd.query.mean)/random.mpd.query.sd
  
  return(
    data.frame(ntaxa = specnumber(species.data),
               mpd.obs.query,
               random.mpd.query.mean,
               random.mpd.query.sd,
               ses.mpd.z.query,
               mpd.obs.p = mpd.obs.rank/(perms + 1),
               nri = -1 * (mpd.obs.query - random.mpd.query.mean)/random.mpd.query.sd,
               perms,
               ntaxa.pool.size = length(phylo.tree$tip.label))
  )
}

mod.ses.mntd.query.sf <- function(phylo.tree = phylo.tree,
                                  perms = 999,
                                  species.data = species.data,
                                  cores = 20){
  
  mntd.obs.query <- mntd.query(tree = phylo.tree,
                               matrix = species.data)
  
  # Replace the mntd of all communities that have 1 or less species for NA
  
  mntd.obs.query <- replace(x = mntd.obs.query,  
                            rowSums(species.data > 0) <= 1, 
                            values = NA)
  
  # system.time({mntd.obs.query.ppn <- mntd.ppn(phylo.dist = cophenetic(phylo.tree),
  #                                          species.data = species.data) })
  # 
  
  taxaShuffle.map <- purrr::map(1:perms, 
                                ~shuffle_tiplabels(phylo.tree))
  
  # format(object.size(taxaShuffle.map.l), 
  #        units = "Gb",
  #        standard = "auto",
  #        digits = 5L)
  
  # build null models in parallel
  
  sfInit(parallel = TRUE, 
         cpus = cores, 
         slaveOutfile = "logAnalysis.txt")
  
  # load libraries in parallel
  
  sfLibrary(PhyloMeasures)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  
  # export objects 
  sfExport('mntd.obs.query', 
           'perms',
           'taxaShuffle.map')
  
  random.mntd.query <- sfLapply(x = taxaShuffle.map,
                                fun = mntd.query,
                                matrix = species.data)
  
  random.mntd.query <- rlist::list.rbind(random.mntd.query)  
  
  random.mntd.query <-  t(apply(random.mntd.query,
                                1,
                                FUN = replace,
                                rowSums(species.data > 0) <= 1, 
                                values = NA))
  
  random.mntd.query.mean <- sfApply(x = random.mntd.query, 
                                    margin = 2,
                                    fun = mean,
                                    na.rm = TRUE)
  
  random.mntd.query.sd <- sfApply(random.mntd.query, 
                                  margin = 2, # by rows
                                  fun = sd,
                                  na.rm = TRUE)
  
  mntd.obs.rank <- sfApply(x = rbind(mntd.obs.query, 
                                     random.mntd.query.mean), 
                           margin = 2, 
                           fun = rank)[1, ]
  
  mntd.obs.rank <- ifelse(is.na(random.mntd.query.mean), 
                          NA, 
                          mntd.obs.rank)
  
  sfStop()
  
  names(mntd.obs.query) <- row.names(species.data)
  
  ses.mntd.z.query <- (mntd.obs.query - random.mntd.query.mean)/random.mntd.query.sd
  
  return(
    data.frame(ntaxa = specnumber(species.data),
               mntd.obs.query,
               random.mntd.query.mean,
               random.mntd.query.sd,
               ses.mntd.z.query,
               mntd.obs.p = mntd.obs.rank/(perms + 1),
               nti = -1 * (mntd.obs.query - random.mntd.query.mean)/random.mntd.query.sd,
               perms,
               ntaxa.pool.size = length(phylo.tree$tip.label))
  )
}

# 
# system.time({mod.ses.mntd.query.sf(phylo.tree = phylo.tree,
#                                    perms = 99,
#                                    species.data = species.data,
#                                    cores = 20)})


# system.time({mod.ses.mpd.query.sf(phylo.tree = phylo.tree,
#                                   perms = 999,
#                                   species.data = species.data,
#                                   cores = 20)})
