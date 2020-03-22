###############################################################
### Modified functions to allow the computation in parallel ###
###############################################################

#### Mean pairwise distances (ses.MPD), heavily based on picante's original function
sf.ses.mpd <- function(samp, dis, abundance.weighted = FALSE, runs = 999, cores = 4){
  dis <- as.matrix(dis) # coerce's the distance matrix to matrix
  
  # calculate the observed mean phylogenetic distance in each observing unit
  mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted) 
  
  ## build null models in parallel
  sfInit(parallel = TRUE, 
         cpus = cores, 
         slaveOutfile="logAnalysis.txt")
  
  # load libraries in parallel
  sfLibrary(picante)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  
  # export objects 
  sfExport('samp', 'dis', 'abundance.weighted', 'runs')
  
  # calculate expected mean phylogenetic distances
  mpd.rand <- t(sfSapply(1:runs, 
                         function(i){ 
                           sfCat(paste("Iteration ", i), sep="\n")
                           mpd(samp, taxaShuffle(dis),
                               abundance.weighted = abundance.weighted)
                         }))
  
  # 
  mpd.rand.mean <- sfApply(x = mpd.rand, margin = 2, fun = mean, 
                           na.rm = TRUE)
  
  mpd.rand.sd <- sfApply(x = mpd.rand, margin = 2, fun = sd, 
                         na.rm = TRUE)
  
  mpd.obs.rank <- sfApply(x = rbind(mpd.obs, mpd.rand), margin = 2, 
                          fun = rank)[1, ]
  
  sfStop()
  
  ###
  
  mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
  
  mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), 
             mpd.obs, 
             mpd.rand.mean, 
             mpd.rand.sd, 
             mpd.obs.rank, 
             mpd.obs.z, 
             mpd.obs.p = mpd.obs.rank/(runs + 1), 
             runs = runs, 
             row.names = row.names(samp))
}

sf.ses.mntd <- function(samp, dis, abundance.weighted = FALSE, runs = 999, cores = 4){
  dis <- as.matrix(dis)
  mntd.obs <- mntd(samp, dis, abundance.weighted = abundance.weighted)
  
  ###
  sfInit(parallel = TRUE, cpus = cores, slaveOutfile="logAnalysis.txt")
  
  sfLibrary(picante)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  
  sfExport('samp', 'dis', 'abundance.weighted', 'runs')
  
  mntd.rand <- t(sfSapply(1:runs, 
                          function(i){ 
                            sfCat(paste("Iteration ", i), sep="\n")
                            mntd(samp, taxaShuffle(dis),
                                 abundance.weighted = abundance.weighted)
                          }))
  
  mntd.rand.mean <- sfApply(x = mntd.rand, margin = 2, fun = mean, 
                            na.rm = TRUE)
  
  mntd.rand.sd <- sfApply(x = mntd.rand, margin = 2, fun = sd, 
                          na.rm = TRUE)
  
  mntd.obs.rank <- sfApply(x = rbind(mntd.obs, mntd.rand), margin = 2, 
                           fun = rank)[1, ]
  
  sfStop()
  
  ###
  
  mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
  
  mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), 
             mntd.obs, 
             mntd.rand.mean, 
             mntd.rand.sd, 
             mntd.obs.rank, 
             mntd.obs.z, 
             mntd.obs.p = mntd.obs.rank/(runs + 1), 
             runs = runs, 
             row.names = row.names(samp))
}


#system.time(sf.ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), runs = 1000))
#system.time(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), null.model= "taxa.labels", runs = 1000))