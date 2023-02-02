######################################################################################
##  Code to compute average local diversification rates for each community based on ##
# rates estimated for each tip in the phylogeny across many samples from a          ##
# phylogenetic hypothesis                                                            #
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2020-03-21"                                                          #
#                                                                                    # 
######################################################################################

library(coda)
library(tidyverse)
library(foreach)
library(doSNOW)

n.cores = 20

nw <- n.cores  # number of workers
cl <- makeSOCKcluster(nw)
registerDoSNOW(cl)

values.rarefac <- foreach(phy.tree.sample = phylo.sample.pool,
                          #.noexport = ls()
                          .packages = c("BAMMtools", "tibble", "dplyr", 
                                        "ape", "phylobase", "coda")
) %dopar% { 
  
  #### Extract the output from BAMM files #####
  
  dir.sample <- paste0("data/BAMM/", phy.tree.sample, "/")
  
  #### Preparing the phyologenetic tree ####
  # 
  # set.seed(15145562)
  # 
  # phylo.sample.pool <- sample(1:1000, 50)[1:50]
  # 
  # phy.tree.sample <-  phylo.sample.pool[2]
  # 
  
  paste0("data/BAMM/", phy.tree.sample, "/")
  
  #### Extract the output from BAMM files #####
  
  read.tree(
    paste0(dir.sample,
           "Chiroptera.FaurSven.ultra.tree.sampled_",
           phy.tree.sample,
           ".tree")
  )
  
  # Reorder the tree by a series of contiguous rows
  Chiroptera.FaurSven.tree.ultra.cw.sampled <- 
    reorder(
      read.tree(
        paste0(dir.sample,
               "Chiroptera.FaurSven.ultra.tree.sampled_",
               phy.tree.sample,
               ".tree")
      ), 
      "cladewise")
  
  
  ### Get the Event data ###
  
  edata.sampled <- getEventData(Chiroptera.FaurSven.tree.ultra.cw.sampled, 
                                eventdata = paste0(getwd(), "/", dir.sample,
                                                   "event_data_",
                                                   phy.tree.sample,
                                                   ".txt"), 
                                burnin = 0.1)
  
  ### Assess MCMC convergence ####
  mcmc.out.sampled <- read.csv(paste0(getwd(), "/", dir.sample,
                                      "mcmc_out_",
                                      phy.tree.sample,
                                      ".txt"), 
                               header = TRUE)
  
  
  # Discard the first 10% samples as burn-in:
  
  burnstart.sampled <- floor(0.1 * nrow(mcmc.out.sampled))
  postburn.sampled <- mcmc.out.sampled[burnstart.sampled:nrow(mcmc.out.sampled), ]
  
  # Check the effective sample sizes of the log-likelihood and the number of shift
  # events present in each sample
  
  effectiveSize(postburn.sampled$N_shifts)
  effectiveSize(postburn.sampled$logLik)
  
  plot(mcmc.out.sampled$logLik ~ mcmc.out.sampled$generation)
  
  # In general, we want these to be at least 200 (and 200 is on the low side, but
  # might be reasonable for very large datasets).
  
  ####  Estimate individual tip-specific evolutionary rates #### 
  
  ## Obtain speciation, extinction, net diversification for the tips of the tree
  ## from all posterior samples in the output
  
  tip.rates.sampled <- getTipRates(edata.sampled)
  
  str(tip.rates.sampled)
  
  
  #### Calculating CWM for Net Diversification Rates ####
  
  ### CWM for Speciation Rates
  
  lambda.CWM.Chiroptera.Comm.sampled <- CWM_Std_TW(Trait = tip.rates.sampled$lambda.avg, 
                                                   Distrib = Chiroptera.FaurSven.comm)
  
  colnames(lambda.CWM.Chiroptera.Comm.sampled) <- paste("lambda", 
                                                        colnames(lambda.CWM.Chiroptera.Comm.sampled), 
                                                        sep = "_")
  
  mu.CWM.Chiroptera.Comm.sampled <- CWM_Std_TW(Trait = tip.rates.sampled$mu.avg, 
                                               Distrib = Chiroptera.FaurSven.comm)
  
  colnames(mu.CWM.Chiroptera.Comm.sampled) <- paste("mu", 
                                                    colnames(mu.CWM.Chiroptera.Comm.sampled), 
                                                    sep = "_")
  
  ## Subtract mu from lambda to obtain per community net diversification rates
  # (note that this is different then subtracting mu from lambda directly at the
  # tip-rates)
  
  netDiv.CWM.Chiroptera.Comm.sampled <- lambda.CWM.Chiroptera.Comm.sampled - mu.CWM.Chiroptera.Comm.sampled
  
  colnames(netDiv.CWM.Chiroptera.Comm.sampled) <- c("netDiv_CWM", 
                                                    "netDiv_CWM_std_tw", 
                                                    "netDiv_CWM_std_w" )
  
  netDiv.CWM.Chiroptera.sampled <- data.frame(lambda.CWM.Chiroptera.Comm.sampled,
                                              mu.CWM.Chiroptera.Comm.sampled, 
                                              netDiv.CWM.Chiroptera.Comm.sampled) %>% 
    rownames_to_column("ID")
  
  netDiv.CWM.Chiroptera.sampled$ID <- as.numeric(netDiv.CWM.Chiroptera.sampled$ID)
  
  write.csv(netDiv.CWM.Chiroptera.sampled,
            paste0(getwd(), "/", dir.sample,
                   "netDiv_CWM_std_tw_",
                   phy.tree.sample,
                   ".csv")
  )
  
}

# Access and read all files

netDiv.CWM.Chiroptera.loaded_sampled <- data.frame(
  ID = as.numeric(netDiv.CWM.Chiroptera.sampled$ID)
)

netDiv_CWM_std_tw_files <- list.files(pattern="netDiv_CWM_std_tw_", 
                                      recursive = TRUE)

netDiv_CWM_std_tw_files.sample <- netDiv_CWM_std_tw_files[1]

for(netDiv_CWM_std_tw_files.sample in netDiv_CWM_std_tw_files){
  
  netDiv.CWM.Chiroptera.sampled <- read.csv(netDiv_CWM_std_tw_files.sample, 
                                            row.names = 1)
  
  netDiv.CWM.Chiroptera.loaded_sampled <- netDiv.CWM.Chiroptera.loaded_sampled %>% 
    left_join(
      netDiv.CWM.Chiroptera.sampled %>%
        dplyr::select(ID, netDiv_CWM_std_tw), 
      by = "ID")
}

netDiv.CWM.Chiroptera.mean <- netDiv.CWM.Chiroptera.loaded_sampled %>%
  mutate(netDiv_CWM_std_tw_mean = rowMeans(
    dplyr::select(netDiv.CWM.Chiroptera.loaded_sampled, 
                  starts_with("netDiv_CWM_std")), 
    na.rm = TRUE)
  ) %>%
  dplyr::select(ID, netDiv_CWM_std_tw_mean)

colnames(netDiv.CWM.Chiroptera.mean)[2] <- "netDiv_CWM_std_tw_rs_1_mean"

hist(netDiv.CWM.Chiroptera.mean$netDiv_CWM_std_tw_rs_1_mean)

write.csv(netDiv.CWM.Chiroptera.mean,
          "data/matrices/netDiv.CWM.Chiroptera.mean.csv")

# Merging data frames ####

### Robustness to phylogenetic uncertainty ###

netDiv.CWM.Chiroptera.rs_1.mean <- read.csv("data/matrices/netDiv.CWM.Chiroptera.rs_1.mean.csv", 
                                            row.names = 1)

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  left_join(netDiv.CWM.Chiroptera.rs_1.mean, 
            by = "ID")


netDiv.CWM.Chiroptera.rs_5.mean <- read.csv("data/matrices/netDiv.CWM.Chiroptera.rs_5.mean.csv", 
                                            row.names = 1)

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  left_join(netDiv.CWM.Chiroptera.rs_5.mean, 
            by = "ID")


# A few verifications

left_join(
  MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    dplyr::select(ID, netDiv_CWM_std_tw_rs_1) %>%
    drop_na,
  MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    dplyr::select(ID, netDiv_CWM_std_tw_rs_1_mean) %>%
    drop_na,
  by = "ID"
)[1000:1600, ]

head(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div)

dim(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div)


## Exploration


MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2) %>%
  dplyr::select(netDiv_CWM_std_tw_rs_1_mean, netDiv_CWM_std_tw_rs_5_mean) %>%
  plot

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2) %>%
  dplyr::select(netDiv_CWM_std_tw_rs_1, 
                netDiv_CWM_std_tw_rs_1_mean) %>%
  plot

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.rarefaction.relative %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2) %>%
  #  select(nri.sec) %>%
  ggplot() +
  geom_histogram(aes(x = nti.sec),
                 bins  = 50)

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.rarefaction.relative.sampling %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2) %>%
  #  select(nri.sec) %>%
  ggplot() +
  geom_histogram(aes(x = nti.sec),
                 bins  = 50)



# write.csv(netDiv.CWM.Chiroptera.Comm, "data/matrices/netDiv.CWM.Chiroptera.Comm.csv")

# Explore units that have more than 2 co-occurring taxa.

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2)

# write.csv(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div, "data/matrices/MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div.sqrt.csv")


# MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- read.csv("data/matrices/MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div.csv", h = T)

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$ID_Realm, 
                                                                                  levels = c('Neotropical',
                                                                                             'Nearctic',
                                                                                             'Afrotropical',
                                                                                             'Palearctic', 
                                                                                             'Indomalay', 
                                                                                             'Australasian'))

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$SamplingPool <-  factor(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$SamplingPool,  
                                                                                       levels = c("Global sampling",
                                                                                                  "Hemispheric sampling",
                                                                                                  "Realm sampling",
                                                                                                  "Plate sampling",
                                                                                                  "Biome sampling",
                                                                                                  "Ecoregion sampling"))


saveRDS(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div,
        "data/matrices/MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div.RDS")


