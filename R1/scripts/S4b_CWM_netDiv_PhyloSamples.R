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
