######################################################################################
## Utility wrapper function to estimate extinction and speciation rates for bats    ##
# using Bayesian Analyses for Macroevolutionary Mixtures                            ##
######################################################################################

# The BAMM.diversification.rate.estimation() function requires:
# a multiPhylo object with containing multiple phylogenetic trees; 
# a number indicating the position of a tree of interest within the multiPhylo object
# a community composition matrix that will be used to prune the tree to the same
# species as the community.
#
# This routine begins by creating a directory indexed and named after the position of
# the tree within the pool of phylogenetic trees available for the computation. The
# routine will then prune the phylogenetic tree to the community composition, force
# branch lengths below or equal to zero to zero and make the tree ultrametric.
#
# Then, the function will estimate the prior block for speciation and extinction rates
# using BAMMtools::setBAMMpriors().
#
# The BAMM control file will the be generated with BAMMtools::generateControlFile().
#
# After this step, the BAMM model will be passed to the system command terminal
# (Linux) directly from R.
#
# All control and output files will be exported to the directory named after the
# position of the tree used in the pool phylogenetic trees.
# 
# This function can be applied in a loop, used with lapply() or in parallel. 

######################################################################################
# Author: Pedro Henrique Pereira Braga                                               # 
# Last Update: "2022-01-18"                                                          #
######################################################################################


# library(BAMMtools)
# library(coda)

#### Preparing the phyologenetic tree ####
# 
# set.seed(15145562)
# 
# phylo.sample.pool <- sample(1:1000, 50)[1:25]
# 
# phy.tree.sample <-  phylo.sample.pool[2]
# 

BAMM.diversification.rate.estimation <- function(Mammal.FaurSven.trees = Mammal.FaurSven.trees,
                                                 Chiroptera.Comm = Chiroptera.Comm,
                                                 phy.tree.sample = phy.tree.sample) {
  
  dir.create(paste0("data/BAMM/", phy.tree.sample))
  
  dir.sample <- paste0("data/BAMM/", phy.tree.sample, "/")
  
  Chiroptera.FaurSven.tree.sampled <- match.phylo.comm(Mammal.FaurSven.trees[[phy.tree.sample]], 
                                                       comm = Chiroptera.Comm)[["phy"]]
  
  # Check if the tree is ultrametric:
  is.ultrametric(Chiroptera.FaurSven.tree.sampled)
  
  # Check if it is dichotomous, i.e. whether every node (including the root node)
  # has exactly two descendant nodes.
  is.binary(Chiroptera.FaurSven.tree.sampled)
  
  # Count branches with negative lengths
  sum(Chiroptera.FaurSven.tree.sampled$edge.length < 0)
  
  # Count branches with lengths equal to zero
  sum(Chiroptera.FaurSven.tree.sampled$edge.length == 0)
  
  
  Chiroptera.FaurSven.ultra.tree.sampled <- force.ultrametric(Chiroptera.FaurSven.tree.sampled)
  
  Chiroptera.FaurSven.ultra.tree.sampled %>% is.ultrametric()
  Chiroptera.FaurSven.ultra.tree.sampled %>% is.binary()
  
  # For negative branches to be zero
  Chiroptera.FaurSven.ultra.tree.sampled$edge.length[Chiroptera.FaurSven.ultra.tree.sampled$edge.length < 0] <- 0
  
  # Force zero lengths to be a very small number
  Chiroptera.FaurSven.ultra.tree.sampled$edge.length[Chiroptera.FaurSven.ultra.tree.sampled$edge.length == 0] <- 0.000001
  
  # Count branches with negative lengths
  sum(Chiroptera.FaurSven.ultra.tree.sampled$edge.length < 0)
  
  # Count branches with lengths equal to zero
  sum(Chiroptera.FaurSven.ultra.tree.sampled$edge.length == 0)
  
  # Export tree
  
  write.tree(Chiroptera.FaurSven.ultra.tree.sampled, 
             paste0(dir.sample,
                    "Chiroptera.FaurSven.ultra.tree.sampled_",
                    phy.tree.sample,
                    ".tree"))
  
  #### Preparing BAMM control file
  
  # Estimate the prior block using BAMMtools::setBAMMpriors()
  
  # In a nutshell, et the main lambdaInit and muInit priors, setBAMMpriors first
  # estimates the rate of speciation for your full tree under a pure birth model
  # of diversification.  a reasonable prior distribution for the initial rate
  # parameters is an exponential distribution with a mean five times greater than
  # this pure birth value. For betaInitPrior and betaInitRootPrior,  rather than
  # fitting a pure-birth model, we find the maximum likelihood estimate of the
  # variance parameter under a Brownian motion model.
  
  ### Correcting for Sampling Fractions
  
  # Copying the file from the MCC dir
  
  file.copy("data/BAMM/MCC/samplingProbsBats.txt", 
            paste0(dir.sample,
                   "samplingProbsBats.txt"))
  
  
  #### Preparing BAMM control file
  
  # Estimate the prior block using BAMMtools::setBAMMpriors()
  
  # In a nutshell, et the main lambdaInit and muInit priors, setBAMMpriors first
  # estimates the rate of speciation for your full tree under a pure birth model
  # of diversification.  a reasonable prior distribution for the initial rate
  # parameters is an exponential distribution with a mean five times greater than
  # this pure birth value. For betaInitPrior and betaInitRootPrior,  rather than
  # fitting a pure-birth model, we find the maximum likelihood estimate of the
  # variance parameter under a Brownian motion model.
  
  (BAMMpriors <- setBAMMpriors(Chiroptera.FaurSven.ultra.tree.sampled,	
                               outfile = NULL))
  
  ## Generate the control file
  generateControlFile(paste0(dir.sample,
                             "chiroptera.divcontrol.",
                             phy.tree.sample,
                             ".txt"), 
                      type = "diversification", 
                      params = list(treefile = paste0(getwd(), "/", dir.sample,
                                                      "Chiroptera.FaurSven.ultra.tree.sampled_",
                                                      phy.tree.sample,
                                                      ".tree"),
                                    runInfoFilename = paste0(getwd(), "/", dir.sample,
                                                             "run_info_",
                                                             phy.tree.sample,
                                                             ".txt"),
                                    eventDataOutfile = paste0(getwd(), "/", dir.sample,
                                                              "event_data_",
                                                              phy.tree.sample,
                                                              ".txt"),
                                    mcmcOutfile = paste0(getwd(), "/", dir.sample,
                                                         "mcmc_out_",
                                                         phy.tree.sample,
                                                         ".txt"),
                                    chainSwapFileName = paste0(getwd(), "/", dir.sample,
                                                               "chain_swap_",
                                                               phy.tree.sample,
                                                               ".txt"),
                                    numberOfGenerations = "20000000",
                                    useGlobalSamplingProbability = 0, 
                                    
                                    sampleProbsFilename =  paste0(getwd(), "/", dir.sample,
                                                                  "samplingProbsBats.txt"),
                                    overwrite = "1", 
                                    lambdaInitPrior = as.numeric(BAMMpriors["lambdaInitPrior"]), 
                                    lambdaShiftPrior = as.numeric(BAMMpriors["lambdaShiftPrior"]), 
                                    muInitPrior = as.numeric(BAMMpriors["muInitPrior"]), 
                                    numberOfChains = 8,
                                    seed = 15145562,
                                    # Results were similar for "1", "10" and "25"
                                    expectedNumberOfShifts = "1", 
                                    # Allow shifts to occur in all branches; thus for
                                    # estimates for terminal branches, in which we are
                                    # interested.
                                    minCladeSizeForShift = 1) 
  )
  
  
  #### Running the BAMM model in the shell ####
  
  # In the shell, run the following code. Note that it will take at least 24 hours
  # under 38 cores and 2.2 GHz clock.
  
  # I highly recommend the use of 'screen' if done remotely.
  
  # system(paste('renice', 19, Sys.getpid()), ignore.stdout = T)
  
  system(paste0("nice -n 10 bamm -c ",
                getwd(), "/", dir.sample,
                "chiroptera.divcontrol.",
                phy.tree.sample,
                ".txt"))
  
  if(file.exists(paste0(getwd(), "/", dir.sample,
                        "mcmc_out_",
                        phy.tree.sample,
                        ".txt"))){
    mcmcout <- read.csv(paste0(getwd(), "/", dir.sample,
                               "mcmc_out_",
                               phy.tree.sample,
                               ".txt"), 
                        header = TRUE)
    
    # Discard the first 10% samples as burn-in:
    
    burnstart <- floor(0.1 * nrow(mcmcout))
    postburn <- mcmcout[burnstart:nrow(mcmcout), ]
    
    # Check the effective sample sizes of the log-likelihood and the number of shift
    # events present in each sample
    
    ess_postburn_n_shifts <- coda::effectiveSize(postburn$N_shifts)
    ess_postburn_logLik <- coda::effectiveSize(postburn$logLik)
    
    ess_results <- data.frame(phy.tree.sample = phy.tree.sample,
                              ess_postburn_logLik = ess_postburn_logLik,
                              ess_postburn_n_shifts = ess_postburn_n_shifts)
  }
}

