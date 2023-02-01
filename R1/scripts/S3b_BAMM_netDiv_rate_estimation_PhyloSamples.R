######################################################################################
### Code to estimate extinction and speciation rates for bats from trees sampled   ###
### from the posterior distribution of phylogenetic relationships                  ###
#
# This code is dependent on the function BAMM.diversification.rate.estimation()      #
# available within the file S0_fun_BAMM.diversification.rate.estimation.R. 
# The estimation process below uses snowfall to apply the
# BAMM.diversification.rate.estimation() to each phylogenetic tree in the pool of
# sampled phylogenetic trees from the posterior distribution made available by Faurby
# and Svenning 2015.                                                                 #
# 
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2020-03-21"                                                          #
#                                                                                    # 
######################################################################################

# Setting seed each time one uses sample is very important so the trees that are
# sampled are always the same when this procedure can be reproduced to output
# the same results as the ones in the submitted document.

set.seed(15145562)

# Sampled 50 numbers from 1 to 1000, which will be used as indices to pull the
# phylogenetic trees from the posterior distributions

phylo.sample.pool <- sample(1:1000, 100)[1:100]

# For testing: phy.tree.sample <-  phylo.sample.pool[2]

# Begin of the parallel application of BAMM.diversification.rate.estimation()
# use ?snowfall to understand how it works.

# Note that BAMM.diversification.rate.estimation() already has a parallelized
# component within it (for the MCMC chains). The desired number needs to be set
# up there, and one must consider that the number of cores in sfInit() (the cpus
# argument) needs to be multiplied by the number of MCMC chains that are running
# in parallel within BAMM.diversification.rate.estimation().

# If sfInit(cpus = 4) and the numberOfChains within
# BAMM.diversification.rate.estimation() is 6, then 24 cores will used.

sfInit(parallel = TRUE, 
       cpus = 8, 
       slaveOutfile = "logAnalysis.txt")

# Load libraries in parallel

sfLibrary(BAMMtools)
sfLibrary(base)
sfLibrary(phytools)
sfLibrary(snowfall)
sfLibrary(picante)
sfLibrary(dplyr)
sfLibrary(coda)

# Export objects to the CPUs
sfExport('Mammal.FaurSven.trees', 
         'Chiroptera.Comm',
         'phylo.sample.pool',
         'BAMM.diversification.rate.estimation')

# Apply the BAMM.diversification.rate.estimation

# ess.BAMM.diversification.rate.estimation outputs a list with effective sample
# sizes calculated from the MCMC output of the diversification rate estimation.
# Adequate estimations tend to have these ESS values above 200 units.

ess.BAMM.diversification.rate.estimation <- sfLapply(
  phylo.sample.pool[1:100], 
  fun = BAMM.diversification.rate.estimation,
  Mammal.FaurSven.trees = Mammal.FaurSven.trees,
  Chiroptera.Comm = Chiroptera.Comm
)

sfStop()






