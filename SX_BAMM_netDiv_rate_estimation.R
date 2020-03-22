#################################################################################
### Code to estimate extinction and speciation rates for bats using           ###
### Bayesian Analyses for Macroevolutionary Mixtures                          ###
#                                                                               #
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2020-03-21"                                                     #
#                                                                               # 
#################################################################################

# library(BAMMtools)
# library(coda)

#### Preparing the phyologenetic tree ####

# Force the tree to be ultrametric by computing the set of edge lengths that
# result in a minimized sum-of-squares distance between the patristic distance
# of the output and input trees (method="nnls");

Chiroptera.FaurSven.tree.ultra <- Chiroptera.FaurSven.tree

# Check if the tree is ultrametric:
is.ultrametric(Chiroptera.FaurSven.tree.ultra)

# Check if it is dichotomous, i.e. whether every node (including the root node)
# has exactly two descendant nodes.
is.binary.tree(Chiroptera.FaurSven.tree.ultra)

# Count branches with negative lengths
sum(Chiroptera.FaurSven.tree.ultra$edge.length < 0)

# Count branches with lengths equal to zero
sum(Chiroptera.FaurSven.tree.ultra$edge.length == 0)

## Force changes to the tree

# For negative branches to be zero
Chiroptera.FaurSven.tree.ultra$edge.length[Chiroptera.FaurSven.tree.ultra$edge.length < 0] <- 0

# Force zero lengths to be a very small number
Chiroptera.FaurSven.tree.ultra$edge.length[Chiroptera.FaurSven.tree.ultra$edge.length == 0] <- 0.000001

# Alternative: convert zero lengths to 1/365
Chiroptera.FaurSven.tree.ultra$edge.length <- pmax(Chiroptera.FaurSven.tree.ultra$edge.length, 1/365)


Chiroptera.FaurSven.tree.ultra <- force.ultrametric(Chiroptera.FaurSven.tree.ultra); is.ultrametric(Chiroptera.FaurSven.tree.ultra)

# Count branches with negative lengths
sum(Chiroptera.FaurSven.tree.ultra$edge.length < 0)

# Count branches with lengths equal to zero
sum(Chiroptera.FaurSven.tree.ultra$edge.length == 0)


# Export tree
write.tree(Chiroptera.FaurSven.tree.ultra, "data/BAMM/Chiroptera.FaurSven.tree.ultra.tree")

# cat(readLines("data/BAMM/Chiroptera.FaurSven.tree.ultra.tree"))

Chiroptera.FaurSven.tree.ultra <- read.tree("data/BAMM/Chiroptera.FaurSven.tree.ultra.tre")
is.ultrametric(Chiroptera.FaurSven.tree.ultra)

#### Preparing BAMM control file

# Estimate the prior block using BAMMtools::setBAMMpriors()

# In a nutshell, et the main lambdaInit and muInit priors, setBAMMpriors first
# estimates the rate of speciation for your full tree under a pure birth model
# of diversification.  a reasonable prior distribution for the initial rate
# parameters is an exponential distribution with a mean five times greater than
# this pure birth value. For betaInitPrior and betaInitRootPrior,  rather than
# fitting a pure-birth model, we find the maximum likelihood estimate of the
# variance parameter under a Brownian motion model.
(BAMMpriors <- setBAMMpriors(Chiroptera.FaurSven.tree.ultra,	
                             outfile = NULL))

### Correcting for Sampling Fractions

separateGenusSp <- stringr::str_split_fixed(Chiroptera.FaurSven.tree.ultra$tip.label, "_", 2)

sppTable <- as.data.frame(matrix(nrow=length(Chiroptera.FaurSven.tree.ultra$tip.label),
                                 ncol=2))
sppTable[ , 1] <- Chiroptera.FaurSven.tree.ultra$tip.label
sppTable[ , 2] <- separateGenusSp[, 1]

ASM_MammalsData <- read.csv("/home/pedro.braga/chapter-BiogHistPhyloRelatedSpatScales/chapter-BiogHistPhyloRelatedSpatScales/data/matrices/asm-all-species-2020-03-20.csv", 
                            header = TRUE)

# colnames(ASM_MammalsData)

ASM_BatsData <- ASM_MammalsData %>%
  filter(Linnean.Order == "Chiroptera") %>%
  select(Genus, Species, Linnean.Family) %>%
  droplevels()

## Correct for nomenclature changes
ASM_BatsData <- ASM_BatsData %>%
  mutate(Genus = recode(Genus, 
                        "Neoplatymops" = "Molossops",
                        "Aeorestes" = "Lasiurus",
                        "Dasypterus" = "Lasiurus",
                        "Cabreramops" = "Molossops",
                        "Cistugo" = "Myotis",
                        "Desmalopex" = "Pteropus",
                        "Doryrhina" = "Hipposideros",
                        "Gardnerycteris" = "Mimon",
                        "Hsunycteris" = "Lonchophylla",
                        "Hypsugo" = "Pipistrellus",
                        "Koopmania" = "Artibeus",
                        "Macronycteris" = "Hipposideros",
                        "Mirimiri" = "Pteralopex",
                        "Parastrellus" = "Pipistrellus",
                        "Perimyotis" = "Pipistrellus",
                        "Rhyneptesicus" = "Eptesicus",
                        "Stenonycteris" = "Rousettus",
                        "Submyotodon" = "Myotis")) %>%
  add_row(Genus = "Lissonycteris", 
          Species = "angolensis", 
          Linnean.Family = "Pteropodidae") %>%
  filter(!(Genus %in% setdiff(unique(Genus), sppTable[ , 2]))) %>%
  droplevels()

(ASM_speciesRichnessGenus <- ASM_BatsData %>% 
    group_by(Genus) %>% 
    summarize(count=n()) %>%
    deframe() 
)

(FaurSven_speciesRichnessGenus <- sppTable %>%
    group_by(V2) %>% 
    summarize(count=n()) %>%
    deframe()
)

# setdiff(names(speciesRichnessGenus), sppTable[ , 2])
# setdiff(sppTable[ , 2], names(speciesRichnessGenus))

# Create a sampling fraction file
samplingProbsBats <- samplingProbs(Chiroptera.FaurSven.tree.ultra, 
                                   cladeTable = sppTable, 
                                   cladeRichness = ASM_speciesRichnessGenus,
                                   globalSampling = 0.8,
                                   writeToDisk = TRUE,
                                   output = "/home/pedro.braga/chapter-BiogHistPhyloRelatedSpatScales/chapter-BiogHistPhyloRelatedSpatScales/data/BAMM/samplingProbsBats.txt")

# To output it to the drive
samplingProbsBats <- samplingProbs(Chiroptera.FaurSven.tree.ultra, 
                                   cladeTable = sppTable, 
                                   cladeRichness = ASM_speciesRichnessGenus,
                                   globalSampling = 0.8,
                                   writeToDisk = FALSE)
#### Preparing BAMM control file

# Estimate the prior block using BAMMtools::setBAMMpriors()

# In a nutshell, et the main lambdaInit and muInit priors, setBAMMpriors first
# estimates the rate of speciation for your full tree under a pure birth model
# of diversification.  a reasonable prior distribution for the initial rate
# parameters is an exponential distribution with a mean five times greater than
# this pure birth value. For betaInitPrior and betaInitRootPrior,  rather than
# fitting a pure-birth model, we find the maximum likelihood estimate of the
# variance parameter under a Brownian motion model.
(BAMMpriors <- setBAMMpriors(Chiroptera.FaurSven.tree.ultra,	
                             outfile = NULL))

## Generate the control file
generateControlFile("data/BAMM/SamplingCorrected/chiroptera.divcontrol.txt", 
                    type = "diversification", 
                    params = list(treefile = paste0(getwd(), "/data/BAMM/", "Chiroptera.FaurSven.tree.ultra.tre"),
                                  runInfoFilename = paste0(getwd(), "/data/BAMM/SamplingCorrected/", "run_info.txt"),
                                  sampleProbsFilename = paste0(getwd(), "/data/BAMM/SamplingCorrected/", "sample_probs.txt"),
                                  eventDataOutfile = paste0(getwd(), "/data/BAMM/SamplingCorrected/", "event_data.txt"),
                                  mcmcOutfile = paste0(getwd(), "/data/BAMM/SamplingCorrected/", "mcmc_out.txt"),
                                  chainSwapFileName = paste0(getwd(), "/data/BAMM/SamplingCorrected/", "chain_swap.txt"),
                                  sampleProbsFilename = paste0(getwd(), "/data/BAMM/", "samplingProbsBats.txt"),
                                  numberOfGenerations = "50000000",
                                  overwrite = "1", 
                                  lambdaInitPrior = as.numeric(BAMMpriors["lambdaInitPrior"]), 
                                  lambdaShiftPrior = as.numeric(BAMMpriors["lambdaShiftPrior"]), 
                                  muInitPrior = as.numeric(BAMMpriors["muInitPrior"]), 
                                  numberOfChains = 39,
                                  expectedNumberOfShifts = "1", # Results are similar for "1", "10" and "25"
                                  minCladeSizeForShift = 1) # Allow shifts to occur in all branches; thus for estimates for terminal branches, in which we are interested.
)

#### Running the BAMM model in the shell ####

# In the shell, run the following code. Note that it will take at least 24 hours under 38 cores and 2.2 GHz clock. 
# I highly recommend the use of 'screen' if done remotely.

# bamm -c /home/pedro.braga/chapter-BiogHistPhyloRelatedSpatScales/chapter-BiogHistPhyloRelatedSpatScales/data/BAMM/SamplingCorrected/chiroptera.divcontrol.txt

#### Extract the output from BAMM files #####

# Reorder the tree by a series of contiguous rows
Chiroptera.FaurSven.tree.ultra.cw <- reorder(Chiroptera.FaurSven.tree.ultra, "cladewise")

### Get the Event data ###
edata <- getEventData(Chiroptera.FaurSven.tree.ultra.cw, 
                      eventdata = "data/BAMM/SamplingCorrected/event_data.txt", 
                      burnin = 0.1)

### Assess MCMC convergence ####
mcmcout <- read.csv("data/BAMM/SamplingCorrected/mcmc_out.txt", 
                    header = TRUE)

plot(mcmcout$logLik ~ mcmcout$generation)

# Discard the first 10% samples as burn-in:

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# Check the effective sample sizes of the log-likelihood and the number of shift
# events present in each sample

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# In general, we want these to be at least 200 (and 200 is on the low side, but
# might be reasonable for very large datasets).

### Analysis of rate shifts ####

# Compute the posterior probabilities of models sampled using BAMM:
post_probs <- table(postburn$N_shifts) / nrow(postburn)

# Which models are part of the set that were sampled?

names(post_probs)

## In general, any model that is not included in names(post_probs) was so lousy
## that it was never even sampled

# You can compute the odds-ratio of two models: post_probs['5'] / post_probs['6']

# Summarize the posterior distribution of the number of shifts
shift_probs <- summary(edata)

# Compute a pairwise matrix of Bayes factors omparing all models with a
# posterior/prior greater than zero

(bfmat <- computeBayesFactors(mcmcout, 
                              expectedNumberOfShifts = 1, 
                              burnin=0.1))

# Visualize the prior and posterior simultaneously
plotPrior(mcmcout, expectedNumberOfShifts=1)

### Visualize mean model-averaged diversificaton rates
# along the tree

bamm.chiroptera <- plot.bammdata(edata, lwd=2, 
                                 labels = F, 
                                 legend = T, 
                                 cex = 0.5)

addBAMMshifts(edata, 
              cex=2)



## Bayesian credible sets of shift configurations

# Identify the 95% credible set of distinct shift configurations
css <- credibleShiftSet(edata, 
                        expectedNumberOfShifts = 1, 
                        threshold = 5, 
                        set.limit = 0.95)

css$number.distinct
summary(css)

# Generate phylorate plots for each of the N shift configurations with the
# highest posterior probabilities
plot.credibleshiftset(css)

dev.off()

### Finding the single best shift configuration ####

# Alternative 1: Obtain the one with the maximum a posteriori (MAP) probability
best <- getBestShiftConfiguration(edata, 
                                  expectedNumberOfShifts=25)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)

# Alternative 2: extract the shift configuration that maximizes the marginal
# probability of rate shifts along individual branches
msc.set <- maximumShiftCredibility(edata, 
                                   maximize='product')

# Pull out a single representative and plot it
msc.config <- subsetEventData(edata, 
                              index = msc.set$sampleindex)
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)

# It is recommended to use credibleShiftSet and not maximumShiftCredibility

## Plotting rate shifts using `plot.phylo`
shiftnodes <- getShiftNodesFromIndex(edata, 
                                     index = 1) # e.g., for the first one

plot.phylo(Chiroptera.FaurSven.tree.ultra.cw, 
           show.tip.label = F)

nodelabels(node = shiftnodes, 
           pch=21, 
           col="red", 
           cex=1.5)

### Rate-through-time analysis ####

# plotRateThroughTime() produces a plot with density shading on confidence
# regions for your speciation-through-time curve.

# This is the RTT plot for speciation rates:
plotRateThroughTime(edata,
                    intervalCol="blue", 
                    avgCol="blue", 
                    ylim=c(0, 0.3))

# RTT plot for extinction rates:
plotRateThroughTime(edata,
                    ratetype = "extinction",
                    intervalCol="red", 
                    avgCol="red", 
                    ylim=c(0, 0.3))

# RTT for net diversification rates
plotRateThroughTime(edata,
                    ratetype = "netdiv",
                    intervalCol="red", 
                    avgCol="blue", 
                    ylim=c(0, 0.3))

### Cumulative shift probabilities ####

# The cumulative shift probability tree shows the cumulative probability, on
# each branch, that a shift occurred somewhere between the focal branch and the
# root of the tree

# The occurrence of such a shift implies that evolutionary dynamics on the focal
# branch are decoupled from the “background” diversification or trait
# evolutionary process at the root of the tree.

cst <- cumulativeShiftProbsTree(edata)

# The marginal shift probability tree is a copy of the target phylogeny, but
# where each branch length is equal to the branch-specific marginal probability
# that a rate-shift occurred on the focal branch. For example, a branch length
# of 0.333 implies that 1/3 of all samples from the posterior had a rate shift
# on the focal branch.

plot.phylo(cst,
           show.tip.label = F)

edgecols <- rep('black', 
                length(Chiroptera.FaurSven.tree.ultra.cw$edge.length))

# And this should plot your tree such that all branches with cumulative
# shift probabilities of 0.95 or higher are identified in red.

is_highprobshift <- cst$edge.length >= 0.95
edgecols[is_highprobshift] <- "red"

plot.phylo(Chiroptera.FaurSven.tree.ultra.cw, 
           edge.color = edgecols,
           show.tip.label = F)

# computing the marginal shift probs tree:
mst <- marginalShiftProbsTree(edata)

#compare the two types of shift trees side-by-side:
plot.new()

par(mfrow=c(1,3))

plot.phylo(mst, 
           no.margin=TRUE, 
           show.tip.label=FALSE)

plot.phylo(cst, no.margin=TRUE, 
           show.tip.label=FALSE)

plot.phylo(Chiroptera.FaurSven.tree.ultra.cw, 
           edge.color = edgecols,
           show.tip.label = F)

dev.off()

####  Estimate individual tip-specific evolutionary rates #### 

## Obtain speciation, extinction, net diversification for the tips of the tree
## from all posterior samples in the output
tip.rates <- getTipRates(edata)

str(tip.rates)

# Should we subset the event data so tip rates are obtained only from the
# credible subset?