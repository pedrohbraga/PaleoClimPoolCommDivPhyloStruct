######################################################################################
#### Code to apply a null-model framework to compute the phylogenetic structure   ####
#### of communities across a gradient of restrictive sampling pools               ####
#                                                                                    #
# Description: To assess how patterns of phylogenetic community structure change     #
# as a function of spatial scale, we applied commonly used null-models to            #
# estimate standardized effect sizes for both metrics (MPDSES and MNTDSES,           #
# respectively) (Kembel et al., 2010). Each null-model simulated random              #
# assemblages by permuting species names across the phylogeny tips 999 times for     #
# a given species pool (i.e., the sampling pool in which species were sampled to     #
# compose random assemblages): global (all species in the phylogeny), east-west      #
# hemispherical (i.e., Old World and New World species pools), biogeographical       #
# realms, tectonic plates, within-realm biomes, and within-realm                     #
# terrestrial-ecoregion scales. NRI and NTI were obtained by multiplying minus       #
# one to MPDSES and MNTDSES, respectively.                                           #
#                                                                                    #
# The framework used here uses the sf.ses.mntd() function, which modifies the        #
# picante::ses.mntd() function to allow for parallel computation using SNOW          #
# (snowfall).                                                                        #
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2019-01-28"                                                          #
#                                                                                    # 
######################################################################################

# Begin of the computation of SES.MNTD values across the sampling pool
# restriction gradient

########################
### Global sampling ####
########################

# Calculate mntd values for each community
SES.MNTD.Chiroptera.World <- sf.ses.mntd(Chiroptera.FaurSven.comm,
                                         cophenetic(Chiroptera.FaurSven.tree), # phylogenetic distance matrix
                                         abundance.weighted = FALSE, 
                                         runs = 999, cores = 40)

# Add everything to a data.frame
MNTD.LatLong.Env.World <- data.frame(cbind(regions.LatLong,
                                           current.Env,
                                           SES.MNTD.Chiroptera.World,
                                           SamplingPool = rep("Global sampling", nrow(regions.LatLong))))

##########################################
#### New World and Old World sampling ####
##########################################

# Defining the New world and the Old World
NW.LatLong <- subset(regions.LatLong, 
                     ID_Realm == "Nearctic" | ID_Realm == "Neotropical")

OW.LatLong <- subset(regions.LatLong, 
                     !(ID_Realm %in% c("Nearctic", "Neotropical")))

nrow(NW.LatLong) + nrow(OW.LatLong); nrow(regions.LatLong) # Checking the number of rows

# Extracting NW communities
NW.Comm <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(NW.LatLong), ]
NW.Comm <- NW.Comm[ , colSums(NW.Comm) >= 1]

# Extracting OW communities
OW.Comm <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(OW.LatLong), ]
OW.Comm <- OW.Comm[ , colSums(OW.Comm) >= 1]

dim(NW.Comm); nrow(NW.LatLong)
dim(OW.Comm); nrow(OW.LatLong)

### Subset the community and phylogenetic data for each hemispheric region
Chiroptera.NW.Comm <- match.phylo.comm(Chiroptera.FaurSven.tree, NW.Comm)$comm
Chiroptera.NW.Comm <- Chiroptera.NW.Comm[ , colSums(Chiroptera.NW.Comm) >= 1]
Chiroptera.NW.Tree <- match.phylo.comm(Chiroptera.FaurSven.tree, NW.Comm)$phy

Chiroptera.OW.Comm <- match.phylo.comm(Chiroptera.FaurSven.tree, OW.Comm)$comm
Chiroptera.OW.Comm <- Chiroptera.OW.Comm[ , colSums(Chiroptera.NW.Comm) >= 1]
Chiroptera.OW.Tree <- match.phylo.comm(Chiroptera.FaurSven.tree, OW.Comm)$phy


# Calculate the ses.mntd for the communities accounting for the hemispheric sampling
# pool

#### New World ####
SES.MNTD.Chiroptera.NW <- sf.ses.mntd(Chiroptera.NW.Comm,
                                      cophenetic(Chiroptera.NW.Tree), # Convert phylogeny to distance matrix
                                      abundance.weighted = FALSE, 
                                      runs = 999, cores = 40)

#### Old World ####
SES.MNTD.Chiroptera.OW <- sf.ses.mntd(Chiroptera.OW.Comm,
                                      cophenetic(Chiroptera.OW.Tree), # Convert phylogeny to distance matrix
                                      abundance.weighted = FALSE, 
                                      runs = 999, cores = 40)

# Export results to data.frames
MNTD.LatLong.Env.NW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.NW.Comm)), ], 
                                        current.Env[(rownames(current.Env) %in% row.names(Chiroptera.NW.Comm)), ],
                                        SES.MNTD.Chiroptera.NW,
                                        SamplingPool = rep("Hemispheric sampling", nrow(SES.MNTD.Chiroptera.NW))))

MNTD.LatLong.Env.OW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.OW.Comm)), ], 
                                        current.Env[(rownames(current.Env) %in% row.names(Chiroptera.OW.Comm)), ],
                                        SES.MNTD.Chiroptera.OW,
                                        SamplingPool = rep("Hemispheric sampling", nrow(SES.MNTD.Chiroptera.OW))))

# Bind both hemispheres in one data.frame
MNTD.LatLong.Env.NW.OW <- rbind(MNTD.LatLong.Env.NW, MNTD.LatLong.Env.OW)

# dim(MNTD.LatLong.Env.NW.OW)

###########################
### Realm sampling pool ###  
###########################

# Create a list
SubRegion.Realm.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Realm))) {
  
  SubRegion.Realm.LatLong <- subset(regions.LatLong, ID_Realm == toString(i))
  SubRegion.Realm.Comm. <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(SubRegion.Realm.LatLong), ]
  SubRegion.Realm.Comm. <- SubRegion.Realm.Comm.[ , colSums(SubRegion.Realm.Comm.) >= 1]
  
  ### Subset the data sets
  SubRegion.Realm.CommPhylo <- match.phylo.comm(Chiroptera.FaurSven.tree, SubRegion.Realm.Comm.)
  SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                  Phylo = SubRegion.Realm.CommPhylo$phy)
}

SubRegion.Realm.Comm.Phylo["NA"] <- NULL

MNTD.LatLong.Env.Realm <- NULL

for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
  # Calculate mntd values for each community
  SES.MNTD.Chiroptera.Realm <- sf.ses.mntd(SubRegion.Realm.Comm.Phylo[[j]]$Comm,
                                           cophenetic(SubRegion.Realm.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                           abundance.weighted = FALSE, 
                                           runs = 999, cores = 40)
  
  MNTD.LatLong.Env.SubRegion.Realm <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ], 
                                                       current.Env[(rownames(current.Env) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ],
                                                       SES.MNTD.Chiroptera.Realm,
                                                       SamplingPool = rep("Realm sampling", nrow(SES.MNTD.Chiroptera.Realm))))
  
  MNTD.LatLong.Env.Realm <- rbind(MNTD.LatLong.Env.Realm, MNTD.LatLong.Env.SubRegion.Realm)
}

########################
#### Plate Sampling ####
########################

# Create a list
SubRegion.Plate.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_PlateName))) {
  SubRegion.Plate.LatLong <- subset(regions.LatLong, ID_PlateName == toString(i))
  SubRegion.Plate.Comm. <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(SubRegion.Plate.LatLong), , drop = FALSE]
  SubRegion.Plate.Comm. <- SubRegion.Plate.Comm.[ , colSums(SubRegion.Plate.Comm.) > 1, drop = FALSE]
  
  ### Subset the datasets
  SubRegion.Plate.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree, SubRegion.Plate.Comm.)
  if(is.null(SubRegion.Plate.CommPhylo$phy) == FALSE){
    SubRegion.Plate.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Plate.CommPhylo$comm,
                                                    Phylo = SubRegion.Plate.CommPhylo$phy)
  }
}

SubRegion.Plate.Comm.Phylo["NA"] <- NULL

MNTD.LatLong.Env.Plate <- NULL

# setdiff(as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"]),
#        unique(MNTD.LatLong.Env.Plate$ID_PlateName))

for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){
  # Calculate mntd values for each community
  SES.MNTD.Chiroptera.Plate <- sf.ses.mntd(SubRegion.Plate.Comm.Phylo[[j]]$Comm,
                                           cophenetic(SubRegion.Plate.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                           abundance.weighted = FALSE, 
                                           runs = 999, cores = 40)
  
  MNTD.LatLong.Env.SubRegion.Plate <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Plate.Comm.Phylo[[j]]$Comm)), ], 
                                                       current.Env[(rownames(current.Env) %in% row.names(SubRegion.Plate.Comm.Phylo[[j]]$Comm)), ],
                                                       SES.MNTD.Chiroptera.Plate,
                                                       SamplingPool = rep("Plate sampling", nrow(SES.MNTD.Chiroptera.Plate))))
  
  MNTD.LatLong.Env.Plate <- rbind(MNTD.LatLong.Env.Plate, MNTD.LatLong.Env.SubRegion.Plate)
}

##########################################
#### Within-Realm Biome Sampling pool ####
##########################################

# Create a list
SubRegion.Biome.Comm.Phylo <- list()

for(i in droplevels(unique(as.factor(regions.LatLong$ID_Biome_Realm)))) {
  
  SubRegion.Biome.LatLong <- subset(regions.LatLong, ID_Biome_Realm == toString(i))
  SubRegion.Biome.Comm. <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(SubRegion.Biome.LatLong), ]
  SubRegion.Biome.Comm. <- SubRegion.Biome.Comm.[ , colSums(SubRegion.Biome.Comm.) >= 1, drop = FALSE]
  
  ### Subset the datasets
  SubRegion.Biome.CommPhylo <- match.phylo.comm(Chiroptera.FaurSven.tree, SubRegion.Biome.Comm.)
  if(is.null(SubRegion.Biome.CommPhylo$phy) == FALSE){
    SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                    Phylo = SubRegion.Biome.CommPhylo$phy)
  }
}

# ls(SubRegion.Biome.Comm.Phylo)[grepl("NA", ls(SubRegion.Biome.Comm.Phylo))]
# SubRegion.Biome.Comm.Phylo["NA__NA"] <- NULL

SubRegion.Biome.Comm.Phylo[ls(SubRegion.Biome.Comm.Phylo)[grepl("NA", 
                                                                ls(SubRegion.Biome.Comm.Phylo))]] <- NULL
SubRegion.Biome.Comm.Phylo$Tropical_Subtropical_Grasslands_Savannas_Shrublands__Oceanic <- NULL

# Applying is.null function to a list of lists
# lapply(SubRegion.Biome.Comm.Phylo, sapply, is.null)
# ls(SubRegion.Biome.Comm.Phylo)

MNTD.LatLong.Env.Biome = NULL

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate standardized versions under the uniform model
  SES.MNTD.Chiroptera.Biome <-   sf.ses.mntd(SubRegion.Biome.Comm.Phylo[[j]]$Comm,
                                             cophenetic(SubRegion.Biome.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                             abundance.weighted = FALSE, 
                                             runs = 999, cores = 40)
  
  MNTD.LatLong.Env.SubRegion.Biome <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ], 
                                                       current.Env[(rownames(current.Env) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ],
                                                       SES.MNTD.Chiroptera.Biome,
                                                       SamplingPool = rep("Biome sampling", 
                                                                          nrow(SES.MNTD.Chiroptera.Biome))))
  
  MNTD.LatLong.Env.Biome <- rbind(MNTD.LatLong.Env.Biome, MNTD.LatLong.Env.SubRegion.Biome)
}


##########################################
############# Ecoregion Sampling #########
##########################################

# Create a list
SubRegion.Ecoregion.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$NAME_Ecoregion))) {
  SubRegion.Ecoregion.LatLong <- subset(regions.LatLong, NAME_Ecoregion == toString(i))
  SubRegion.Ecoregion.Comm. <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(SubRegion.Ecoregion.LatLong), , drop = FALSE]
  SubRegion.Ecoregion.Comm. <- SubRegion.Ecoregion.Comm.[ , colSums(SubRegion.Ecoregion.Comm.) > 1, drop = FALSE]
  
  ### Subset the datasets
  SubRegion.Ecoregion.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree,
                                                     SubRegion.Ecoregion.Comm.)
  
  if(is.null(SubRegion.Ecoregion.CommPhylo$phy) == FALSE){
    if(ncol(SubRegion.Ecoregion.CommPhylo$comm) > 1){
      SubRegion.Ecoregion.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Ecoregion.CommPhylo$comm,
                                                          Phylo = SubRegion.Ecoregion.CommPhylo$phy)
    }
  }
}

SubRegion.Ecoregion.Comm.Phylo["NA"] <- NULL

MNTD.LatLong.Env.Ecoregion <- NULL

for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){
  # Calculate MNTD values for each community
  SES.MNTD.Chiroptera.Ecoregion <- sf.ses.mntd(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm,
                                               cophenetic(SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                               abundance.weighted = FALSE, 
                                               runs = 999, cores = 40)
  
  MNTD.LatLong.Env.SubRegion.Ecoregion <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm)), ], 
                                                           current.Env[(rownames(current.Env) %in% row.names(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm)), ],
                                                           SES.MNTD.Chiroptera.Ecoregion,
                                                           SamplingPool = rep("Ecoregion sampling", 
                                                                              nrow(SES.MNTD.Chiroptera.Ecoregion))))
  
  MNTD.LatLong.Env.Ecoregion <- rbind(MNTD.LatLong.Env.Ecoregion, MNTD.LatLong.Env.SubRegion.Ecoregion)
}


# End of the computation of SES.MNTD values across the sampling pool
# restriction gradient

#################################
#### Combine all data frames ####
#################################

MNTD.LatLong.Env.AllScales <- rbind(MNTD.LatLong.Env.World,
                                    MNTD.LatLong.Env.NW.OW,
                                    MNTD.LatLong.Env.Realm,
                                    MNTD.LatLong.Env.Plate,
                                    MNTD.LatLong.Env.Biome,
                                    MNTD.LatLong.Env.Ecoregion)

head(MNTD.LatLong.Env.AllScales)

# Calculate NTI by multiplying minus one to MNTD.SES

MNTD.LatLong.Env.AllScales$NTI <- - 1 * ((MNTD.LatLong.Env.AllScales$mntd.obs - MNTD.LatLong.Env.AllScales$mntd.rand.mean) / MNTD.LatLong.Env.AllScales$mntd.rand.sd)

MNTD.LatLong.Env.AllScales$PhyloMetricScale <- rep("MNTD", 
                                                   nrow(MNTD.LatLong.Env.AllScales)) 

# Refactor data set

MNTD.LatLong.Env.AllScales$ID_Realm <- factor(MNTD.LatLong.Env.AllScales$ID_Realm, 
                                              levels = c('Neotropical',
                                                         'Nearctic',
                                                         'Afrotropical',
                                                         'Palearctic', 
                                                         'Indomalay', 
                                                         'Australasian'))

MNTD.LatLong.Env.AllScales$SamplingPool <-  factor(MNTD.LatLong.Env.AllScales$SamplingPool,  
                                                   levels = c("Global sampling",
                                                              "Hemispheric sampling",
                                                              "Realm sampling",
                                                              "Plate sampling",
                                                              "Biome sampling",
                                                              "Ecoregion sampling"))
# Save data set

write.table(MNTD.LatLong.Env.AllScales, 
            "MNTD.LatLong.Env.AllScales.50K.FaurSven.3.txt")

# Read data set (if needed)

# MNTD.LatLong.Env.AllScales <- read.table("MNTD.LatLong.Env.AllScales.50K.FaurSven.3.txt", h = T)

