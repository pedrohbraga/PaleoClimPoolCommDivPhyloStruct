########################
### Global sampling ####
########################

# Calculate mntd values for each community
SES.MNTD.Chiroptera.World <- sf.ses.mntd(Chiroptera.FaurSven.comm,
                                       cophenetic(Chiroptera.FaurSven.tree), # Convert phylogeny to distance matrix
                                       abundance.weighted = FALSE, 
                                       runs = 999, cores = 24)


MNTD.LatLong.Env.World <- data.frame(cbind(regions.LatLong, 
                                          Env,
                                          SES.MNTD.Chiroptera.World,
                                          SamplingPool = rep("Global sampling", nrow(regions.LatLong))))

##########################################
#### New World and Old World Sampling ####
##########################################

NW.LatLong <- subset(regions.LatLong, ID_Realm == "Nearctic" | ID_Realm == "Neotropical")
OW.LatLong <- subset(regions.LatLong, !(ID_Realm %in% c("Nearctic", "Neotropical")))
nrow(NW.LatLong) + nrow(OW.LatLong)

NW.Comm <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(NW.LatLong), ]
NW.Comm <- NW.Comm[ , colSums(NW.Comm) >= 1]

OW.Comm <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(OW.LatLong), ]
OW.Comm <- OW.Comm[ , colSums(OW.Comm) >= 1]

dim(NW.Comm) ; nrow(NW.LatLong)
dim(OW.Comm) ; nrow(OW.LatLong)

### Subset the datasets

Chiroptera.NW.Comm <- match.phylo.comm(Chiroptera.FaurSven.tree, NW.Comm)$comm
Chiroptera.NW.Tree <- match.phylo.comm(Chiroptera.FaurSven.tree, NW.Comm)$phy

Chiroptera.OW.Comm <- match.phylo.comm(Chiroptera.FaurSven.tree, OW.Comm)$comm
Chiroptera.OW.Tree <- match.phylo.comm(Chiroptera.FaurSven.tree, OW.Comm)$phy

#############################################################

#### New World ####
# Calculate mntd values for each community
SES.MNTD.Chiroptera.NW <- sf.ses.mntd(Chiroptera.NW.Comm,
                                    cophenetic(Chiroptera.NW.Tree), # Convert phylogeny to distance matrix
                                    abundance.weighted = FALSE, 
                                    runs = 999, cores = 24)


#### Old World ####

SES.MNTD.Chiroptera.OW <- sf.ses.mntd(Chiroptera.OW.Comm,
                                    cophenetic(Chiroptera.OW.Tree), # Convert phylogeny to distance matrix
                                    abundance.weighted = FALSE, 
                                    runs = 999, cores = 24)

MNTD.LatLong.Env.NW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.NW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.NW.Comm)), ],
                                       SES.MNTD.Chiroptera.NW,
                                       SamplingPool = rep("Hemispheric sampling", nrow(SES.MNTD.Chiroptera.NW))))
dim(MNTD.LatLong.Env.NW)

MNTD.LatLong.Env.OW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.OW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.OW.Comm)), ],
                                       SES.MNTD.Chiroptera.OW,
                                       SamplingPool = rep("Hemispheric sampling", nrow(SES.MNTD.Chiroptera.OW))))
dim(MNTD.LatLong.Env.OW)

MNTD.LatLong.Env.NW.OW <- rbind(MNTD.LatLong.Env.NW, MNTD.LatLong.Env.OW)

dim(MNTD.LatLong.Env.NW.OW)

##########################################
############# Realm Sampling #############
##########################################

# Create a list
SubRegion.Realm.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Realm))) {
  
  SubRegion.Realm.LatLong <- subset(regions.LatLong, ID_Realm == toString(i))
  SubRegion.Realm.Comm. <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(SubRegion.Realm.LatLong), ]
  SubRegion.Realm.Comm. <- SubRegion.Realm.Comm.[ , colSums(SubRegion.Realm.Comm.) >= 1]
  
  ### Subset the datasets
  SubRegion.Realm.CommPhylo <- match.phylo.comm(Chiroptera.FaurSven.tree, SubRegion.Realm.Comm.)
  SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                  Phylo = SubRegion.Realm.CommPhylo$phy)
}

SubRegion.Realm.Comm.Phylo["NA"] <- NULL

MNTD.LatLong.Env.Realm = NULL

for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mntd values for each community
  SES.MNTD.Chiroptera.Realm <- sf.ses.mntd(SubRegion.Realm.Comm.Phylo[[j]]$Comm,
                                         cophenetic(SubRegion.Realm.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                         abundance.weighted = FALSE, 
                                         runs = 999, cores = 24)
  
  MNTD.LatLong.Env.SubRegion.Realm <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ],
                                                      SES.MNTD.Chiroptera.Realm,
                                                      SamplingPool = rep("Realm sampling", nrow(SES.MNTD.Chiroptera.Realm))))
  
  MNTD.LatLong.Env.Realm <- rbind(MNTD.LatLong.Env.Realm, MNTD.LatLong.Env.SubRegion.Realm)
}

#######################################################
############# Within Realm Biome Sampling #############
#######################################################

# Create a list
SubRegion.Biome.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Biome_Realm))) {
  
  SubRegion.Biome.LatLong <- subset(regions.LatLong, ID_Biome_Realm == toString(i))
  SubRegion.Biome.Comm. <- Chiroptera.FaurSven.comm[rownames(Chiroptera.FaurSven.comm) %in% row.names(SubRegion.Biome.LatLong), ]
  SubRegion.Biome.Comm. <- SubRegion.Biome.Comm.[ , colSums(SubRegion.Biome.Comm.) >= 1]
  
  ### Subset the datasets
  SubRegion.Biome.CommPhylo <- match.phylo.comm(Chiroptera.FaurSven.tree, SubRegion.Biome.Comm.)
  SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                  Phylo = SubRegion.Biome.CommPhylo$phy)
}

SubRegion.Biome.Comm.Phylo["NA"] <- NULL

MNTD.LatLong.Env.Biome = NULL

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate standardized versions under the uniform model
  SES.MNTD.Chiroptera.Biome <-   sf.ses.mntd(SubRegion.Biome.Comm.Phylo[[j]]$Comm,
                                           cophenetic(SubRegion.Biome.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                           abundance.weighted = FALSE, 
                                           runs = 999, cores = 24)
  
  MNTD.LatLong.Env.SubRegion.Biome <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ],
                                                      SES.MNTD.Chiroptera.Biome,
                                                      SamplingPool = rep("Biome sampling", nrow(SES.MNTD.Chiroptera.Biome))))
  
  MNTD.LatLong.Env.Biome <- rbind(MNTD.LatLong.Env.Biome, MNTD.LatLong.Env.SubRegion.Biome)
}

#################################
#### Combine all data frames ####
#################################

MNTD.LatLong.Env.AllScales <- rbind(MNTD.LatLong.Env.World,
                                   MNTD.LatLong.Env.NW.OW,
                                   MNTD.LatLong.Env.Realm,
                                   MNTD.LatLong.Env.Biome)

head(MNTD.LatLong.Env.AllScales)


# Calculate NRI
MNTD.LatLong.Env.AllScales$NRI <- - 1 * ((MNTD.LatLong.Env.AllScales$mntd.obs - MNTD.LatLong.Env.AllScales$mntd.rand.mean) / MNTD.LatLong.Env.AllScales$mntd.rand.sd)

# Rename dataset

MNTD.LatLong.Env.AllScales$ID_Realm. <- factor(MNTD.LatLong.Env.AllScales$ID_Realm, levels=c('Neotropical','Nearctic','Afrotropical','Palearctic', 'Indomalay', 'Australasian'))

# Save dataset

write.table(MNTD.LatLong.Env.AllScales, "MNTD.LatLong.Env.AllScales.FaurSven.txt")
