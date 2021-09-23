########################
### Global sampling ####
########################

# Calculate mpd values for each community
SES.MPD.Chiroptera.World <- sf.ses.mpd(Chiroptera.Jones.comm,
                                       cophenetic(Chiroptera.Jones.tree), # Convert phylogeny to distance matrix
                                       abundance.weighted = FALSE, 
                                       runs = 999, cores = 24)


MPD.LatLong.Env.World <- data.frame(cbind(regions.LatLong, 
                                          Env,
                                          SES.MPD.Chiroptera.World,
                                          SamplingPool = rep("Global sampling", nrow(regions.LatLong))))

##########################################
#### New World and Old World Sampling ####
##########################################

NW.LatLong <- subset(regions.LatLong, ID_Realm == "Nearctic" | ID_Realm == "Neotropical")
OW.LatLong <- subset(regions.LatLong, !(ID_Realm %in% c("Nearctic", "Neotropical")))
nrow(NW.LatLong) + nrow(OW.LatLong)

NW.Comm <- Chiroptera.Jones.comm[rownames(Chiroptera.Jones.comm) %in% row.names(NW.LatLong), ]
OW.Comm <- Chiroptera.Jones.comm[rownames(Chiroptera.Jones.comm) %in% row.names(OW.LatLong), ]

dim(NW.Comm) ; nrow(NW.LatLong)
dim(OW.Comm) ; nrow(OW.LatLong)

### Subset the datasets

Chiroptera.NW.Comm <- match.phylo.comm(Chiroptera.Jones.tree, NW.Comm)$comm
Chiroptera.NW.Tree <- match.phylo.comm(Chiroptera.Jones.tree, NW.Comm)$phy

Chiroptera.OW.Comm <- match.phylo.comm(Chiroptera.Jones.tree, OW.Comm)$comm
Chiroptera.OW.Tree <- match.phylo.comm(Chiroptera.Jones.tree, OW.Comm)$phy

#############################################################

#### New World ####
# Calculate mpd values for each community
SES.MPD.Chiroptera.NW <- sf.ses.mpd(Chiroptera.NW.Comm,
                                    cophenetic(Chiroptera.NW.Tree), # Convert phylogeny to distance matrix
                                    abundance.weighted = FALSE, 
                                    runs = 999, cores = 24)


#### Old World ####

SES.MPD.Chiroptera.OW <- sf.ses.mpd(Chiroptera.OW.Comm,
                                    cophenetic(Chiroptera.OW.Tree), # Convert phylogeny to distance matrix
                                    abundance.weighted = FALSE, 
                                    runs = 999, cores = 24)

MPD.LatLong.Env.NW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.NW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.NW.Comm)), ],
                                       SES.MPD.Chiroptera.NW,
                                       SamplingPool = rep("Hemispheric sampling", nrow(SES.MPD.Chiroptera.NW))))
dim(MPD.LatLong.Env.NW)

MPD.LatLong.Env.OW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.OW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.OW.Comm)), ],
                                       SES.MPD.Chiroptera.OW,
                                       SamplingPool = rep("Hemispheric sampling", nrow(SES.MPD.Chiroptera.OW))))
dim(MPD.LatLong.Env.OW)

MPD.LatLong.Env.NW.OW <- rbind(MPD.LatLong.Env.NW, MPD.LatLong.Env.OW)

dim(MPD.LatLong.Env.NW.OW)

##########################################
############# Realm Sampling #############
##########################################

# Create a list
SubRegion.Realm.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Realm))) {
  
  SubRegion.Realm.LatLong <- subset(regions.LatLong, ID_Realm == toString(i))
  SubRegion.Realm.Comm. <- Chiroptera.Jones.comm[rownames(Chiroptera.Jones.comm) %in% row.names(SubRegion.Realm.LatLong), ]
  
  ### Subset the datasets
  SubRegion.Realm.CommPhylo <- match.phylo.comm(Chiroptera.Jones.tree, SubRegion.Realm.Comm.)
  SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                  Phylo = SubRegion.Realm.CommPhylo$phy)
}

MPD.LatLong.Env.Realm = NULL

for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mpd values for each community
  SES.MPD.Chiroptera.Realm <- sf.ses.mpd(SubRegion.Realm.Comm.Phylo[[j]]$Comm,
                                         cophenetic(SubRegion.Realm.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                         abundance.weighted = FALSE, 
                                         runs = 999, cores = 24)
  
  MPD.LatLong.Env.SubRegion.Realm <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ],
                                                      SES.MPD.Chiroptera.Realm,
                                                      SamplingPool = rep("Realm sampling", nrow(SES.MPD.Chiroptera.Realm))))
  
  MPD.LatLong.Env.Realm <- rbind(MPD.LatLong.Env.Realm, MPD.LatLong.Env.SubRegion.Realm)
}

##########################################
############# Biome Sampling #############
##########################################

# Create a list
SubRegion.Biome.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Biome))) {
  
  SubRegion.Biome.LatLong <- subset(regions.LatLong, ID_Biome == toString(i))
  SubRegion.Biome.Comm. <- Chiroptera.Jones.comm[rownames(Chiroptera.Jones.comm) %in% row.names(SubRegion.Biome.LatLong), ]
  
  ### Subset the datasets
  SubRegion.Biome.CommPhylo <- match.phylo.comm(Chiroptera.Jones.tree, SubRegion.Biome.Comm.)
  SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                  Phylo = SubRegion.Biome.CommPhylo$phy)
}

MPD.LatLong.Env.Biome = NULL

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate standardized versions under the uniform model
  SES.MPD.Chiroptera.Biome <-   sf.ses.mpd(SubRegion.Biome.Comm.Phylo[[j]]$Comm,
                                           cophenetic(SubRegion.Biome.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                           abundance.weighted = FALSE, 
                                           runs = 999, cores = 24)
  
  MPD.LatLong.Env.SubRegion.Biome <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ],
                                                      SES.MPD.Chiroptera.Biome,
                                                      SamplingPool = rep("Biome sampling", nrow(SES.MPD.Chiroptera.Biome))))
  
  MPD.LatLong.Env.Biome <- rbind(MPD.LatLong.Env.Biome, MPD.LatLong.Env.SubRegion.Biome)
}


#################################
#### Combine all data frames ####
#################################

MPD.LatLong.Env.AllScales <- rbind(MPD.LatLong.Env.World,
                                   MPD.LatLong.Env.NW.OW,
                                   MPD.LatLong.Env.Realm,
                                   MPD.LatLong.Env.Biome)

head(MPD.LatLong.Env.AllScales)

# Calculate NRI
MPD.LatLong.Env.AllScales$NRI <- - 1 * ((MPD.LatLong.Env.AllScales$mpd.obs - MPD.LatLong.Env.AllScales$mpd.rand.mean) / MPD.LatLong.Env.AllScales$mpd.rand.sd)

# Rename dataset
MPD.LatLong.Env.AllScales$ID_Realm. <- factor(MPD.LatLong.Env.AllScales$ID_Realm, levels=c('Neotropical','Nearctic','Afrotropical','Palearctic', 'Indomalay', 'Australasian'))

# Save dataset
write.table(MPD.LatLong.Env.AllScales, "MPD.LatLong.Env.AllScales.txt")