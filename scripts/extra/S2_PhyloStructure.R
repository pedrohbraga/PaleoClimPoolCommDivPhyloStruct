length(Chiroptera.FaurbySven.Phylo.ultra$tip.label); ncol(Chiroptera.FaurbySven.comm)

########################
### Global sampling ####
########################

# Calculate mpd values for each community
MPD.Chiroptera.World <- mpd.query(Chiroptera.FaurbySven.Phylo.ultra, Chiroptera.FaurbySven.comm)

# Calculate standardized versions under the uniform model
SES.MPD.Chiroptera.World <- mpd.query(Chiroptera.FaurbySven.Phylo.ultra, Chiroptera.FaurbySven.comm, 
                                      standardize = TRUE)

# Create random abundance weights
weights = runif(length(Chiroptera.FaurbySven.Phylo.ultra$tip.label))
names(weights) = Chiroptera.FaurbySven.Phylo.ultra$tip.label

#Use query function to calculate standardized versions under the sequential model
SES.MPD.Seq.Chiroptera.World <- mpd.query(Chiroptera.FaurbySven.Phylo.ultra, Chiroptera.FaurbySven.comm, 
                                          standardize = TRUE, 
                                          null.model= "sequential",
                                          abundance.weights = weights, 
                                          reps=1000)

MPD.LatLong.Env.World <- data.frame(cbind(regions.LatLong, 
                                          Env,
                                          MPD = MPD.Chiroptera.World,
                                          SES.MPD.Seq = SES.MPD.Seq.Chiroptera.World,
                                          SES.MPD = SES.MPD.Chiroptera.World,
                                          SamplingPool = rep("Global sampling", nrow(regions.LatLong))))

head(MPD.LatLong.Env.World)


##########################################
#### New World and Old World Sampling ####
##########################################

NW.LatLong <- subset(regions.LatLong, ID_Realm == "Nearctic" | ID_Realm == "Neotropical")
OW.LatLong <- subset(regions.LatLong, !(ID_Realm %in% c("Nearctic", "Neotropical")))
nrow(NW.LatLong) + nrow(OW.LatLong)

NW.Comm <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(NW.LatLong), ]
OW.Comm <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(OW.LatLong), ]

dim(NW.Comm) ; nrow(NW.LatLong)
dim(OW.Comm) ; nrow(OW.LatLong)

### Subset the datasets

Chiroptera.NW.Comm <- match.phylo.comm(Mammal.FaurbySven.tree.1, NW.Comm)$comm
Chiroptera.NW.Tree <- match.phylo.comm(Mammal.FaurbySven.tree.1, NW.Comm)$phy

Chiroptera.OW.Comm <- match.phylo.comm(Mammal.FaurbySven.tree.1, OW.Comm)$comm
Chiroptera.OW.Tree <- match.phylo.comm(Mammal.FaurbySven.tree.1, OW.Comm)$phy

#############################################################

#### New World ####
# Calculate mpd values for each community
MPD.Chiroptera.NW <- mpd.query(Chiroptera.NW.Tree, Chiroptera.NW.Comm)

# Calculate standardized versions under the uniform model
SES.MPD.Chiroptera.NW <- mpd.query(Chiroptera.NW.Tree, Chiroptera.NW.Comm, 
                                      standardize = TRUE)

# Create random abundance weights
weights = runif(length(Chiroptera.NW.Tree$tip.label))
names(weights) = Chiroptera.NW.Tree$tip.label

#Use query function to calculate standardized versions under the sequential model
SES.MPD.Seq.Chiroptera.NW <- mpd.query(Chiroptera.NW.Tree, Chiroptera.NW.Comm, 
                                          standardize = TRUE, 
                                          null.model= "sequential",
                                          abundance.weights = weights, 
                                          reps=1000)

#### Old World ####

# Calculate mpd values for each community
MPD.Chiroptera.OW <- mpd.query(Chiroptera.OW.Tree, Chiroptera.OW.Comm)

# Calculate standardized versions under the uniform model
SES.MPD.Chiroptera.OW <- mpd.query(Chiroptera.OW.Tree, Chiroptera.OW.Comm, 
                                   standardize = TRUE)

# Create random abundance weights
weights = runif(length(Chiroptera.OW.Tree$tip.label))
names(weights) = Chiroptera.OW.Tree$tip.label

#Use query function to calculate standardized versions under the sequential model
SES.MPD.Seq.Chiroptera.OW <- mpd.query(Chiroptera.OW.Tree, Chiroptera.OW.Comm, 
                                       standardize = TRUE, 
                                       null.model= "sequential",
                                       abundance.weights = weights, 
                                       reps=1000)

MPD.LatLong.Env.NW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.NW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.NW.Comm)), ],
                                       MPD = MPD.Chiroptera.NW,
                                       SES.MPD.Seq = SES.MPD.Seq.Chiroptera.NW,
                                       SES.MPD = SES.MPD.Chiroptera.NW,
                                       SamplingPool = rep("Hemispheric sampling", length(MPD.Chiroptera.NW))))
dim(MPD.LatLong.Env.NW)

MPD.LatLong.Env.OW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.OW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.OW.Comm)), ],
                                       MPD = MPD.Chiroptera.OW,
                                       SES.MPD.Seq = SES.MPD.Seq.Chiroptera.OW,
                                       SES.MPD = SES.MPD.Chiroptera.OW,
                                       SamplingPool = rep("Hemispheric sampling", length(MPD.Chiroptera.OW))))
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
  SubRegion.Realm.Comm. <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(SubRegion.Realm.LatLong), ]
  
  ### Subset the datasets
  SubRegion.Realm.CommPhylo <- match.phylo.comm(Mammal.FaurbySven.tree.1, SubRegion.Realm.Comm.)
  SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                  Phylo = SubRegion.Realm.CommPhylo$phy)
}

MPD.LatLong.Env.Realm = NULL

for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mpd values for each community
  MPD.Chiroptera.Realm <- mpd.query(SubRegion.Realm.Comm.Phylo[[j]]$Phylo, SubRegion.Realm.Comm.Phylo[[j]]$Comm)
  
  # Calculate standardized versions under the uniform model
  SES.MPD.Chiroptera.Realm <- mpd.query(SubRegion.Realm.Comm.Phylo[[j]]$Phylo, SubRegion.Realm.Comm.Phylo[[j]]$Comm, 
                                        standardize = TRUE)
  
  # Create random abundance weights
  weights = runif(length(SubRegion.Realm.Comm.Phylo[[j]]$Phylo$tip.label))
  names(weights) = SubRegion.Realm.Comm.Phylo[[j]]$Phylo$tip.label
  
  #Use query function to calculate standardized versions under the sequential model
  SES.MPD.Seq.Chiroptera.Realm <- mpd.query(SubRegion.Realm.Comm.Phylo[[j]]$Phylo, SubRegion.Realm.Comm.Phylo[[j]]$Comm, 
                                            standardize = TRUE, 
                                            null.model= "sequential",
                                            abundance.weights = weights, 
                                            reps = 1000)
  
  MPD.LatLong.Env.SubRegion.Realm <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ],
                                                      MPD = MPD.Chiroptera.Realm,
                                                      SES.MPD.Seq = SES.MPD.Seq.Chiroptera.Realm,
                                                      SES.MPD = SES.MPD.Chiroptera.Realm,
                                                      SamplingPool = rep("Realm sampling", length(MPD.Chiroptera.Realm))))
  
  MPD.LatLong.Env.Realm <- rbind(MPD.LatLong.Env.Realm, MPD.LatLong.Env.SubRegion.Realm)
}

##########################################
############# Biome Sampling #############
##########################################

# Create a list
SubRegion.Biome.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Biome))) {
  
  SubRegion.Biome.LatLong <- subset(regions.LatLong, ID_Biome == toString(i))
  SubRegion.Biome.Comm. <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(SubRegion.Biome.LatLong), ]
  
  ### Subset the datasets
  SubRegion.Biome.CommPhylo <- match.phylo.comm(Mammal.FaurbySven.tree.1, SubRegion.Biome.Comm.)
  SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                  Phylo = SubRegion.Biome.CommPhylo$phy)
}

MPD.LatLong.Env.Biome = NULL

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate mpd values for each community
  MPD.Chiroptera.Biome <- mpd.query(SubRegion.Biome.Comm.Phylo[[j]]$Phylo, SubRegion.Biome.Comm.Phylo[[j]]$Comm)
  
  # Calculate standardized versions under the uniform model
  SES.MPD.Chiroptera.Biome <- mpd.query(SubRegion.Biome.Comm.Phylo[[j]]$Phylo, SubRegion.Biome.Comm.Phylo[[j]]$Comm, 
                                        standardize = TRUE)
  
  # Create random abundance weights
  weights = runif(length(SubRegion.Biome.Comm.Phylo[[j]]$Phylo$tip.label))
  names(weights) = SubRegion.Biome.Comm.Phylo[[j]]$Phylo$tip.label
  
  #Use query function to calculate standardized versions under the sequential model
  SES.MPD.Seq.Chiroptera.Biome <- mpd.query(SubRegion.Biome.Comm.Phylo[[j]]$Phylo, SubRegion.Biome.Comm.Phylo[[j]]$Comm, 
                                            standardize = TRUE, 
                                            null.model= "sequential",
                                            abundance.weights = weights, 
                                            reps = 1000)
  
  MPD.LatLong.Env.SubRegion.Biome <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ],
                                                      MPD = MPD.Chiroptera.Biome,
                                                      SES.MPD.Seq = SES.MPD.Seq.Chiroptera.Biome,
                                                      SES.MPD = SES.MPD.Chiroptera.Biome,
                                                      SamplingPool = rep("Biome sampling", length(MPD.Chiroptera.Biome))))
  
  MPD.LatLong.Env.Biome <- rbind(MPD.LatLong.Env.Biome, MPD.LatLong.Env.SubRegion.Biome)
}


#################################
#### Combine all data frames ####
#################################

MPD.LatLong.Env.AllScales <- rbind(MPD.LatLong.Env.World,
                                   MPD.LatLong.Env.NW.OW,
                                   MPD.LatLong.Env.Realm,
                                   MPD.LatLong.Env.Biome)


SES.MPD.Global <- subset(MPD.LatLong.Env.AllScales, SamplingPool == "Global sampling", select = c("SES.MPD"))
SES.MPD.Realm <- subset(MPD.LatLong.Env.AllScales, SamplingPool == "Realm sampling", select = c("SES.MPD"))

mean(SES.MPD.Global$SES.MPD)
mean(SES.MPD.Realm$SES.MPD)

#########################################################

gg <- ggplot(MPD.LatLong.Env.AllScales, aes(x=SamplingPool, y=SES.MPD)) + 
  geom_boxplot(aes(fill=SamplingPool)) + 
  facet_wrap(~ID_Realm) + 
  labs(x="") +
  theme_bw() + 
  theme(strip.background=element_rect(fill="black")) + 
  theme(strip.text=element_text(color="white", face="bold"))
gg




##############################################

Bats.SES.MPD.Region <- lm(MPD.Chiroptera.World ~ ID_Realm, data = MPD.regions.LatLong.Env)
summary(Bats.SES.MPD.Region)

#######################################################################################################

Bats.MPD.AnMnTemp_Bio1 <- lm(MPD.Chiroptera.World ~ AnMnTemp_Bio1, data = MPD.regions.LatLong.Env)
summary(Bats.MPD.AnMnTemp_Bio1)
plot(Bats.MPD.AnMnTemp_Bio1)

Bats.MPD.AnMnTemp_Bio1.LMM <- lmer(MPD.Chiroptera.World ~ AnMnTemp_Bio1 + (1|ID_Realm) + (1|ID_Biome), 
                                   data = MPD.regions.LatLong.Env, REML=TRUE)
summary(Bats.MPD.AnMnTemp_Bio1.LMM)


Bats.MPD.Latitude <- lm(MPD.Chiroptera.World ~ Latitude, data = MPD.regions.LatLong.Env)
summary(Bats.MPD.Latitude)


theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(NRI.LatLong.Temp, aes(MPD.AT, Env.AnMnTemp_Bio1)) + 
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  
  labs(subtitle="mpg: city vs highway mileage", 
       y="MPD.AT", 
       x="Env.AnMnTemp_Bio1", 
       title="Jittered Points")

#  scale_y_continuous(trans='log10') +
g
