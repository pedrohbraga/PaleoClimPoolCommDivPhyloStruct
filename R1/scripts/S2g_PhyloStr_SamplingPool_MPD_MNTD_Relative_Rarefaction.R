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
# The framework used here uses the sf.ses.mpd() function, which modifies the         #
# picante::ses.mpd() function to allow for parallel computation using SNOW           #
# (snowfall).                                                                        #
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2022-01-15"                                                          #
#                                                                                    # 
######################################################################################

set.seed(15145562)

# Deciding on species pool sizes

# Determine the how many species will compose the new sampling pool   

divide.sampling.pool.size.by <- 3
min.size <- divide.sampling.pool.size.by*2


# dir.create("data/forFuriosa")
# 
# saveRDS(Chiroptera.FaurSven.Phylo.ultra,
#         "data/forFuriosa/Chiroptera.FaurSven.Phylo.ultra.RDS")
# 
# saveRDS(Chiroptera.FaurSven.comm,
#         "data/forFuriosa/Chiroptera.FaurSven.comm.RDS")

saveRDS(MPD.LatLong.World,
        "data/forFuriosa/MPD.LatLong.World.RDS")

saveRDS(Chiroptera.NW.Tree,
        "data/forFuriosa/Chiroptera.NW.Tree.RDS")

saveRDS(Chiroptera.NW.Tree,
        "data/forFuriosa/Chiroptera.NW.Tree.RDS")



#########################
#### Global sampling ####
#########################

# Calculate mpd values for each community
SES.MPD.Chiroptera.World.rarefaction.relative <- opt.rarefaction.mpd(phylo.tree = Chiroptera.FaurSven.Phylo.ultra, 
                                                                     species.data = Chiroptera.FaurSven.comm, 
                                                                     n.species.rarefac = round(ncol(Chiroptera.FaurSven.comm)/divide.sampling.pool.size.by),
                                                                     n.rep = 999, 
                                                                     n.cores = 12)

# Add everything to a data.frame
SES.MPD.Chiroptera.World.rarefaction.relative$ses.mpd.z.query.rarefac$ID <- as.numeric(row.names(SES.MPD.Chiroptera.World.rarefaction.relative$ses.mpd.z.query.rarefac))

MPD.LatLong.World.rarefaction.relative <- right_join(SES.MPD.Chiroptera.World.rarefaction.relative$ses.mpd.z.query.rarefac,
                                                     MPD.LatLong.World,
                                                     by = c("ID" = "ID")
)

##########################################
#### New World and Old World sampling ####
##########################################

#### New World ####
SES.MPD.Chiroptera.NW.rarefaction.relative <- opt.rarefaction.mpd(phylo.tree = Chiroptera.NW.Tree, 
                                                                  species.data = Chiroptera.NW.Comm, 
                                                                  n.species.rarefac = round(ncol(Chiroptera.NW.Comm)/divide.sampling.pool.size.by),
                                                                  n.rep = 999, n.cores = 4)

#### Old World ####
SES.MPD.Chiroptera.OW.rarefaction.relative <-  opt.rarefaction.mpd(phylo.tree = Chiroptera.OW.Tree, 
                                                                   species.data = Chiroptera.OW.Comm, 
                                                                   n.species.rarefac = round(ncol(Chiroptera.OW.Comm)/divide.sampling.pool.size.by),
                                                                   n.rep = 999, n.cores = 4)


MPD.LatLong.NW.OW.rarefaction.relative <- rbind(SES.MPD.Chiroptera.NW.rarefaction.relative$ses.mpd.z.query.rarefac, 
                                                SES.MPD.Chiroptera.OW.rarefaction.relative$ses.mpd.z.query.rarefac) %>%
  rownames_to_column("ID") %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Hemispheric sampling") %>%
  full_join(MPD.LatLong.NW.OW,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Realm sampling pool ###  
###########################

MPD.LatLong.Realm.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mpd values for each community
  if(ncol(SubRegion.Realm.Comm.Phylo[[j]]$Comm) >= min.size){
    SES.MPD.Chiroptera.Realm.rarefaction.relative <- opt.rarefaction.mpd(phylo.tree = SubRegion.Realm.Comm.Phylo[[j]]$Phylo, 
                                                                         species.data = SubRegion.Realm.Comm.Phylo[[j]]$Comm, 
                                                                         n.species.rarefac = round(ncol(SubRegion.Realm.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                         n.rep = 999, n.cores = 4)
    
    MPD.LatLong.Realm.rarefaction.relative[[paste0(j)]] <- SES.MPD.Chiroptera.Realm.rarefaction.relative
    
    
  }
}

MPD.LatLong.Realm.rarefaction.relative <- do.call("rbind", 
                                                  lapply(MPD.LatLong.Realm.rarefaction.relative, 
                                                         "[[", 
                                                         2)) %>%
  rownames_to_column("Realm.ID") %>%
  separate(Realm.ID, into = c("ID_Realm", "ID"), sep = "\\.") %>%
  select(-ID_Realm) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Hemispheric sampling") %>%
  full_join(MPD.LatLong.Realm,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)


###########################
### Plate sampling pool ###  
###########################

MPD.LatLong.Plate.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo))){
  # Calculate mpd values for each community
  
  if(ncol(SubRegion.Plate.Comm.Phylo[[j]]$Comm) >= min.size){
    
    SES.MPD.Chiroptera.Plate.rarefaction.relative <- opt.rarefaction.mpd(phylo.tree = SubRegion.Plate.Comm.Phylo[[j]]$Phylo, 
                                                                         species.data = SubRegion.Plate.Comm.Phylo[[j]]$Comm, 
                                                                         n.species.rarefac = round(ncol(SubRegion.Plate.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                         n.rep = 999, n.cores = 4)
    
    MPD.LatLong.Plate.rarefaction.relative[[paste0(j)]] <- SES.MPD.Chiroptera.Plate.rarefaction.relative
    
  }
}


MPD.LatLong.Plate.rarefaction.relative <- do.call("rbind", 
                                                  lapply(MPD.LatLong.Plate.rarefaction.relative, 
                                                         "[[", 
                                                         2)) %>%
  rownames_to_column("PlateName.ID") %>%
  separate(PlateName.ID, into = c("ID_PlateName", "ID"), sep = "\\.") %>%
  select(-ID_PlateName) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Hemispheric sampling") %>%
  full_join(MPD.LatLong.Plate,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Biome sampling pool ###  
###########################

MPD.LatLong.Biome.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo)[!ls(SubRegion.Biome.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate mpd values for each community
  
  if(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm) >= min.size){
    
    SES.MPD.Chiroptera.Biome.rarefaction.relative <- opt.rarefaction.mpd(phylo.tree = SubRegion.Biome.Comm.Phylo[[j]]$Phylo, 
                                                                         species.data = SubRegion.Biome.Comm.Phylo[[j]]$Comm, 
                                                                         n.species.rarefac = round(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                         n.rep = 999, n.cores = 4)
    
    MPD.LatLong.Biome.rarefaction.relative[[paste0(j)]] <- SES.MPD.Chiroptera.Biome.rarefaction.relative
    
  }
}


MPD.LatLong.Biome.rarefaction.relative <- do.call("rbind", 
                                                  lapply(MPD.LatLong.Biome.rarefaction.relative, 
                                                         "[[", 
                                                         2)) %>%
  rownames_to_column("Biome.ID") %>%
  separate(Biome.ID, into = c("ID_Biome", "ID"), sep = "\\.") %>%
  select(-ID_Biome) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Biome sampling") %>%
  full_join(MPD.LatLong.Biome,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Ecoregion sampling pool ###  
###########################

MPD.LatLong.Ecoregion.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))){
  # Calculate mpd values for each community
  
  if(ncol(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm) >= min.size){
    
    SES.MPD.Chiroptera.Ecoregion.rarefaction.relative <- opt.rarefaction.mpd(phylo.tree = SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo, 
                                                                             species.data = SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm, 
                                                                             n.species.rarefac = round(ncol(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                             n.rep = 999, n.cores = 4)
    
    MPD.LatLong.Ecoregion.rarefaction.relative[[paste0(j)]] <- SES.MPD.Chiroptera.Ecoregion.rarefaction.relative
    
  }
}


MPD.LatLong.Ecoregion.rarefaction.relative <- do.call("rbind", 
                                                      lapply(MPD.LatLong.Ecoregion.rarefaction.relative, 
                                                             "[[", 
                                                             2)) %>%
  rownames_to_column("ID_Ecoregion_ID") %>%
  mutate("ID_Ecoregion_ID" = stringi::stri_sub(ID_Ecoregion_ID, -8, -1)) %>%
  separate(ID_Ecoregion_ID, into = c("ID_Ecoregion", "ID"), sep = "\\.") %>%
  select(-ID_Ecoregion) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Ecoregion sampling") %>%
  full_join(MPD.LatLong.Ecoregion,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

################################
#### Combine all data frames ###
################################

# MPD.LatLong.World.rarefaction %>% dim()
# MPD.LatLong.NW.OW.rarefaction %>% dim()
# MPD.LatLong.Realm.rarefaction %>% dim()
# MPD.LatLong.Plate.rarefaction %>% dim()
# MPD.LatLong.Biome.rarefaction %>% dim()
# MPD.LatLong.Ecoregion.rarefaction %>% dim()

MPD.LatLong.AllScales.rarefaction.relative <- rbind(MPD.LatLong.World.rarefaction.relative,
                                                    MPD.LatLong.NW.OW.rarefaction.relative,
                                                    MPD.LatLong.Realm.rarefaction.relative,
                                                    MPD.LatLong.Plate.rarefaction.relative,
                                                    MPD.LatLong.Biome.rarefaction.relative,
                                                    MPD.LatLong.Ecoregion.rarefaction.relative)


##########################################
#### Global sampling ####
##########################################

# Calculate mntd values for each community
SES.MNTD.Chiroptera.World.rarefaction.relative <- opt.rarefaction.mntd(phylo.tree = Chiroptera.FaurSven.Phylo.ultra, 
                                                                       species.data = Chiroptera.FaurSven.comm, 
                                                                       n.species.rarefac = round(ncol(Chiroptera.FaurSven.comm)/divide.sampling.pool.size.by),
                                                                       n.rep = 999, 
                                                                       n.cores = 20)

# Add everything to a data.frame
SES.MNTD.Chiroptera.World.rarefaction.relative$ses.mntd.z.query.rarefac$ID <- as.numeric(row.names(SES.MNTD.Chiroptera.World.rarefaction.relative$ses.mntd.z.query.rarefac))

MNTD.LatLong.World.rarefaction.relative <- right_join(SES.MNTD.Chiroptera.World.rarefaction.relative$ses.mntd.z.query.rarefac,
                                                      MNTD.LatLong.World,
                                                      by = c("ID" = "ID")
)


##########################################
#### New World and Old World sampling ####
##########################################

#### New World ####
SES.MNTD.Chiroptera.NW.rarefaction.relative <- opt.rarefaction.mntd(phylo.tree = Chiroptera.NW.Tree, 
                                                                    species.data = Chiroptera.NW.Comm, 
                                                                    n.species.rarefac = round(ncol(Chiroptera.NW.Comm)/divide.sampling.pool.size.by),
                                                                    n.rep = 999, n.cores = 8)

#### Old World ####
SES.MNTD.Chiroptera.OW.rarefaction.relative <-  opt.rarefaction.mntd(phylo.tree = Chiroptera.OW.Tree, 
                                                                     species.data = Chiroptera.OW.Comm, 
                                                                     n.species.rarefac = round(ncol(Chiroptera.OW.Comm)/divide.sampling.pool.size.by),
                                                                     n.rep = 999, n.cores = 8)


MNTD.LatLong.NW.OW.rarefaction.relative <- rbind(SES.MNTD.Chiroptera.NW.rarefaction.relative$ses.mntd.z.query.rarefac, 
                                                 SES.MNTD.Chiroptera.OW.rarefaction.relative$ses.mntd.z.query.rarefac) %>%
  rownames_to_column("ID") %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Hemispheric sampling") %>%
  full_join(MNTD.LatLong.NW.OW,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Realm sampling pool ###  
###########################

MNTD.LatLong.Realm.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mntd values for each community
  
  SES.MNTD.Chiroptera.Realm.rarefaction <- opt.rarefaction.mntd(phylo.tree = SubRegion.Realm.Comm.Phylo[[j]]$Phylo, 
                                                                species.data = SubRegion.Realm.Comm.Phylo[[j]]$Comm, 
                                                                n.species.rarefac = round(ncol(SubRegion.Realm.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                n.rep = 999, n.cores = 4)
  
  MNTD.LatLong.Realm.rarefaction.relative[[paste0(j)]] <- SES.MNTD.Chiroptera.Realm.rarefaction
  
  
}

MNTD.LatLong.Realm.rarefaction.relative <- do.call("rbind", 
                                                   lapply(MNTD.LatLong.Realm.rarefaction.relative, 
                                                          "[[", 
                                                          2)) %>%
  rownames_to_column("Realm.ID") %>%
  separate(Realm.ID, into = c("ID_Realm", "ID"), sep = "\\.") %>%
  select(-ID_Realm) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Realm sampling") %>%
  full_join(MNTD.LatLong.Realm,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)


###########################
### Plate sampling pool ###  
###########################

MNTD.LatLong.Plate.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo))){
  # Calculate mntd values for each community
  
  if(ncol(SubRegion.Plate.Comm.Phylo[[j]]$Comm) >= min.size){
    
    SES.MNTD.Chiroptera.Plate.rarefaction.relative <- opt.rarefaction.mntd(phylo.tree = SubRegion.Plate.Comm.Phylo[[j]]$Phylo, 
                                                                           species.data = SubRegion.Plate.Comm.Phylo[[j]]$Comm, 
                                                                           n.species.rarefac = round(ncol(SubRegion.Plate.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                           n.rep = 999, n.cores = 4)
    
    MNTD.LatLong.Plate.rarefaction.relative[[paste0(j)]] <- SES.MNTD.Chiroptera.Plate.rarefaction.relative
    
  }
}


MNTD.LatLong.Plate.rarefaction.relative <- do.call("rbind", 
                                                   lapply(MNTD.LatLong.Plate.rarefaction.relative, 
                                                          "[[", 
                                                          2)) %>%
  rownames_to_column("PlateName.ID") %>%
  separate(PlateName.ID, into = c("ID_PlateName", "ID"), sep = "\\.") %>%
  select(-ID_PlateName) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Plate sampling") %>%
  full_join(MNTD.LatLong.Plate,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Biome sampling pool ###  
###########################

MNTD.LatLong.Biome.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo)[!ls(SubRegion.Biome.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate mntd values for each community
  
  if(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm) >= min.size){
    
    SES.MNTD.Chiroptera.Biome.rarefaction.relative <- opt.rarefaction.mntd(phylo.tree = SubRegion.Biome.Comm.Phylo[[j]]$Phylo, 
                                                                           species.data = SubRegion.Biome.Comm.Phylo[[j]]$Comm, 
                                                                           n.species.rarefac = round(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                           n.rep = 999, n.cores = 4)
    
    MNTD.LatLong.Biome.rarefaction.relative[[paste0(j)]] <- SES.MNTD.Chiroptera.Biome.rarefaction.relative
    
  }
}


MNTD.LatLong.Biome.rarefaction.relative <- do.call("rbind", 
                                                   lapply(MNTD.LatLong.Biome.rarefaction.relative, 
                                                          "[[", 
                                                          2)) %>%
  rownames_to_column("Biome.ID") %>%
  separate(Biome.ID, into = c("ID_Biome", "ID"), sep = "\\.") %>%
  select(-ID_Biome) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Biome sampling") %>%
  full_join(MNTD.LatLong.Biome,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Ecoregion sampling pool ###  
###########################

MNTD.LatLong.Ecoregion.rarefaction.relative <- list()

# for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))){
  # Calculate mntd values for each community
  
  if(ncol(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm) >= min.size){
    
    SES.MNTD.Chiroptera.Ecoregion.rarefaction.relative <- opt.rarefaction.mntd(phylo.tree = SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo, 
                                                                               species.data = SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm, 
                                                                               n.species.rarefac = round(ncol(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                               n.rep = 999, n.cores = 4)
    
    MNTD.LatLong.Ecoregion.rarefaction.relative[[paste0(j)]] <- SES.MNTD.Chiroptera.Ecoregion.rarefaction.relative
    
  }
}


MNTD.LatLong.Ecoregion.rarefaction.relative <- do.call("rbind", 
                                                       lapply(MNTD.LatLong.Ecoregion.rarefaction.relative, 
                                                              "[[", 
                                                              2)) %>%
  rownames_to_column("ID_Ecoregion_ID") %>%
  mutate("ID_Ecoregion_ID" = stringi::stri_sub(ID_Ecoregion_ID, -8, -1)) %>%
  separate(ID_Ecoregion_ID, into = c("ID_Ecoregion", "ID"), sep = "\\.") %>%
  select(-ID_Ecoregion) %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Ecoregion sampling") %>%
  full_join(MNTD.LatLong.Ecoregion,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

################################
#### Combine all data frames ###
################################

# MNTD.LatLong.World.rarefaction %>% dim()
# MNTD.LatLong.NW.OW.rarefaction %>% dim()
# MNTD.LatLong.Realm.rarefaction %>% dim()
# MNTD.LatLong.Plate.rarefaction %>% dim()
# MNTD.LatLong.Biome.rarefaction %>% dim()
# MNTD.LatLong.Ecoregion.rarefaction.relative %>% dim()

MNTD.LatLong.AllScales.rarefaction.relative <- rbind(MNTD.LatLong.World.rarefaction.relative,
                                                     MNTD.LatLong.NW.OW.rarefaction.relative,
                                                     MNTD.LatLong.Realm.rarefaction.relative,
                                                     MNTD.LatLong.Plate.rarefaction.relative,
                                                     MNTD.LatLong.Biome.rarefaction.relative,
                                                     MNTD.LatLong.Ecoregion.rarefaction.relative)

MNTD.LatLong.AllScales.rarefaction.relative %>% filter(SamplingPool == "Realm sampling")

######################################################################################
#### Code to merge, represent the MPD and MNTD, and test for the effects of       ####
#### sampling pool restriction in the phylogenetic structuring of communities     ####
#                                                                                    #
# Description: This code merges the MPD.LatLong.AllScales.rarefaction and
# MNTD.LatLong.AllScales.rarefaction data sets into a single one. It then
# represents NRI and NTI across the gradient of sampling pool restrictions. It
# finally tests the effects of sampling pool (or species pool) restriction on
# the phylogenetic community structure of bats.
# #
#                                                                                    #
# We performed separate robust heteroscedastic repeated measurement                  #
# analyses of variance based on 20% trimmed-means (Wilcox, 1997, 2012). For each     #
# realm, we estimated whether NRI and NTI (response variables, in separate           #       
# tests) of bat communities (dependent-groups) differed across the discrete          #
# sampling pool restriction gradient (group-levels). Finally, we performed           #
# post-hoc comparisons using Hochberg’s approach to control for family-wise          #
# error (Hochberg, 1988; Wilcox, 2012).                                              #
#                                                                                    #
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2021-09-16"                                                          #
#                                                                                    # 
######################################################################################

#########################################
### Merge data into a single data set ###
#########################################

MPD.MNTD.LatLong.AllScales.rarefaction.relative <- MPD.LatLong.AllScales.rarefaction.relative %>%
  unite(ID_SamplingPool, 
        ID, 
        SamplingPool, 
        remove = FALSE) %>%
  filter(!is.na(SamplingPool)) %>%
  left_join(MNTD.LatLong.AllScales.rarefaction.relative %>% 
              unite(ID_SamplingPool, 
                    ID, 
                    SamplingPool, 
                    remove = FALSE) %>%
              filter(!is.na(SamplingPool)) %>%
              select(-c("ntaxa", "ID", "perms", 
                        "ntaxa.pool.size", "rarefac.pool.size",
                        LONB:SamplingPool)),
            by = "ID_SamplingPool")

# Refactor data

MPD.MNTD.LatLong.AllScales.rarefaction.relative$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales.rarefaction.relative$ID_Realm, 
                                                                   levels=c('Neotropical',
                                                                            'Nearctic',
                                                                            'Afrotropical',
                                                                            'Palearctic', 
                                                                            'Indomalay', 
                                                                            'Australasian'))

MPD.MNTD.LatLong.AllScales.rarefaction.relative$SamplingPool <-  factor(MPD.MNTD.LatLong.AllScales.rarefaction.relative$SamplingPool,  
                                                                        levels = c("Global sampling",
                                                                                   "Hemispheric sampling",
                                                                                   "Realm sampling",
                                                                                   "Plate sampling",
                                                                                   "Biome sampling",
                                                                                   "Ecoregion sampling"))

MPD.MNTD.LatLong.AllScales.rarefaction.relative$ID_Biome_Acronym = factor(MPD.MNTD.LatLong.AllScales.rarefaction.relative$ID_Biome, 
                                                                          labels = abbreviate(
                                                                            gsub("_", 
                                                                                 " ",
                                                                                 levels(MPD.LatLong.AllScales.rarefaction.relative$ID_Biome))))

write.csv(MPD.MNTD.LatLong.AllScales.rarefaction.relative, 
          "data/matrices/MPD.MNTD.LatLong.AllScales.rarefaction.relative.sqrt.csv")

 
# MPD.MNTD.LatLong.AllScales.rarefaction.relative <- read.csv("data/matrices/MPD.MNTD.LatLong.AllScales.rarefaction.relative.csv", h = T, row.names = 1)

# MPD.MNTD.LatLong.AllScales.rarefaction.relative.sqrt <- read.csv("data/matrices/MPD.MNTD.LatLong.AllScales.rarefaction.relative.sqrt.csv", h = T)

############################################################################ 
### Representing the phylogenetic community structure of bats across the ###
### gradient of sampling pool (or species pool) restriction              ###
############################################################################

# Representing NRI wrapping across Realms

(NRI.Realm.rarefaction.relative.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.rarefaction.relative,
                                                         is.na(ID_Realm) == FALSE),
                                                  aes(x = SamplingPool, 
                                                      y = nri.rarefac.mean)) +
   geom_boxplot(aes(fill = SamplingPool)) +
   facet_wrap(~ID_Realm,
              nrow = 1,
              strip.position = "bottom") +
   ggsci::scale_fill_uchicago() +
   # scale_fill_nord("baie_mouton", reverse = TRUE) +
   # scale_fill_okabeito(reverse = TRUE) +
   # scale_fill_brewer(palette = "Dark2") +
   # scale_fill_viridis(discrete = TRUE,
   #                    name="Sampling Pool") +
   scale_y_continuous(breaks = round(
     #pretty(MPD.LatLong.AllScales.rarefaction$NRI, n = 7)
     c(min(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nri.rarefac.mean, na.rm = T),
       0,
       5,
       10,
       15,
       max(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nri.rarefac.mean, na.rm = T))),
     limits = c(min(MPD.LatLong.AllScales.rarefaction.relative$nri.rarefac.mean, na.rm = T) - 0.1,
                max(MPD.LatLong.AllScales.rarefaction.relative$nri.rarefac.mean, na.rm = T) + 0.1)
   ) +
   geom_hline(yintercept = 0,
              alpha = 0.4) +
   labs(y="NRI", x = "") +
   theme_minimal() +
   theme_minimal(base_size = 20) +
   theme(strip.background = element_rect(fill = "white",
                                         linetype = NULL,
                                         color = "white"),
         strip.text = element_text(color = "black",
                                   face = "bold",
                                   size = 15),
         axis.text.y = element_text(face = "bold"),
         axis.text.x = element_blank(),
         axis.title = element_text(size = 16, face="bold"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = "bottom",
         legend.text = element_text(size = 14),
         legend.title = element_text(size = 16),
         legend.title.align = 0.5,
         legend.margin=margin(0, 0, 0, 0),
         legend.box.margin=margin(-5, -5, 2, -5)
   ) +   
   guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                              title.position = "bottom",
                              title = "Lower  ⬅  Geographical Extent Restriction  ➡  Higher",)
   )
)


(NTI.Realm.rarefaction.relative.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.rarefaction.relative,
                                                         is.na(ID_Realm) == FALSE),
                                                  aes(x = SamplingPool, y = nti.rarefac.mean)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    facet_wrap(~ID_Realm,
               nrow = 1,
               strip.position = "bottom") +
    ggsci::scale_fill_uchicago() +
    scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nti.rarefac.mean, n = 5),
                       # round(c(min(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nti, na.rm = T),
                       #   0,
                       #   5,
                       #   10,
                       #   15,
                       #   max(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nti, na.rm = T)))),
                       limits = c(min(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nti.rarefac.mean, na.rm = T) - 0.1,
                                  max(MPD.MNTD.LatLong.AllScales.rarefaction.relative$nti.rarefac.mean, na.rm = T) + 0.1)
    ) +
    geom_hline(yintercept = 0,
               alpha = 0.4) +
    labs(y = "NTI", x = "") +
    theme_minimal(base_size = 20) +
    theme(strip.background = element_rect(fill = "white",
                                          linetype = NULL,
                                          color = "white"),
          strip.text = element_text(color = "black",
                                    face = "bold",
                                    size = 15),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15),
          legend.title.align = 0.5,
          legend.margin=margin(0, 0, 0, 0),
          legend.box.margin=margin(-5, -5, -5, -5)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                               title.position = "bottom",
                               title = "Lower  ⬅  Geographical Extent Restriction  ➡  Higher",)
    )
)

### Combining plots into a single figure ###

(fig.NRI.NTI.Realm.rarefaction.relative.boxplot <- ggarrange(
  NRI.Realm.rarefaction.relative.boxplot +
    theme(axis.title.x = element_blank()), 
  NTI.Realm.rarefaction.relative.boxplot,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  widths = c(1, 1),
  heights = c(1, 1.1),
  align = "v",
  common.legend = TRUE,
  legend = "bottom")
)

### Exporting the box-plot to PNG ###

# This is the Figure 2 in the manuscript

ggsave(filename = "figures/fig.NRI.NTI.Realm.rarefaction.relative.boxplot_1_3.sampling.pool.png", 
       dpi = 300, 
       width = 18.5, height = 10, 
       units = "in")


####

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  select(-c("nri.rarefac.mean", "nti.rarefac.mean",  "ses.mpd.z.query.rarefac.mean",  "ses.mpd.z.query.rarefac.mean" )) %>%
  left_join(MPD.MNTD.LatLong.AllScales.rarefaction.relative %>%
              select(c("ID_SamplingPool", "nri.rarefac.mean",
                       "nti.rarefac.mean",  "ses.mpd.z.query.rarefac.mean",  "ses.mpd.z.query.rarefac.mean" )),
            by = "ID_SamplingPool"
            )

saveRDS(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div, 
        "data/matrices/MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div.RDS")
