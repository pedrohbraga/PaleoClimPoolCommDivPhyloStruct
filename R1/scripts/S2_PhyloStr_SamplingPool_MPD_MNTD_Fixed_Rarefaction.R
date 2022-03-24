# Deciding on species pool sizes

# Hemispheric sampling

MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Hemispheric sampling") %>%
  select("ntaxa.pool.size") %>%
  distinct()

MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Realm sampling") %>%
  select("ID_Realm" , "ntaxa.pool.size") %>%
  distinct() %>%
  arrange(desc(ntaxa.pool.size))

MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Plate sampling") %>%
  select("ID_PlateName", "ntaxa.pool.size") %>%
  distinct() %>%
  arrange(desc(ntaxa.pool.size))

MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Biome sampling") %>%
  select("ID_Biome_Realm", "ntaxa.pool.size") %>%
  distinct() %>%
  arrange(desc(ntaxa.pool.size))

MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Ecoregion sampling") %>%
  select("ID_Ecoregion", "ntaxa.pool.size") %>%
  distinct() %>%
  arrange(desc(ntaxa.pool.size))


MPD.LatLong.AllScales %>%
  group_by(SamplingPool) %>%
  filter(SamplingPool == "Hemispheric sampling") %>%
  select("ntaxa.pool.size") %>%
  distinct() %>%
  group_by("ntaxa.pool.size")


# Determine the how many species will compose the new sampling pool   

divide.sampling.pool.size.by <- 2

n.species.rarefac.no = 20

#########################
#### Global sampling ####
#########################


# Calculate mpd values for each community
SES.MPD.Chiroptera.World.rarefaction.fixed <- opt.rarefaction.mpd(phylo.tree = Chiroptera.FaurSven.Phylo.ultra, 
                                                                  species.data = Chiroptera.FaurSven.comm, 
                                                                  n.species.rarefac = n.species.rarefac.no,
                                                                  n.rep = 99, n.cores = 8)

# Add everything to a data.frame
SES.MPD.Chiroptera.World.rarefaction.fixed$ses.mpd.z.query.rarefac$ID <- as.numeric(row.names(SES.MPD.Chiroptera.World.rarefaction.fixed$ses.mpd.z.query.rarefac))

MPD.LatLong.World.rarefaction.fixed <- right_join(SES.MPD.Chiroptera.World.rarefaction.fixed$ses.mpd.z.query.rarefac,
                                                  MPD.LatLong.World,
                                                  by = c("ID" = "ID")
)


##########################################
#### New World and Old World sampling ####
##########################################

#### New World ####
SES.MPD.Chiroptera.NW.rarefaction.fixed <- opt.rarefaction.mpd(phylo.tree = Chiroptera.NW.Tree, 
                                                               species.data = Chiroptera.NW.Comm, 
                                                               n.species.rarefac = n.species.rarefac.no,
                                                               n.rep = 99, n.cores = 8)

#### Old World ####
SES.MPD.Chiroptera.OW.rarefaction.fixed <-  opt.rarefaction.mpd(phylo.tree = Chiroptera.OW.Tree, 
                                                                species.data = Chiroptera.OW.Comm, 
                                                                n.species.rarefac = n.species.rarefac.no,
                                                                n.rep = 99, n.cores = 4)


MPD.LatLong.NW.OW.rarefaction.fixed <- rbind(SES.MPD.Chiroptera.NW.rarefaction.fixed$ses.mpd.z.query.rarefac, 
                                             SES.MPD.Chiroptera.OW.rarefaction.fixed$ses.mpd.z.query.rarefac) %>%
  rownames_to_column("ID") %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Hemispheric sampling") %>%
  full_join(MPD.LatLong.NW.OW,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Realm sampling pool ###  
###########################

MPD.LatLong.Realm.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mpd values for each community
  
  SES.MPD.Chiroptera.Realm.rarefaction.fixed <- opt.rarefaction.mpd(phylo.tree = SubRegion.Realm.Comm.Phylo[[j]]$Phylo, 
                                                                    species.data = SubRegion.Realm.Comm.Phylo[[j]]$Comm, 
                                                                    n.species.rarefac = n.species.rarefac.no,
                                                                    n.rep = 99, n.cores = 8)
  
  MPD.LatLong.Realm.rarefaction.fixed[[paste0(j)]] <- SES.MPD.Chiroptera.Realm.rarefaction.fixed
  
  
}

MPD.LatLong.Realm.rarefaction.fixed <- do.call("rbind", 
                                               lapply(MPD.LatLong.Realm.rarefaction.fixed, 
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

MPD.LatLong.Plate.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo))){
  # Calculate mpd values for each community
  
  if(ncol(SubRegion.Plate.Comm.Phylo[[j]]$Comm) >= n.species.rarefac.no){
    
    SES.MPD.Chiroptera.Plate.rarefaction.fixed <- opt.rarefaction.mpd(phylo.tree = SubRegion.Plate.Comm.Phylo[[j]]$Phylo, 
                                                                      species.data = SubRegion.Plate.Comm.Phylo[[j]]$Comm, 
                                                                      n.species.rarefac = n.species.rarefac.no,
                                                                      n.rep = 99, n.cores = 8)
    
    MPD.LatLong.Plate.rarefaction.fixed[[paste0(j)]] <- SES.MPD.Chiroptera.Plate.rarefaction.fixed
    
  }
}


MPD.LatLong.Plate.rarefaction.fixed <- do.call("rbind", 
                                               lapply(MPD.LatLong.Plate.rarefaction.fixed, 
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

MPD.LatLong.Biome.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo)[!ls(SubRegion.Biome.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate mpd values for each community
  
  if(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm) >= n.species.rarefac.no){
    
    SES.MPD.Chiroptera.Biome.rarefaction.fixed <- opt.rarefaction.mpd(phylo.tree = SubRegion.Biome.Comm.Phylo[[j]]$Phylo, 
                                                                      species.data = SubRegion.Biome.Comm.Phylo[[j]]$Comm, 
                                                                      n.species.rarefac = round(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm)/divide.sampling.pool.size.by),
                                                                      n.rep = 99, n.cores = 8)
    
    MPD.LatLong.Biome.rarefaction.fixed[[paste0(j)]] <- SES.MPD.Chiroptera.Biome.rarefaction.fixed
    
  }
}


MPD.LatLong.Biome.rarefaction.fixed <- do.call("rbind", 
                                               lapply(MPD.LatLong.Biome.rarefaction.fixed, 
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

MPD.LatLong.Ecoregion.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))){
  # Calculate mpd values for each community
  
  if(ncol(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm) >= n.species.rarefac.no){
    
    SES.MPD.Chiroptera.Ecoregion.rarefaction.fixed <- opt.rarefaction.mpd(phylo.tree = SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo, 
                                                                          species.data = SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm, 
                                                                          n.species.rarefac = n.species.rarefac.no,
                                                                          n.rep = 99, n.cores = 4)
    
    MPD.LatLong.Ecoregion.rarefaction.fixed[[paste0(j)]] <- SES.MPD.Chiroptera.Ecoregion.rarefaction.fixed
    
  }
}


MPD.LatLong.Ecoregion.rarefaction.fixed <- do.call("rbind", 
                                                   lapply(MPD.LatLong.Ecoregion.rarefaction.fixed, 
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

MPD.LatLong.AllScales.rarefaction.fixed <- rbind(MPD.LatLong.World.rarefaction.fixed,
                                                 MPD.LatLong.NW.OW.rarefaction.fixed,
                                                 MPD.LatLong.Realm.rarefaction.fixed,
                                                 MPD.LatLong.Plate.rarefaction.fixed,
                                                 MPD.LatLong.Biome.rarefaction.fixed,
                                                 MPD.LatLong.Ecoregion.rarefaction.fixed)


##########################################
#### Global sampling ####
##########################################

# Calculate mntd values for each community
SES.MNTD.Chiroptera.World.rarefaction.fixed <- opt.rarefaction.mntd(phylo.tree = Chiroptera.FaurSven.Phylo.ultra, 
                                                                    species.data = Chiroptera.FaurSven.comm, 
                                                                    n.species.rarefac = n.species.rarefac.no,
                                                                    n.rep = 99, n.cores = 8)

# Add everything to a data.frame
SES.MNTD.Chiroptera.World.rarefaction.fixed$ses.mntd.z.query.rarefac$ID <- as.numeric(row.names(SES.MNTD.Chiroptera.World.rarefaction.fixed$ses.mntd.z.query.rarefac))

MNTD.LatLong.World.rarefaction.fixed <- right_join(SES.MNTD.Chiroptera.World.rarefaction.fixed$ses.mntd.z.query.rarefac,
                                                   MNTD.LatLong.World,
                                                   by = c("ID" = "ID")
)


##########################################
#### New World and Old World sampling ####
##########################################

#### New World ####
SES.MNTD.Chiroptera.NW.rarefaction.fixed <- opt.rarefaction.mntd(phylo.tree = Chiroptera.NW.Tree, 
                                                                 species.data = Chiroptera.NW.Comm, 
                                                                 n.species.rarefac = n.species.rarefac.no,
                                                                 n.rep = 99, n.cores = 8)

#### Old World ####
SES.MNTD.Chiroptera.OW.rarefaction.fixed <-  opt.rarefaction.mntd(phylo.tree = Chiroptera.OW.Tree, 
                                                                  species.data = Chiroptera.OW.Comm, 
                                                                  n.species.rarefac = n.species.rarefac.no,
                                                                  n.rep = 99, n.cores = 4)


MNTD.LatLong.NW.OW.rarefaction.fixed <- rbind(SES.MNTD.Chiroptera.NW.rarefaction.fixed$ses.mntd.z.query.rarefac, 
                                              SES.MNTD.Chiroptera.OW.rarefaction.fixed$ses.mntd.z.query.rarefac) %>%
  rownames_to_column("ID") %>%
  mutate_at("ID", as.numeric) %>%
  #  mutate(SamplingPool = "Hemispheric sampling") %>%
  full_join(MNTD.LatLong.NW.OW,
            by = c("ID" = "ID")) %>%
  set_rownames(.$ID)

###########################
### Realm sampling pool ###  
###########################

MNTD.LatLong.Realm.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mntd values for each community
  
  SES.MNTD.Chiroptera.Realm.rarefaction <- opt.rarefaction.mntd(phylo.tree = SubRegion.Realm.Comm.Phylo[[j]]$Phylo, 
                                                                species.data = SubRegion.Realm.Comm.Phylo[[j]]$Comm, 
                                                                n.species.rarefac = n.species.rarefac.no,
                                                                n.rep = 99, n.cores = 8)
  
  MNTD.LatLong.Realm.rarefaction.fixed[[paste0(j)]] <- SES.MNTD.Chiroptera.Realm.rarefaction
  
  
}

MNTD.LatLong.Realm.rarefaction.fixed <- do.call("rbind", 
                                                lapply(MNTD.LatLong.Realm.rarefaction.fixed, 
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

MNTD.LatLong.Plate.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo))){
  # Calculate mntd values for each community
  
  if(ncol(SubRegion.Plate.Comm.Phylo[[j]]$Comm) >= n.species.rarefac.no){
    
    SES.MNTD.Chiroptera.Plate.rarefaction.fixed <- opt.rarefaction.mntd(phylo.tree = SubRegion.Plate.Comm.Phylo[[j]]$Phylo, 
                                                                        species.data = SubRegion.Plate.Comm.Phylo[[j]]$Comm, 
                                                                        n.species.rarefac = n.species.rarefac.no,
                                                                        n.rep = 99, n.cores = 8)
    
    MNTD.LatLong.Plate.rarefaction.fixed[[paste0(j)]] <- SES.MNTD.Chiroptera.Plate.rarefaction.fixed
    
  }
}


MNTD.LatLong.Plate.rarefaction.fixed <- do.call("rbind", 
                                                lapply(MNTD.LatLong.Plate.rarefaction.fixed, 
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

MNTD.LatLong.Biome.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo)[!ls(SubRegion.Biome.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate mntd values for each community
  
  if(ncol(SubRegion.Biome.Comm.Phylo[[j]]$Comm) >= n.species.rarefac.no){
    
    SES.MNTD.Chiroptera.Biome.rarefaction.fixed <- opt.rarefaction.mntd(phylo.tree = SubRegion.Biome.Comm.Phylo[[j]]$Phylo, 
                                                                        species.data = SubRegion.Biome.Comm.Phylo[[j]]$Comm, 
                                                                        n.species.rarefac = n.species.rarefac.no,
                                                                        n.rep = 99, n.cores = 8)
    
    MNTD.LatLong.Biome.rarefaction.fixed[[paste0(j)]] <- SES.MNTD.Chiroptera.Biome.rarefaction.fixed
    
  }
}


MNTD.LatLong.Biome.rarefaction.fixed <- do.call("rbind", 
                                                lapply(MNTD.LatLong.Biome.rarefaction.fixed, 
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

MNTD.LatLong.Ecoregion.rarefaction.fixed <- list()

# for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){

for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))){
  # Calculate mntd values for each community
  
  if(ncol(SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm) >= n.species.rarefac.no){
    
    SES.MNTD.Chiroptera.Ecoregion.rarefaction.fixed <- opt.rarefaction.mntd(phylo.tree = SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo, 
                                                                            species.data = SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm, 
                                                                            n.species.rarefac = n.species.rarefac.no,
                                                                            n.rep = 99, n.cores = 4)
    
    MNTD.LatLong.Ecoregion.rarefaction.fixed[[paste0(j)]] <- SES.MNTD.Chiroptera.Ecoregion.rarefaction.fixed
    
  }
}


MNTD.LatLong.Ecoregion.rarefaction.fixed <- do.call("rbind", 
                                                    lapply(MNTD.LatLong.Ecoregion.rarefaction.fixed, 
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
# MNTD.LatLong.Ecoregion.rarefaction.fixed %>% dim()

MNTD.LatLong.AllScales.rarefaction.fixed <- rbind(MNTD.LatLong.World.rarefaction.fixed,
                                                  MNTD.LatLong.NW.OW.rarefaction.fixed,
                                                  MNTD.LatLong.Realm.rarefaction.fixed,
                                                  MNTD.LatLong.Plate.rarefaction.fixed,
                                                  MNTD.LatLong.Biome.rarefaction.fixed,
                                                  MNTD.LatLong.Ecoregion.rarefaction.fixed)

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

MPD.MNTD.LatLong.AllScales.rarefaction.fixed <- MPD.LatLong.AllScales.rarefaction.fixed %>%
  unite(ID_SamplingPool, 
        ID, 
        SamplingPool, 
        remove = FALSE) %>%
  filter(!is.na(SamplingPool)) %>%
  left_join(MNTD.LatLong.AllScales.rarefaction.fixed %>% 
              unite(ID_SamplingPool, 
                    ID, 
                    SamplingPool, 
                    remove = FALSE) %>%
              filter(!is.na(SamplingPool)) %>%
              select(-c("ntaxa", "ID", "perms", "ntaxa.pool.size", "rarefac.pool.size",
                        LONB:SamplingPool)),
            by = "ID_SamplingPool")

# Refactor data

MPD.MNTD.LatLong.AllScales.rarefaction.fixed$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$ID_Realm, 
                                                                levels=c('Neotropical',
                                                                         'Nearctic',
                                                                         'Afrotropical',
                                                                         'Palearctic', 
                                                                         'Indomalay', 
                                                                         'Australasian'))

MPD.MNTD.LatLong.AllScales.rarefaction.fixed$SamplingPool <-  factor(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$SamplingPool,  
                                                                     levels = c("Global sampling",
                                                                                "Hemispheric sampling",
                                                                                "Realm sampling",
                                                                                "Plate sampling",
                                                                                "Biome sampling",
                                                                                "Ecoregion sampling"))

MPD.MNTD.LatLong.AllScales.rarefaction.fixed$ID_Biome_Acronym = factor(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$ID_Biome, 
                                                                       labels = abbreviate(gsub("_", 
                                                                                                " ",
                                                                                                levels(MPD.LatLong.AllScales.rarefaction.fixed$ID_Biome))))

write.csv(MPD.MNTD.LatLong.AllScales.rarefaction.fixed,
          "data/matrices/MPD.MNTD.LatLong.AllScales.rarefaction.fixed.csv")

############################################################################ 
### Representing the phylogenetic community structure of bats across the ###
### gradient of sampling pool (or species pool) restriction              ###
############################################################################

# Representing NRI wrapping across Realms

(NRI.Realm.rarefaction.fixed.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.rarefaction.fixed,
                                                      is.na(ID_Realm) == FALSE),
                                               aes(x = SamplingPool, y = nri.rarefac.mean)) +
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
     c(min(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T),
       0,
       5,
       10,
       15,
       max(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T))),
     limits = c(min(MPD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T) - 0.1,
                max(MPD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T) + 0.1)
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
                              title="Lower  ⬅  Sampling Pool Restriction  ➡  Higher",)
   )
)


(NTI.Realm.rarefaction.fixed.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.rarefaction.fixed,
                                                      is.na(ID_Realm) == FALSE),
                                               aes(x = SamplingPool, y = nti.rarefac.mean)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    facet_wrap(~ID_Realm,
               nrow = 1,
               strip.position = "bottom") +
    ggsci::scale_fill_uchicago() +
    scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nti.rarefac.mean, n = 5),
                       # round(c(min(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nti, na.rm = T),
                       #   0,
                       #   5,
                       #   10,
                       #   15,
                       #   max(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nti, na.rm = T)))),
                       limits = c(min(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nti.rarefac.mean, na.rm = T) - 0.1,
                                  max(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nti.rarefac.mean, na.rm = T) + 0.1)
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
                               title="Lower  ⬅  Sampling Pool Restriction  ➡  Higher",)
    )
)

### Combining plots into a single figure ###

(fig.NRI.NTI.Realm.rarefaction.fixed.boxplot <- ggarrange(
  NRI.Realm.rarefaction.fixed.boxplot +
    theme(axis.title.x = element_blank()), 
  NTI.Realm.rarefaction.fixed.boxplot,
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

ggsave(filename = "figures/fig.NRI.NTI.Realm.rarefaction.fixed.boxplot_90.sampling.pool.png", 
       dpi = 300, 
       width = 18.5, height = 10, 
       units = "in")


(NRI.rarefaction.fixed.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.rarefaction.fixed,
                                                      is.na(ID_Realm) == FALSE),
                                               aes(x = SamplingPool, y = nri.rarefac.mean)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    # facet_wrap(~ID_Realm,
    #            nrow = 1,
    #            strip.position = "bottom") +
    ggsci::scale_fill_uchicago() +
    # scale_fill_nord("baie_mouton", reverse = TRUE) +
    # scale_fill_okabeito(reverse = TRUE) +
    # scale_fill_brewer(palette = "Dark2") +
    # scale_fill_viridis(discrete = TRUE,
    #                    name="Sampling Pool") +
    scale_y_continuous(breaks = round(
      #pretty(MPD.LatLong.AllScales.rarefaction$NRI, n = 7)
      c(min(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T),
        0,
        5,
        10,
        15,
        max(MPD.MNTD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T))),
      limits = c(min(MPD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T) - 0.1,
                 max(MPD.LatLong.AllScales.rarefaction.fixed$nri.rarefac.mean, na.rm = T) + 0.1)
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
                               title="Lower  ⬅  Sampling Pool Restriction  ➡  Higher",)
    )
)
