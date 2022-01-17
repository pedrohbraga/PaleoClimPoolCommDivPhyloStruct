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

phylo.sample.pool <- sample(1:1000, 50)[3:50]

starting.time <- Sys.time()

# phylo.sampled.pos =  phylo.sample.pool[1]

for(phylo.sampled.pos in phylo.sample.pool){
  
  message("Estimating phylogenetic structure using the phylogenetic tree: ", phylo.sampled.pos,
          ".\nPhylogeny ",
          which(phylo.sampled.pos == phylo.sample.pool),
          " from ",
          length(phylo.sample.pool),
          ".\n")  
  
  Chiroptera.FaurSven.tree.sampled <- match.phylo.comm.(Mammal.FaurSven.trees[[phylo.sampled.pos]], 
                                                        Chiroptera.Comm)$phy %>%
    force.ultrametric()
  
  Chiroptera.FaurSven.comm.sampled <- match.phylo.comm.(Mammal.FaurSven.trees[[phylo.sampled.pos]], 
                                                        Chiroptera.Comm)$comm
  
  # Begin of the computation of SES.MPD values across the sampling pool
  # restriction gradient
  
  ########################
  ### Global sampling ####
  ########################
  
  # Calculate mpd values for each community
  SES.MPD.Chiroptera.World <- mod.ses.mpd.query.sf(species.data = Chiroptera.FaurSven.comm.sampled,
                                                   phylo.tree = Chiroptera.FaurSven.tree.sampled, 
                                                   perms = 999, cores = 25)
  # Add everything to a data.frame
  SES.MPD.Chiroptera.World$ID <- as.numeric(SES.MPD.Chiroptera.World$ID)
  
  MPD.LatLong.World <- right_join(SES.MPD.Chiroptera.World,
                                  world_grid_50km_cat_df,
                                  by = c("ID" = "ID")) %>%
    mutate(SamplingPool = rep("Global sampling", 
                              nrow(SES.MPD.Chiroptera.World))
    )
  
  ##########################################
  #### New World and Old World sampling ####
  ##########################################
  
  # Defining the New world and the Old World
  NW.LatLong <- subset(world_grid_50km_cat_df, 
                       ID_Realm == "Nearctic" | ID_Realm == "Neotropical")
  
  OW.LatLong <- subset(world_grid_50km_cat_df, 
                       (ID_Realm %in% c("Palearctic", "Australasian",
                                        "Afrotropical", "Indomalay",
                                        "Oceanic")
                       )
  )
  
  nrow(NW.LatLong) + nrow(OW.LatLong); nrow(world_grid_50km_cat_df) # Checking the number of rows
  
  # Extracting NW communities
  NW.Comm <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(NW.LatLong), ]
  NW.Comm <- NW.Comm[ , colSums(NW.Comm) >= 1]
  
  # Extracting OW communities
  OW.Comm <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(OW.LatLong), ]
  OW.Comm <- OW.Comm[ , colSums(OW.Comm) >= 1]
  
  dim(NW.Comm); nrow(NW.LatLong)
  dim(OW.Comm); nrow(OW.LatLong)
  
  ### Subset the community and phylogenetic data for each hemispheric region
  Chiroptera.NW.Comm <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, NW.Comm)$comm
  Chiroptera.NW.Comm <- Chiroptera.NW.Comm[ , colSums(Chiroptera.NW.Comm) >= 1]
  Chiroptera.NW.Tree <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, NW.Comm)$phy
  
  Chiroptera.OW.Comm <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, OW.Comm)$comm
  Chiroptera.OW.Comm <- Chiroptera.OW.Comm[ , colSums(Chiroptera.NW.Comm) >= 1]
  Chiroptera.OW.Tree <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, OW.Comm)$phy
  
  # Calculate the ses.mpd for the communities accounring for the hemispheric sampling
  # pool
  
  #### New World ####
  SES.MPD.Chiroptera.NW <- mod.ses.mpd.query.sf(species.data = Chiroptera.NW.Comm,
                                                phylo.tree = Chiroptera.NW.Tree, 
                                                perms = 999, cores = 20)
  
  #### Old World + Australasian ####
  SES.MPD.Chiroptera.OW <-  mod.ses.mpd.query.sf(species.data = Chiroptera.OW.Comm,
                                                 phylo.tree = Chiroptera.OW.Tree, 
                                                 perms = 999, cores = 20)
  
  # Export results to data.frames
  # Bind both hemispheres in one data.frame
  
  MPD.LatLong.NW.OW <- rbind(SES.MPD.Chiroptera.NW, 
                             SES.MPD.Chiroptera.OW) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Hemispheric sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  # MPD.LatLong.NW.OW %>%
  #   count(SamplingPool)
  
  # dim(MPD.LatLong.NW.OW)
  
  ###########################
  ### Realm sampling pool ###  
  ###########################
  
  chosen_Realms <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_Realm) %>%
    unique() %>%
    filter(!is.na(ID_Realm)) %>%
    filter(!(ID_Realm %in% c("Antarctic")
    )
    ) %>%
    droplevels() %>%
    pull(ID_Realm)
  
  # i = chosen_Realms[1]
  
  # Create a list
  SubRegion.Realm.Comm.Phylo <- list()
  
  for(i in chosen_Realms) {
    
    SubRegion.Realm.LatLong <- subset(world_grid_50km_cat_df, 
                                      ID_Realm == toString(i))
    SubRegion.Realm.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Realm.LatLong), ]
    SubRegion.Realm.Comm. <- SubRegion.Realm.Comm.[ , colSums(SubRegion.Realm.Comm.) >= 1]
    
    ### Subset the data sets
    SubRegion.Realm.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Realm.Comm.)
    SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                    Phylo = SubRegion.Realm.CommPhylo$phy)
  }
  
  
  # j <- ls(SubRegion.Realm.Comm.Phylo)[1]
  
  MPD.LatLong.Realm <- list()
  
  # for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
  for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
    
    message("Computing MPD for: ", j,
            ".\nRegion ",
            which(j == as.factor(ls(SubRegion.Realm.Comm.Phylo))),
            " from ",
            length(as.factor(ls(SubRegion.Realm.Comm.Phylo))),
            ".\n")
    
    # Calculate mpd values for each community
    SES.MPD.Chiroptera.Realm <- mod.ses.mpd.query.sf(species.data = SubRegion.Realm.Comm.Phylo[[j]]$Comm,
                                                     phylo.tree = SubRegion.Realm.Comm.Phylo[[j]]$Phylo, 
                                                     perms = 999, 
                                                     cores = 20)
    
    
    MPD.LatLong.Realm[[paste0(j)]] <- SES.MPD.Chiroptera.Realm
    
    
  }
  
  
  MPD.LatLong.Realm <- MPD.LatLong.Realm %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Realm sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  ###########################
  ### Plate sampling pool ###  
  ###########################
  
  chosen_Plates <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_PlateName) %>%
    group_by(ID_PlateName) %>%
    count() %>%
    filter(n > 1) %>%
    filter(!is.na(ID_PlateName)) %>%
    droplevels() %>%
    pull(ID_PlateName)
  
  # i = chosen_Plates[1]
  
  # Create a list
  SubRegion.Plate.Comm.Phylo <- list()
  
  for(i in chosen_Plates){
    
    SubRegion.Plate.LatLong <- subset(world_grid_50km_cat_df, 
                                      ID_PlateName == toString(i))
    
    SubRegion.Plate.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Plate.LatLong), , drop = FALSE]
    SubRegion.Plate.Comm. <- SubRegion.Plate.Comm.[ , colSums(SubRegion.Plate.Comm.) >= 1, drop = FALSE]
    
    ### Subset the data sets
    SubRegion.Plate.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Plate.Comm.)
    if(is.null(SubRegion.Plate.CommPhylo$phy)){
      message(i, "had zero species present in the phylogeny, and will be excluded from hereafter.")
      vec <- c(vec, i)
    } else {
      SubRegion.Plate.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Plate.CommPhylo$comm,
                                                      Phylo = SubRegion.Plate.CommPhylo$phy)
    }
  }
  
  # j <- ls(SubRegion.Plate.Comm.Phylo)[1]
  
  MPD.LatLong.Plate <- list()
  
  # for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){
  
  for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo))){
    # Calculate mpd values for each community
    
    SES.MPD.Chiroptera.Plate <- mod.ses.mpd.query.sf(species.data = SubRegion.Plate.Comm.Phylo[[j]]$Comm,
                                                     phylo.tree = SubRegion.Plate.Comm.Phylo[[j]]$Phylo, 
                                                     perms = 999, 
                                                     cores = 20)
    
    
    MPD.LatLong.Plate[[paste0(j)]] <- SES.MPD.Chiroptera.Plate
    
    
  }
  
  MPD.LatLong.Plate <- MPD.LatLong.Plate %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Plate sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  
  ###########################
  ### Biome sampling pool ###  
  ###########################
  
  chosen_Biomes <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_Biome_Realm) %>%
    group_by(ID_Biome_Realm) %>%
    count() %>%
    filter(n > 1) %>%
    filter(!is.na(ID_Biome_Realm)) %>%
    droplevels() %>%
    pull(ID_Biome_Realm)
  
  # i = chosen_Biomes[1]
  
  # Create a list
  SubRegion.Biome.Comm.Phylo <- list()
  
  for(i in chosen_Biomes){
    
    SubRegion.Biome.LatLong <- subset(world_grid_50km_cat_df, 
                                      ID_Biome_Realm == toString(i))
    
    SubRegion.Biome.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Biome.LatLong), , drop = FALSE]
    SubRegion.Biome.Comm. <- SubRegion.Biome.Comm.[ , colSums(SubRegion.Biome.Comm.) >= 1, drop = FALSE]
    
    ### Subset the data sets
    SubRegion.Biome.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Biome.Comm.)
    if(is.null(SubRegion.Biome.CommPhylo$phy)){
      message(i, "had zero species present in the phylogeny, and will be excluded from hereafter.")
      vec <- c(vec, i)
    } else {
      SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                      Phylo = SubRegion.Biome.CommPhylo$phy)
    }
  }
  
  # j <- ls(SubRegion.Biome.Comm.Phylo)[1]
  
  MPD.LatLong.Biome <- list()
  
  # for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo)[!ls(SubRegion.Biome.Comm.Phylo) %in% "Antarctic"])){
  
  for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
    
    message("Computing MPD for: ", j,
            ".\nRegion ",
            which(j == as.factor(ls(SubRegion.Biome.Comm.Phylo))),
            " from ",
            length(as.factor(ls(SubRegion.Biome.Comm.Phylo))),
            ".\n")
    
    # Calculate mpd values for each community
    
    SES.MPD.Chiroptera.Biome <- mod.ses.mpd.query.sf(species.data = SubRegion.Biome.Comm.Phylo[[j]]$Comm,
                                                     phylo.tree = SubRegion.Biome.Comm.Phylo[[j]]$Phylo, 
                                                     perms = 999, 
                                                     cores = 20)
    
    
    MPD.LatLong.Biome[[paste0(j)]] <- SES.MPD.Chiroptera.Biome
    
    
  }
  
  MPD.LatLong.Biome <- MPD.LatLong.Biome %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Biome sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  ###########################
  ### Ecoregion sampling pool ###  
  ###########################
  
  chosen_Ecoregions <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_Ecoregion) %>%
    group_by(ID_Ecoregion) %>%
    count() %>%
    filter(n > 1) %>%
    filter(!is.na(ID_Ecoregion)) %>%
    droplevels() %>%
    pull(ID_Ecoregion)
  
  # i = chosen_Ecoregions[1]
  
  # Create a list
  SubRegion.Ecoregion.Comm.Phylo <- list()
  
  for(i in chosen_Ecoregions){
    
    message("Processing: ", i,  
            ".\nRegion ", 
            which(i == chosen_Ecoregions),
            " from ",
            length(chosen_Ecoregions),
            ".\n")
    
    SubRegion.Ecoregion.LatLong <- subset(world_grid_50km_cat_df, 
                                          ID_Ecoregion == toString(i))
    
    SubRegion.Ecoregion.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Ecoregion.LatLong), , drop = FALSE]
    SubRegion.Ecoregion.Comm. <- SubRegion.Ecoregion.Comm.[ , colSums(SubRegion.Ecoregion.Comm.) >= 1, drop = FALSE]
    
    ### Subset the data sets
    SubRegion.Ecoregion.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Ecoregion.Comm.)
    if(is.null(SubRegion.Ecoregion.CommPhylo$phy)){
      message(i, "had zero species present in the phylogeny, and will be excluded from hereafter.")
    } else {
      SubRegion.Ecoregion.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Ecoregion.CommPhylo$comm,
                                                          Phylo = SubRegion.Ecoregion.CommPhylo$phy)
    }
  }
  
  # j <- ls(SubRegion.Ecoregion.Comm.Phylo)[1]
  
  # MPD.LatLong.Ecoregion <- list()
  # 
  # # for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){
  # 
  # for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))){
  #   
  #   message("Computing MPD for: ", j,  
  #           ".\nRegion ", 
  #           which(j == as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))),
  #           " from ",
  #           length(as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))),
  #           ".\n")
  #   
  #   # Calculate mpd values for each community
  #   
  #   SES.MPD.Chiroptera.Ecoregion <- mod.ses.mpd.query.sf(species.data = SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm,
  #                                                        phylo.tree = SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo, 
  #                                                        perms = 999, 
  #                                                        cores = 20)
  #   
  #   
  #   MPD.LatLong.Ecoregion[[paste0(j)]] <- SES.MPD.Chiroptera.Ecoregion
  #   
  #   
  # }
  
  perms <- 999
  
  sfInit(parallel = TRUE, 
         cpus = 35, 
         slaveOutfile = "logAnalysis.txt")
  
  # load libraries in parallel
  
  sfLibrary(PhyloMeasures)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  sfLibrary(picante)
  
  # export objects 
  sfExport('mod.ses.mpd.query', 
           'ses.mpd.query.call',
           'perms',
           'SubRegion.Ecoregion.Comm.Phylo',
           'shuffle_tiplabels',
           'match.phylo.comm.')
  
  # lapply(ls(SubRegion.Ecoregion.Comm.Phylo)[1:3], 
  #        ses.mpd.query.call)
  
  MPD.LatLong.Ecoregion <- sfLapply(ls(SubRegion.Ecoregion.Comm.Phylo), 
                                    ses.mpd.query.call,
                                    Comm.Phylo.list = SubRegion.Ecoregion.Comm.Phylo)
  
  sfStop()
  
  MPD.LatLong.Ecoregion <- MPD.LatLong.Ecoregion %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Ecoregion sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  
  # End of the computation of SES.MPD values across the sampling pool
  # restriction gradient
  
  #################################
  #### Combine all data frames ####
  ################################
  
  # MPD.LatLong.World %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MPD.LatLong.NW.OW %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MPD.LatLong.Realm %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MPD.LatLong.Plate %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MPD.LatLong.Biome %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MPD.LatLong.Ecoregion %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # 
  # MPD.LatLong.World %>%  filter(is.na(SamplingPool))
  # MPD.LatLong.NW.OW %>%  filter(is.na(SamplingPool))
  # MPD.LatLong.Realm %>%  filter(is.na(SamplingPool))
  # MPD.LatLong.Plate %>%  filter(is.na(SamplingPool))
  # MPD.LatLong.Biome %>%  filter(is.na(SamplingPool))
  # MPD.LatLong.Ecoregion %>%  filter(is.na(SamplingPool))
  
  MPD.LatLong.AllScales <- rbind(MPD.LatLong.World,
                                 MPD.LatLong.NW.OW,
                                 MPD.LatLong.Realm,
                                 MPD.LatLong.Plate,
                                 MPD.LatLong.Biome,
                                 MPD.LatLong.Ecoregion)
  
  
  
  # Begin of the computation of SES.MNTD values across the sampling pool
  # restriction gradient
  
  ########################
  ### Global sampling ####
  ########################
  
  # Calculate mntd values for each community
  SES.MNTD.Chiroptera.World <- mod.ses.mntd.query.sf(species.data = Chiroptera.FaurSven.comm.sampled,
                                                     phylo.tree = Chiroptera.FaurSven.tree.sampled, 
                                                     perms = 999, cores = 20)
  # Add everything to a data.frame
  SES.MNTD.Chiroptera.World$ID <- as.numeric(SES.MNTD.Chiroptera.World$ID)
  
  MNTD.LatLong.World <- right_join(SES.MNTD.Chiroptera.World,
                                   world_grid_50km_cat_df,
                                   by = c("ID" = "ID")) %>%
    mutate(SamplingPool = rep("Global sampling", 
                              nrow(SES.MNTD.Chiroptera.World))
    )
  
  ##########################################
  #### New World and Old World sampling ####
  ##########################################
  
  # Defining the New world and the Old World
  NW.LatLong <- subset(world_grid_50km_cat_df, 
                       ID_Realm == "Nearctic" | ID_Realm == "Neotropical")
  
  OW.LatLong <- subset(world_grid_50km_cat_df, 
                       (ID_Realm %in% c("Palearctic", "Australasian",
                                        "Afrotropical", "Indomalay",
                                        "Oceanic")
                       )
  )
  
  nrow(NW.LatLong) + nrow(OW.LatLong); nrow(world_grid_50km_cat_df) # Checking the number of rows
  
  # Extracting NW communities
  NW.Comm <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(NW.LatLong), ]
  NW.Comm <- NW.Comm[ , colSums(NW.Comm) >= 1]
  
  # Extracting OW communities
  OW.Comm <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(OW.LatLong), ]
  OW.Comm <- OW.Comm[ , colSums(OW.Comm) >= 1]
  
  dim(NW.Comm); nrow(NW.LatLong)
  dim(OW.Comm); nrow(OW.LatLong)
  
  ### Subset the community and phylogenetic data for each hemispheric region
  Chiroptera.NW.Comm <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, NW.Comm)$comm
  Chiroptera.NW.Comm <- Chiroptera.NW.Comm[ , colSums(Chiroptera.NW.Comm) >= 1]
  Chiroptera.NW.Tree <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, NW.Comm)$phy
  
  Chiroptera.OW.Comm <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, OW.Comm)$comm
  Chiroptera.OW.Comm <- Chiroptera.OW.Comm[ , colSums(Chiroptera.NW.Comm) >= 1]
  Chiroptera.OW.Tree <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, OW.Comm)$phy
  
  # Calculate the ses.mntd for the communities accounting for the hemispheric sampling
  # pool
  
  #### New World ####
  SES.MNTD.Chiroptera.NW <- mod.ses.mntd.query.sf(species.data = Chiroptera.NW.Comm,
                                                  phylo.tree = Chiroptera.NW.Tree, 
                                                  perms = 999, cores = 20)
  
  #### Old World ####
  SES.MNTD.Chiroptera.OW <-  mod.ses.mntd.query.sf(species.data = Chiroptera.OW.Comm,
                                                   phylo.tree = Chiroptera.OW.Tree, 
                                                   perms = 999, cores = 20)
  
  # Export results to data.frames
  # Bind both hemispheres in one data.frame
  
  MNTD.LatLong.NW.OW <- rbind(SES.MNTD.Chiroptera.NW, 
                              SES.MNTD.Chiroptera.OW) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Hemispheric sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  # MNTD.LatLong.NW.OW %>%
  #   count(SamplingPool)
  
  # dim(MNTD.LatLong.NW.OW)
  
  ###########################
  ### Realm sampling pool ###  
  ###########################
  
  chosen_Realms <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_Realm) %>%
    unique() %>%
    filter(!is.na(ID_Realm)) %>%
    filter(!(ID_Realm %in% c("Antarctic", "Oceanic"))) %>%
    droplevels() %>%
    pull(ID_Realm)
  
  # i = chosen_Realms[1]
  
  # Create a list
  SubRegion.Realm.Comm.Phylo <- list()
  
  for(i in chosen_Realms) {
    
    SubRegion.Realm.LatLong <- subset(world_grid_50km_cat_df, 
                                      ID_Realm == toString(i))
    SubRegion.Realm.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Realm.LatLong), ]
    SubRegion.Realm.Comm. <- SubRegion.Realm.Comm.[ , colSums(SubRegion.Realm.Comm.) >= 1]
    
    ### Subset the data sets
    SubRegion.Realm.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Realm.Comm.)
    SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                    Phylo = SubRegion.Realm.CommPhylo$phy)
  }
  
  # j <- ls(SubRegion.Realm.Comm.Phylo)[1]
  
  MNTD.LatLong.Realm <- list()
  
  # for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo)[!ls(SubRegion.Realm.Comm.Phylo) %in% "Antarctic"])){
  for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
    
    message("Computing MNTD for: ", j,
            ".\nRegion ",
            which(j == as.factor(ls(SubRegion.Realm.Comm.Phylo))),
            " from ",
            length(as.factor(ls(SubRegion.Realm.Comm.Phylo))),
            ".\n")
    # Calculate mntd values for each community
    
    SES.MNTD.Chiroptera.Realm <- mod.ses.mntd.query.sf(species.data = SubRegion.Realm.Comm.Phylo[[j]]$Comm,
                                                       phylo.tree = SubRegion.Realm.Comm.Phylo[[j]]$Phylo, 
                                                       perms = 999, 
                                                       cores = 20)
    
    
    MNTD.LatLong.Realm[[paste0(j)]] <- SES.MNTD.Chiroptera.Realm
    
    
  }
  
  
  MNTD.LatLong.Realm <- MNTD.LatLong.Realm %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Realm sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  ###########################
  ### Plate sampling pool ###  
  ###########################
  
  chosen_Plates <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_PlateName) %>%
    group_by(ID_PlateName) %>%
    count() %>%
    filter(n > 1) %>%
    filter(!is.na(ID_PlateName)) %>%
    droplevels() %>%
    pull(ID_PlateName)
  
  # i = chosen_Plates[1]
  
  # Create a list
  SubRegion.Plate.Comm.Phylo <- list()
  
  for(i in chosen_Plates){
    
    SubRegion.Plate.LatLong <- subset(world_grid_50km_cat_df, 
                                      ID_PlateName == toString(i))
    
    SubRegion.Plate.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Plate.LatLong), , drop = FALSE]
    SubRegion.Plate.Comm. <- SubRegion.Plate.Comm.[ , colSums(SubRegion.Plate.Comm.) >= 1, drop = FALSE]
    
    ### Subset the data sets
    SubRegion.Plate.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Plate.Comm.)
    if(is.null(SubRegion.Plate.CommPhylo$phy)){
      message(i, "had zero species present in the phylogeny, and will be excluded from hereafter.")
      vec <- c(vec, i)
    } else {
      SubRegion.Plate.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Plate.CommPhylo$comm,
                                                      Phylo = SubRegion.Plate.CommPhylo$phy)
    }
  }
  
  # j <- ls(SubRegion.Plate.Comm.Phylo)[1]
  
  MNTD.LatLong.Plate <- list()
  
  # for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo)[!ls(SubRegion.Plate.Comm.Phylo) %in% "Antarctic"])){
  
  for(j in as.factor(ls(SubRegion.Plate.Comm.Phylo))){
    # Calculate mntd values for each community
    
    SES.MNTD.Chiroptera.Plate <- mod.ses.mntd.query.sf(species.data = SubRegion.Plate.Comm.Phylo[[j]]$Comm,
                                                       phylo.tree = SubRegion.Plate.Comm.Phylo[[j]]$Phylo, 
                                                       perms = 999, 
                                                       cores = 20)
    
    
    MNTD.LatLong.Plate[[paste0(j)]] <- SES.MNTD.Chiroptera.Plate
    
    
  }
  
  MNTD.LatLong.Plate <- MNTD.LatLong.Plate %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Plate sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  
  ###########################
  ### Biome sampling pool ###  
  ###########################
  
  chosen_Biomes <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_Biome_Realm) %>%
    group_by(ID_Biome_Realm) %>%
    count() %>%
    filter(n > 1) %>%
    filter(!is.na(ID_Biome_Realm)) %>%
    droplevels() %>%
    pull(ID_Biome_Realm)
  
  # i = chosen_Biomes[1]
  
  # Create a list
  SubRegion.Biome.Comm.Phylo <- list()
  
  for(i in chosen_Biomes){
    
    SubRegion.Biome.LatLong <- subset(world_grid_50km_cat_df, 
                                      ID_Biome_Realm == toString(i))
    
    SubRegion.Biome.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Biome.LatLong), , drop = FALSE]
    SubRegion.Biome.Comm. <- SubRegion.Biome.Comm.[ , colSums(SubRegion.Biome.Comm.) >= 1, drop = FALSE]
    
    ### Subset the data sets
    SubRegion.Biome.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Biome.Comm.)
    if(is.null(SubRegion.Biome.CommPhylo$phy)){
      message(i, "had zero species present in the phylogeny, and will be excluded from hereafter.")
      vec <- c(vec, i)
    } else {
      SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                      Phylo = SubRegion.Biome.CommPhylo$phy)
    }
  }
  
  # j <- ls(SubRegion.Biome.Comm.Phylo)[1]
  
  MNTD.LatLong.Biome <- list()
  
  # for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo)[!ls(SubRegion.Biome.Comm.Phylo) %in% "Antarctic"])){
  
  for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
    
    message("Computing MNTD for: ", j,
            ".\nRegion ",
            which(j == as.factor(ls(SubRegion.Biome.Comm.Phylo))),
            " from ",
            length(as.factor(ls(SubRegion.Biome.Comm.Phylo))),
            ".\n")
    
    # Calculate mntd values for each community
    
    SES.MNTD.Chiroptera.Biome <- mod.ses.mntd.query.sf(species.data = SubRegion.Biome.Comm.Phylo[[j]]$Comm,
                                                       phylo.tree = SubRegion.Biome.Comm.Phylo[[j]]$Phylo, 
                                                       perms = 999, 
                                                       cores = 20)
    
    
    MNTD.LatLong.Biome[[paste0(j)]] <- SES.MNTD.Chiroptera.Biome
    
    
  }
  
  MNTD.LatLong.Biome <- MNTD.LatLong.Biome %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Biome sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  ###########################
  ### Ecoregion sampling pool ###  
  ###########################
  
  chosen_Ecoregions <- world_grid_50km_cat_df %>%
    as_tibble() %>%
    select(ID_Ecoregion) %>%
    group_by(ID_Ecoregion) %>%
    count() %>%
    filter(n > 1) %>%
    filter(!is.na(ID_Ecoregion)) %>%
    droplevels() %>%
    pull(ID_Ecoregion)
  
  # i = chosen_Ecoregions[1]
  
  # Create a list
  SubRegion.Ecoregion.Comm.Phylo <- list()
  
  for(i in chosen_Ecoregions){
    
    message("Processing: ", i,  
            ".\nRegion ", 
            which(i == chosen_Ecoregions),
            " from ",
            length(chosen_Ecoregions),
            ".\n")
    
    SubRegion.Ecoregion.LatLong <- subset(world_grid_50km_cat_df, 
                                          ID_Ecoregion == toString(i))
    
    SubRegion.Ecoregion.Comm. <- Chiroptera.FaurSven.comm.sampled[rownames(Chiroptera.FaurSven.comm.sampled) %in% row.names(SubRegion.Ecoregion.LatLong), , drop = FALSE]
    SubRegion.Ecoregion.Comm. <- SubRegion.Ecoregion.Comm.[ , colSums(SubRegion.Ecoregion.Comm.) >= 1, drop = FALSE]
    
    ### Subset the data sets
    SubRegion.Ecoregion.CommPhylo <- match.phylo.comm.(Chiroptera.FaurSven.tree.sampled, SubRegion.Ecoregion.Comm.)
    if(is.null(SubRegion.Ecoregion.CommPhylo$phy)){
      message(i, "had zero species present in the phylogeny, and will be excluded from hereafter.")
      vec <- c(vec, i)
    } else {
      SubRegion.Ecoregion.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Ecoregion.CommPhylo$comm,
                                                          Phylo = SubRegion.Ecoregion.CommPhylo$phy)
    }
  }
  
  # j <- ls(SubRegion.Ecoregion.Comm.Phylo)[1]
  
  # MNTD.LatLong.Ecoregion <- list()
  # 
  # # for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo)[!ls(SubRegion.Ecoregion.Comm.Phylo) %in% "Antarctic"])){
  # 
  # for(j in as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))){
  #   
  #   message("Computing MNTD for: ", j,  
  #           ".\nRegion ", 
  #           which(j == as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))),
  #           " from ",
  #           length(as.factor(ls(SubRegion.Ecoregion.Comm.Phylo))),
  #           ".\n")
  #   
  #   # Calculate mntd values for each community
  #   
  #   SES.MNTD.Chiroptera.Ecoregion <- mod.ses.mntd.query.sf(species.data = SubRegion.Ecoregion.Comm.Phylo[[j]]$Comm,
  #                                                        phylo.tree = SubRegion.Ecoregion.Comm.Phylo[[j]]$Phylo, 
  #                                                        perms = 999, 
  #                                                        cores = 20)
  #   
  #   
  #   MNTD.LatLong.Ecoregion[[paste0(j)]] <- SES.MNTD.Chiroptera.Ecoregion
  #   
  #   
  # }
  
  perms <- 999
  
  sfInit(parallel = TRUE, 
         cpus = 35, 
         slaveOutfile = "logAnalysis.txt")
  
  # load libraries in parallel
  
  sfLibrary(PhyloMeasures)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  sfLibrary(picante)
  
  # export objects 
  sfExport('mod.ses.mntd.query', 
           'ses.mntd.query.call',
           'perms',
           'SubRegion.Ecoregion.Comm.Phylo',
           'shuffle_tiplabels')
  
  # lapply(ls(SubRegion.Ecoregion.Comm.Phylo)[1:3], 
  #        ses.mntd.query.call)
  
  MNTD.LatLong.Ecoregion <- sfLapply(ls(SubRegion.Ecoregion.Comm.Phylo), 
                                     ses.mntd.query.call,
                                     Comm.Phylo.list = SubRegion.Ecoregion.Comm.Phylo)
  
  sfStop()
  
  MNTD.LatLong.Ecoregion <- MNTD.LatLong.Ecoregion %>%
    list.rbind() %>%
    set_rownames(.$ID) %>%
    mutate_at("ID", as.numeric) %>%
    mutate(SamplingPool = "Ecoregion sampling") %>%
    full_join(world_grid_50km_cat_df,
              by = c("ID" = "ID")) %>%
    set_rownames(.$ID)
  
  
  # End of the computation of SES.MPD values across the sampling pool
  # restriction gradient
  
  #################################
  #### Combine all data frames ####
  ################################
  
  # MNTD.LatLong.World %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MNTD.LatLong.NW.OW %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MNTD.LatLong.Realm %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MNTD.LatLong.Plate %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MNTD.LatLong.Biome %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # MNTD.LatLong.Ecoregion %>%  group_by(SamplingPool) %>% count(SamplingPool)
  # 
  # MNTD.LatLong.World %>%  filter(is.na(SamplingPool))
  # MNTD.LatLong.NW.OW %>%  filter(is.na(SamplingPool))
  # MNTD.LatLong.Realm %>%  filter(is.na(SamplingPool))
  # MNTD.LatLong.Plate %>%  filter(is.na(SamplingPool))
  # MNTD.LatLong.Biome %>%  filter(is.na(SamplingPool))
  # MNTD.LatLong.Ecoregion %>%  filter(is.na(SamplingPool))
  
  MNTD.LatLong.AllScales <- rbind(MNTD.LatLong.World,
                                  MNTD.LatLong.NW.OW,
                                  MNTD.LatLong.Realm,
                                  MNTD.LatLong.Plate,
                                  MNTD.LatLong.Biome,
                                  MNTD.LatLong.Ecoregion)
  
  
  # Begin of the computation of SES.MPD values across the sampling pool
  # restriction gradient
  
  #########################################
  ### Merge data into a single data set ###
  #########################################
  
  MPD.MNTD.LatLong.AllScales <- MPD.LatLong.AllScales %>%
    unite(ID_SamplingPool, 
          ID, 
          SamplingPool, 
          remove = FALSE) %>%
    filter(!is.na(SamplingPool)) %>%
    left_join(MNTD.LatLong.AllScales %>% 
                unite(ID_SamplingPool, 
                      ID, 
                      SamplingPool, 
                      remove = FALSE) %>%
                filter(!is.na(SamplingPool)) %>%
                select(-c("ntaxa", "ID", "perms", "ntaxa.pool.size",
                          LONB:SamplingPool)),
              by = "ID_SamplingPool")
  
  MPD.MNTD.LatLong.AllScales$phylo.tree.pos <- phylo.sampled.pos
  
  filename.chosen <- paste0("MPD.MNTD.LatLong.AllScales",
                            "_phylo_",
                            phylo.sampled.pos,
                            ".csv")
  
  data.table::fwrite(MPD.MNTD.LatLong.AllScales,
                     file = paste0("data/matrices/phyoSampling_mpd_mntd/", 
                                   filename.chosen),
                     row.names = TRUE, 
                     quote = TRUE)
  
}

ending.time <- Sys.time()
