######################################################################################
#### Code to merge, represent the MPD and MNTD, and test for the effects of       ####
#### sampling pool restriction in the phylogenetic structuring of communities     ####
#                                                                                    #
# Description: This code merges the MPD.LatLong.Env.AllScales and                    #
# MNTD.LatLong.Env.AllScales data sets into a single one. It then represents NRI     #
# and NTI across the gradient of sampling pool restrictions. It finally tests        #
# the effects of sampling pool (or species pool) restriction on the phylogenetic     #
# community structure of bats.                                                       #
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

MPD.MNTD.LatLong.AllScales <- MPD.LatLong.Env.AllScales %>%
  dplyr::select(-c(bio_1:bio_19)) %>% # You may remove this part.
  dplyr::rename(mpd.runs = runs,
                mpd.ntaxa = ntaxa) %>%
  unite(ID_SamplingPool, ID, SamplingPool, remove = FALSE) %>%
  left_join(MNTD.LatLong.Env.AllScales %>%
              dplyr::select(-c(bio_1:bio_19)) %>%
              dplyr::rename(mntd.runs = runs,
                            mntd.ntaxa = ntaxa) %>% 
              unite(ID_SamplingPool, ID, SamplingPool, remove = FALSE),
            by = "ID_SamplingPool")

# Verify

MPD.MNTD.LatLong.AllScales <- MPD.LatLong.Env.AllScales %>%
  dplyr::select(-c(bio_1:bio_19)) %>%
  dplyr::rename(mpd.runs = runs,
                mpd.ntaxa = ntaxa) %>%
  unite(ID_SamplingPool, ID, SamplingPool, remove = FALSE) %>%
  left_join(MNTD.LatLong.Env.AllScales %>%
              dplyr::select(-c(Longitude:bio_19)) %>%
              dplyr::rename(mntd.runs = runs,
                            mntd.ntaxa = ntaxa) %>% 
              unite(ID_SamplingPool, ID, SamplingPool, 
                    remove = TRUE),
            by = "ID_SamplingPool")


# Refactor data

MPD.MNTD.LatLong.AllScales$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales$ID_Realm, 
                                              levels=c('Neotropical',
                                                       'Nearctic',
                                                       'Afrotropical',
                                                       'Palearctic', 
                                                       'Indomalay', 
                                                       'Australasian'))

MPD.MNTD.LatLong.AllScales$SamplingPool <-  factor(MPD.MNTD.LatLong.AllScales$SamplingPool,  
                                                   levels = c("Global sampling",
                                                              "Hemispheric sampling",
                                                              "Realm sampling",
                                                              "Plate sampling",
                                                              "Biome sampling",
                                                              "Ecoregion sampling"))

############################################################################ 
### Representing the phylogenetic community structure of bats across the ###
### gradient of sampling pool (or species pool) restriction              ###
############################################################################

MPD.LatLong.Env.AllScales$ID_Biome_Acronym = factor(MPD.LatLong.Env.AllScales$ID_Biome, 
                                                    labels = abbreviate(gsub("_", 
                                                                             " ",
                                                                             levels(MPD.LatLong.Env.AllScales$ID_Biome))))

# Representing NRI wrapping across Realms

(NRI.Realm.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales,
                                    is.na(ID_Realm) == FALSE),
                             aes(x = SamplingPool, y = NRI)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    facet_wrap(~ID_Realm,
               nrow = 1,
               strip.position = "bottom") +
    scale_fill_uchicago() +
    # scale_fill_nord("baie_mouton", reverse = TRUE) +
    # scale_fill_okabeito(reverse = TRUE) +
    # scale_fill_brewer(palette = "Dark2") +
    # scale_fill_viridis(discrete = TRUE,
    #                    name="Sampling Pool") +
    scale_y_continuous(breaks = round(
      #pretty(MPD.LatLong.Env.AllScales$NRI, n = 7)
      c(min(MPD.MNTD.LatLong.AllScales$NRI, na.rm = T),
        0,
        5,
        10,
        15,
        max(MPD.MNTD.LatLong.AllScales$NRI, na.rm = T)))
    ) +
    geom_hline(yintercept = 0,
               alpha = 0.4) +
    labs(x="") +
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

# If necessary, save figure
# ggsave(filename = "NRI.Realm.boxplot.png", 
#        dpi = 300, 
#        width = 18.5, height = 5, 
#        units = "in")

# Representing NTI wrapping across Realms

(NTI.Realm.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales,
                                    is.na(ID_Realm) == FALSE),
                             aes(x = SamplingPool, y = NTI)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    facet_wrap(~ID_Realm,
               nrow = 1,
               strip.position = "bottom") +
    scale_fill_uchicago() +
    # scale_fill_nord("baie_mouton", reverse = TRUE) +
    # scale_fill_okabeito(reverse = TRUE) +
    # scale_fill_brewer(palette = "Dark2") +
    # scale_fill_viridis(discrete = TRUE,
    #                    name="Sampling Pool") +
    # scale_y_continuous(breaks = pretty(MPD.LatLong.Env.AllScales$NTI, n = 6)
    #   # round(c(min(MPD.MNTD.LatLong.AllScales$NTI, na.rm = T),
    #   #   0,
    #   #   5,
    #   #   10,
    #   #   15,
  #   #   max(MPD.MNTD.LatLong.AllScales$NTI, na.rm = T))
  # )) +
  geom_hline(yintercept = 0,
             alpha = 0.4) +
    labs(x="") +
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

(fig.NRI.NTI.Realm.boxplot <- ggarrange(
  NRI.Realm.boxplot +
    theme(axis.title.x = element_blank()), 
  NTI.Realm.boxplot,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  widths = c(1, 1),
  heights = c(1, 1.1),
  align = "v",
  common.legend = TRUE,
  legend = "bottom"
)
)

### Exporting the box-plot to PNG ###

# This is the Figure 2 in the manuscript

ggsave(filename = "figures/fig.NRI.NTI.Realm.boxplot.png", 
       dpi = 300, 
       width = 18.5, height = 10, 
       units = "in")

##########################################################
####
#

# head(MPD.MNTD.LatLong.AllScales)


MPD.MNTD.LatLong.AllScales %>%
  select(ID, NRI, SamplingPool, ID_Realm) %>%
  filter(ID_Realm == "Neotropical") %>%
  spread(SamplingPool, NRI) %>%
  drop_na() %>%
  gather(key = "SamplingPool", value = "NRI", -ID_Realm, -ID)

friedman.test(NRI ~ ID_Realm | ID, data = MPD.MNTD.LatLong.AllScales %>%
                select(ID, NRI, SamplingPool, ID_Realm) %>%
                filter(ID_Realm == "Neotropical") %>%
                spread(SamplingPool, NRI) %>%
                drop_na() %>%
                gather(key = "SamplingPool", value = "NRI", -ID_Realm, -ID) %>%
                mutate(ID = as.factor(ID)) %>%
                as.data.frame()) 


tabulation <- table(MPD.MNTD.LatLong.AllScales %>%
                      select(ID, NRI, SamplingPool, ID_Realm) %>%
                      filter(ID_Realm == "Neotropical") %>%
                      spread(SamplingPool, NRI) %>%
                      drop_na() %>%
                      gather(key = "SamplingPool", value = "NRI", -ID_Realm, -ID) %>%
                      as.data.frame() %>%
                      select(ID, SamplingPool)) %>% 
  as.data.frame() %>%
  filter_if(is.numeric, all_vars((.) == 1))

rowSums(tabulation)



NRI.SamplingPool.NT <- MPD.MNTD.LatLong.AllScales %>%
  select(ID, NRI, SamplingPool, ID_Realm) %>%
  filter(ID_Realm == "Neotropical") %>%
  spread(SamplingPool, NRI) %>%
  drop_na() %>%
  gather(key = "SamplingPool", value = "NRI", -ID_Realm, -ID) %>%
  mutate(ID = as.factor(ID)) %>%
  as.data.frame()

library(WRS2)

r <- levels(MPD.MNTD.LatLong.AllScales$ID_Realm)

rmmcp.NRI.SamplingPool.summary.r <- data.frame()

# ri <- r[1]

for(ri in r){
  
  NRI.SamplingPool.ri <- MPD.MNTD.LatLong.AllScales %>%
    select(ID, NRI, SamplingPool, ID_Realm) %>%
    filter(ID_Realm == paste(ri)) %>%
    spread(SamplingPool, NRI) %>%
    drop_na() %>%
    gather(key = "SamplingPool", value = "NRI", -ID_Realm, -ID) %>%
    mutate(ID = as.factor(ID)) %>%
    as.data.frame()
  
  rmanova.NRI.SamplingPool.ri <- WRS2::rmanova(NRI.SamplingPool.ri$NRI, 
                                         NRI.SamplingPool.ri$SamplingPool, 
                                         NRI.SamplingPool.ri$ID)
  
  
  rmmcp.NRI.SamplingPool.ri <- WRS2::rmmcp(NRI.SamplingPool.ri$NRI, 
                                     NRI.SamplingPool.ri$SamplingPool, 
                                     NRI.SamplingPool.ri$ID)
  
  rmmcp.NRI.SamplingPool.summary.ri <- rmmcp.NRI.SamplingPool.ri$comp %>%
    as.data.frame() %>%
    rename(Group_1 = "Group") %>% 
    rename(Group_2 = "Group") %>%
    mutate_at(.vars = vars(Group_1:Group_2), 
              .funs = function(x) recode(x, 
                                         `1` = rmmcp.NRI.SamplingPool.ri$fnames[1], 
                                         `2` = rmmcp.NRI.SamplingPool.ri$fnames[2],
                                         `3` = rmmcp.NRI.SamplingPool.ri$fnames[3],
                                         `4` = rmmcp.NRI.SamplingPool.ri$fnames[4],
                                         `5` = rmmcp.NRI.SamplingPool.ri$fnames[5],
                                         `6` = rmmcp.NRI.SamplingPool.ri$fnames[6])) %>%
    mutate(Comparison = paste(Group_1, "vs.", Group_2), .before = 1) %>%
    mutate(ID_Realm = paste(ri), .before = 1) %>%
    mutate(sig = ifelse(p.value >= 0.01 & p.value < 0.05, "*", 
                        ifelse(p.value < 0.01 & p.value >= 0.001, "**",
                               ifelse(p.value < 0.001, "***", "")))) # %>%
  # select(-Group_1, -Group_2)
  
  
  rmmcp.NRI.SamplingPool.summary.r <- bind_rows(rmmcp.NRI.SamplingPool.summary.r, 
                                                rmmcp.NRI.SamplingPool.summary.ri)
}

rmmcp.NRI.SamplingPool.summary.r


###

rmmcp.NTI.SamplingPool.summary.r <- data.frame()

for(ri in r){
  
  NTI.SamplingPool.ri <- MPD.MNTD.LatLong.AllScales %>%
    select(ID, NTI, SamplingPool, ID_Realm) %>%
    filter(ID_Realm == paste(ri)) %>%
    spread(SamplingPool, NTI) %>%
    drop_na() %>%
    gather(key = "SamplingPool", value = "NTI", -ID_Realm, -ID) %>%
    mutate(ID = as.factor(ID)) %>%
    as.data.frame()
  
  rmanova.NTI.SamplingPool.ri <- WRS2::rmanova(NTI.SamplingPool.ri$NTI, 
                                         NTI.SamplingPool.ri$SamplingPool, 
                                         NTI.SamplingPool.ri$ID)
  
  
  rmmcp.NTI.SamplingPool.ri <- WRS2::rmmcp(NTI.SamplingPool.ri$NTI, 
                                     NTI.SamplingPool.ri$SamplingPool, 
                                     NTI.SamplingPool.ri$ID)
  
  rmmcp.NTI.SamplingPool.summary.ri <- rmmcp.NTI.SamplingPool.ri$comp %>%
    as.data.frame() %>%
    rename(Group_1 = "Group") %>% 
    rename(Group_2 = "Group") %>%
    mutate_at(.vars = vars(Group_1:Group_2), 
              .funs = function(x) recode(x, 
                                         `1` = rmmcp.NTI.SamplingPool.ri$fnames[1], 
                                         `2` = rmmcp.NTI.SamplingPool.ri$fnames[2],
                                         `3` = rmmcp.NTI.SamplingPool.ri$fnames[3],
                                         `4` = rmmcp.NTI.SamplingPool.ri$fnames[4],
                                         `5` = rmmcp.NTI.SamplingPool.ri$fnames[5],
                                         `6` = rmmcp.NTI.SamplingPool.ri$fnames[6])) %>%
    mutate(Comparison = paste(Group_1, "vs.", Group_2), .before = 1) %>%
    mutate(ID_Realm = paste(ri), .before = 1) %>%
    mutate(sig = ifelse(p.value >= 0.01 & p.value < 0.05, "*", 
                        ifelse(p.value < 0.01 & p.value >= 0.001, "**",
                               ifelse(p.value < 0.001, "***", "")))) # %>%
  # select(-Group_1, -Group_2)
  
  rmmcp.NTI.SamplingPool.summary.r <- bind_rows(rmmcp.NTI.SamplingPool.summary.r, 
                                                rmmcp.NTI.SamplingPool.summary.ri)
}

# Write table with post-hoc tests

kable(rmmcp.NRI.SamplingPool.summary.r)

write.table(rmmcp.NRI.SamplingPool.summary.r)

write.table(rmmcp.NTI.SamplingPool.summary.r)
