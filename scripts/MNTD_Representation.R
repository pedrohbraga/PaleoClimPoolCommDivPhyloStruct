###########################################################
#### Statistical analyses and graphical representation ####
###########################################################

#########################################################

## How does NTI across different biogeographical realms?
# NTI ~ Realm

# Linear Model
## Global sampling

MNTD.Global <- MNTD.LatLong.Env.AllScales %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(!is.na(ID_Realm))

# Assumption testing
car::leveneTest(NTI ~ ID_Realm, data = MNTD.Global, na.rm = TRUE)
# Assumptions violated

Bats.NTI.Realm <- aov(NTI ~ ID_Realm, data = MNTD.Global)
summary(Bats.NTI.Realm)

TukeyHSD(x=Bats.NTI.Realm, 
         'ID_Realm', 
         conf.level=0.95)

# Non-parametric alternative

kruskal.test(NTI ~ ID_Realm, data = MNTD.Global)

pairwise.wilcox.test(MNTD.Global$NTI, MNTD.Global$ID_Realm,
                     p.adjust.method = "BH")

# Linear Mixed Model

Bats.NTI.Realm.LMM <- lmer(NTI ~ ID_Realm + (1|SamplingPool), 
                           data = MNTD.LatLong.Env.AllScales)
plot(Bats.NTI.Realm.LMM)

summary(Bats.NTI.Realm.LMM)


#################################################################
NTI.SamplingPool.ID_Realm <- list(filter(MNTD.LatLong.Env.AllScales, 
                                         SamplingPool == "Global sampling")[, c(1, 6, 36)],
                                  filter(MNTD.LatLong.Env.AllScales, 
                                         SamplingPool == "Hemispheric sampling")[, c(1, 6, 36)],
                                  filter(MNTD.LatLong.Env.AllScales, 
                                         SamplingPool == "Realm sampling")[, c(1, 6, 36)],
                                  filter(MNTD.LatLong.Env.AllScales, 
                                         SamplingPool == "Biome sampling")[, c(1, 6, 36)]) %>%
  Reduce(function(dtf1, dtf2) left_join(dtf1,dtf2,by="ID"), .)

colnames(NTI.SamplingPool.ID_Realm) <- c("ID",
                                         "ID_Realm", "NTI.Global",
                                         "ID_Realm.", "NTI.Hemispheric",
                                         "ID_Realm..", "NTI.Realm",
                                         "ID_Realm...", "NTI.Biome")

diff.NTI.SamplingPool <- data.frame(ID = NTI.SamplingPool.ID_Realm$ID,
                                    ID_Realm = NTI.SamplingPool.ID_Realm$ID_Realm,
                                    diff.NTI.Global.Hemispheric = NTI.SamplingPool.ID_Realm$NTI.Global - NTI.SamplingPool.ID_Realm$NTI.Hemispheric,
                                    diff.NTI.Hemispheric.Realm = NTI.SamplingPool.ID_Realm$NTI.Hemispheric - NTI.SamplingPool.ID_Realm$NTI.Realm,
                                    diff.NTI.Realm.Biome = NTI.SamplingPool.ID_Realm$NTI.Realm - NTI.SamplingPool.ID_Realm$NTI.Biome)

diff.NTI.SamplingPool <- diff.NTI.SamplingPool %>%
  gather(diff.SamplingPool, diff.NTI, -ID, -ID_Realm)



diff.NTI.SamplingPool.boxplot <- ggplot(filter(diff.NTI.SamplingPool,
                                               is.na(ID_Realm) == FALSE),
                                        aes(x = diff.SamplingPool, y = diff.NTI)) +
  geom_boxplot(aes(fill = diff.SamplingPool)) +
  facet_wrap(~ID_Realm,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(breaks = round(
    #pretty(MNTD.LatLong.Env.AllScales$NTI, n = 7)
    c(min(diff.NTI.SamplingPool$diff.NTI, na.rm = T),
      0,
      5,
      10,
      15,
      max(diff.NTI.SamplingPool$diff.NTI, na.rm = T)))
  ) +
  geom_hline(yintercept=0,
             alpha = 0.4) +
  labs(x="",
       y = expression(Delta*"NTI")) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE,
                     name="Sampling Pool Difference") +
  theme(strip.background=element_rect(fill = "white",
                                      linetype = NULL,
                                      color = "white"),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = ""
  )

diff.NTI.SamplingPool.boxplot

ggsave(filename = "diff.NTI.SamplingPool.boxplot.png", 
       dpi = 300, 
       width = 17.5, height = 5, 
       units = "in")

# Effects of sampling pool scaling on NTI within biogeographic realms. Within
# each panel, the left box plot shows the change in NTI between global and
# hemispheric sampling pools (i.e., NTIglobal – NTIhemispheric). The right box
# plot shows the change in NTI between hemispheric and realm sampling pools
# (NTIhemispheric – NTIrealm). Values > 0 indicate stronger phylogenetic
# clustering with the geographically more extensive sampling pool.

########################################################


MNTD.LatLong.Env.AllScales$ID_Biome_Acronym = factor(MNTD.LatLong.Env.AllScales$ID_Biome, 
                                                    labels = abbreviate(gsub("_", 
                                                                             " ",
                                                                             levels(MNTD.LatLong.Env.AllScales$ID_Biome))))

NTI.Realm.boxplot <- ggplot(filter(MNTD.LatLong.Env.AllScales,
                                   is.na(ID_Realm) == FALSE),
                            aes(x = SamplingPool, y = NTI)) +
  geom_boxplot(aes(fill = SamplingPool)) +
  facet_wrap(~ID_Realm,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(breaks = round(
    #pretty(MNTD.LatLong.Env.AllScales$NTI, n = 7)
    c(min(MNTD.LatLong.Env.AllScales$NTI, na.rm = T),
      0,
      5,
      10,
      15,
      max(MNTD.LatLong.Env.AllScales$NTI, na.rm = T)))
  ) +
  geom_hline(yintercept=0,
             alpha = 0.4) +
  labs(x="") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE,
                     name="Sampling Pool") +
  theme(strip.background=element_rect(fill = "white",
                                      linetype = NULL,
                                      color = "white"),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
  )

NTI.Realm.boxplot

ggsave(filename = "NTI.Realm.boxplot.png", 
       dpi = 300, 
       width = 18.5, height = 5, 
       units = "in")

##

NTI.Biome.boxplot <- ggplot(filter(MNTD.LatLong.Env.AllScales,
                                   is.na(ID_Biome) == FALSE),
                            aes(x = SamplingPool, y = NTI)) +
  geom_boxplot(aes(fill = SamplingPool)) +
  facet_wrap(~ID_Biome_Acronym,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(breaks = round(
    #pretty(MNTD.LatLong.Env.AllScales$NTI, n = 7)
    c(min(MNTD.LatLong.Env.AllScales$NTI, na.rm = T),
      0,
      5,
      10,
      15,
      max(MNTD.LatLong.Env.AllScales$NTI, na.rm = T)))
  ) +
  geom_hline(yintercept=0,
             alpha = 0.4) +
  labs(x="") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE,
                     name="Sampling Pool") +
  theme(strip.background=element_rect(fill = "white",
                                      linetype = NULL,
                                      color = "white"),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
  )

NTI.Biome.boxplot

ggsave(filename = "NTI.Biome.boxplot.png", 
       dpi = 300, 
       width = 18.5, height = 5, 
       units = "in")

#######################################################################################################
