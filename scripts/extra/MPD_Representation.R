###########################################################
#### Statistical analyses and graphical representation ####
###########################################################

#########################################################

## How does NRI across different biogeographical realms?
# NRI ~ Realm

# Linear Model
## Global sampling

MPD.Global <- MPD.LatLong.Env.AllScales %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(!is.na(ID_Realm))

# Assumption testing
car::leveneTest(NRI ~ ID_Realm, data = MPD.Global, na.rm = TRUE)
# Assumptions violated

Bats.NRI.Realm <- aov(NRI ~ ID_Realm, data = MPD.Global)
summary(Bats.NRI.Realm)

TukeyHSD(x = Bats.NRI.Realm, 
         'ID_Realm', 
         conf.level=0.95)

# Non-parametric alternative

kruskal.test(NRI ~ ID_Realm, data = MPD.Global)

pairwise.wilcox.test(MPD.Global$NRI, MPD.Global$ID_Realm,
                     p.adjust.method = "BH",
                     paired = TRUE)

x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
y <- c(3.8, 2.7, 4.0, 2.4) # with obstructive airway disease
z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
dunn.test(x=list(x,y,z))
x <- c(x, y, z)
g <- factor(rep(1:3, c(5, 4, 5)),
            labels = c("Normal",
                       "COPD",
                       "Asbestosis"))
dunn.test(x, g)

dunn.test(MPD.Global$NRI, MPD.Global$ID_Realm., 
          kw = TRUE, 
          method = "by",
          rmc = TRUE,
          altp = TRUE,
          list = TRUE)

dunn.test(MNTD.Global$NTI, MNTD.Global$ID_Realm., 
          kw = TRUE, 
          method = "by",
          rmc = TRUE,
          altp = TRUE,
          list = TRUE)




# Linear Mixed Model

Bats.NRI.Realm.LMM <- lmer(NRI ~ ID_Realm + (1|SamplingPool), 
                           data = MPD.LatLong.Env.AllScales)
plot(Bats.NRI.Realm.LMM)

summary(Bats.NRI.Realm.LMM)


#################################################################
#### 

NRI.SamplingPool.ID_Realm <- list(filter(MPD.LatLong.Env.AllScales, 
                                         SamplingPool == "Global sampling")[, c(1, 6, 36)],
                                  filter(MPD.LatLong.Env.AllScales, 
                                         SamplingPool == "Hemispheric sampling")[, c(1, 6, 36)],
                                  filter(MPD.LatLong.Env.AllScales, 
                                         SamplingPool == "Realm sampling")[, c(1, 6, 36)],
                                  filter(MPD.LatLong.Env.AllScales, 
                                         SamplingPool == "Biome sampling")[, c(1, 6, 36)]) %>%
  Reduce(function(dtf1, dtf2) left_join(dtf1,dtf2,by="ID"), .)

colnames(NRI.SamplingPool.ID_Realm) <- c("ID",
                                         "ID_Realm", "NRI.Global",
                                         "ID_Realm.", "NRI.Hemispheric",
                                         "ID_Realm..", "NRI.Realm",
                                         "ID_Realm...", "NRI.Biome")

diff.NRI.SamplingPool <- data.frame(ID = NRI.SamplingPool.ID_Realm$ID,
           ID_Realm = NRI.SamplingPool.ID_Realm$ID_Realm,
           diff.NRI.Global.Hemispheric = NRI.SamplingPool.ID_Realm$NRI.Global - NRI.SamplingPool.ID_Realm$NRI.Hemispheric,
           diff.NRI.Hemispheric.Realm = NRI.SamplingPool.ID_Realm$NRI.Hemispheric - NRI.SamplingPool.ID_Realm$NRI.Realm,
           diff.NRI.Realm.Biome = NRI.SamplingPool.ID_Realm$NRI.Realm - NRI.SamplingPool.ID_Realm$NRI.Biome)

diff.NRI.SamplingPool <- diff.NRI.SamplingPool %>%
  gather(diff.SamplingPool, diff.NRI, -ID, -ID_Realm)


#### Representation of the delta NRI

diff.NRI.SamplingPool.boxplot <- ggplot(filter(diff.NRI.SamplingPool,
                                   is.na(ID_Realm) == FALSE),
                            aes(x = diff.SamplingPool, y = diff.NRI)) +
  geom_boxplot(aes(fill = diff.SamplingPool)) +
  facet_wrap(~ID_Realm,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(breaks = round(
    #pretty(MPD.LatLong.Env.AllScales$NRI, n = 7)
    c(min(diff.NRI.SamplingPool$diff.NRI, na.rm = T),
      0,
      5,
      10,
      15,
      max(diff.NRI.SamplingPool$diff.NRI, na.rm = T)))
  ) +
  geom_hline(yintercept=0,
             alpha = 0.4) +
  labs(x="",
       y = expression(Delta*"NRI")) +
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

diff.NRI.SamplingPool.boxplot

ggsave(filename = "diff.NRI.SamplingPool.boxplot.png", 
       dpi = 300, 
       width = 17.5, height = 5, 
       units = "in")

# Effects of sampling pool scaling on NRI within biogeographic realms. Within
# each panel, the left box plot shows the change in NRI between global and
# hemispheric sampling pools (i.e., NRIglobal – NRIhemispheric). The right box
# plot shows the change in NRI between hemispheric and realm sampling pools
# (NRIhemispheric – NRIrealm). Values > 0 indicate stronger phylogenetic
# clustering with the geographically more extensive sampling pool.

########################################################

MPD.LatLong.Env.AllScales$ID_Biome_Acronym = factor(MPD.LatLong.Env.AllScales$ID_Biome, 
                                                    labels = abbreviate(gsub("_", 
                                                                             " ",
                                                                             levels(MPD.LatLong.Env.AllScales$ID_Biome))))

NRI.Realm.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales,
                                   is.na(ID_Realm) == FALSE),
                            aes(x = SamplingPool, y = NRI)) +
  geom_boxplot(aes(fill = SamplingPool)) +
  facet_wrap(~ID_Realm,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(breaks = round(
    #pretty(MPD.LatLong.Env.AllScales$NRI, n = 7)
    c(min(MPD.LatLong.Env.AllScales$NRI, na.rm = T),
      0,
      5,
      10,
      15,
      max(MPD.LatLong.Env.AllScales$NRI, na.rm = T)))
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

NRI.Realm.boxplot

ggsave(filename = "NRI.Realm.boxplot.png", 
       dpi = 300, 
       width = 18.5, height = 5, 
       units = "in")

##
NRI.Biome.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales,
                                   is.na(ID_Biome) == FALSE),
                            aes(x = SamplingPool, y = NRI)) +
  geom_boxplot(aes(fill = SamplingPool)) +
  facet_wrap(~ID_Biome_Acronym,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(breaks = round(
    #pretty(MPD.LatLong.Env.AllScales$NRI, n = 7)
    c(min(MPD.LatLong.Env.AllScales$NRI, na.rm = T),
      0,
      5,
      10,
      15,
      max(MPD.LatLong.Env.AllScales$NRI, na.rm = T)))
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
        # axis.text.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14), #
        axis.title = element_text(size = 16, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
  )

NRI.Biome.boxplot

ggsave(filename = "NRI.Biome.boxplot.png", 
       dpi = 300, 
       width = 18.5, height = 5, 
       units = "in")
