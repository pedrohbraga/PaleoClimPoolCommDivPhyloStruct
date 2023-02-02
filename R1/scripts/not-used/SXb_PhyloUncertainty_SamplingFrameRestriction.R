# Recursively read all MPD MNTD files from the phylo sampling procedure

phyloSampling_mpd_mntd_dir <- dir("data/matrices/phyoSampling_mpd_mntd/", full.names = TRUE)

MPD.MNTD.LatLong.AllScales_phyloSampling <- data.frame()

for(full.file.name.i in phyloSampling_mpd_mntd_dir){
  MPD.MNTD.LatLong.AllScales_sampled <- fread(full.file.name.i) %>%
    dplyr::select(-"V1")
  
  MPD.MNTD.LatLong.AllScales_phyloSampling <- MPD.MNTD.LatLong.AllScales_phyloSampling %>%
    bind_rows(MPD.MNTD.LatLong.AllScales_sampled)
}



MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Realm, 
                                                            levels=c('Neotropical',
                                                                     'Nearctic',
                                                                     'Afrotropical',
                                                                     'Palearctic', 
                                                                     'Indomalay', 
                                                                     'Australasian'))

MPD.MNTD.LatLong.AllScales_phyloSampling$SamplingPool <-  factor(MPD.MNTD.LatLong.AllScales_phyloSampling$SamplingPool,  
                                                                 levels = c("Global sampling",
                                                                            "Hemispheric sampling",
                                                                            "Realm sampling",
                                                                            "Plate sampling",
                                                                            "Biome sampling",
                                                                            "Ecoregion sampling"))

MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Biome <- as.factor(MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Biome)

MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Biome_Acronym = factor(MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Biome, 
                                                                   labels = abbreviate(gsub("_", 
                                                                                            " ",
                                                                                            levels(MPD.MNTD.LatLong.AllScales_phyloSampling$ID_Biome))))

MPD.MNTD.LatLong.AllScales_phyloSampling %>%
  filter(phylo.tree.pos %in% c(138, 887, 996))


MPD.MNTD.LatLong.AllScales_phyloSampling_summary <- MPD.MNTD.LatLong.AllScales_phyloSampling %>%
  group_by(ID_SamplingPool) %>%
  summarise(nri.mean.samp = mean(nri, na.rm = T), 
            nti.mean.samp = mean(nti, na.rm = T),
            n.samp = n()) %>%
  right_join(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div, 
             by = "ID_SamplingPool")

MPD.MNTD.LatLong.AllScales_phyloSampling_summary


############################################################################ 
### Representing the phylogenetic community structure of bats across the ###
### gradient of sampling pool (or species pool) restriction              ###
############################################################################

# Representing NRI wrapping across Realms

NRI.Realm.sampling.boxplot <- ggplot(MPD.MNTD.LatLong.AllScales_phyloSampling_summary %>%
                                        filter(!is.na(ID_Realm)),
                                      aes(x = SamplingPool, y = nri.mean.samp, fill = SamplingPool)) +
   geom_boxplot(aes(x = SamplingPool, y = nri.mean.samp)) +
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
     #pretty(MPD.LatLong.Env.AllScales$NRI, n = 7)
     c(min(MPD.MNTD.LatLong.AllScales_phyloSampling_summary$nri.mean.samp, na.rm = T),
       0,
       5,
       10,
       15,
       max(MPD.MNTD.LatLong.AllScales_phyloSampling_summary$nri.mean.samp, na.rm = T)))
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

############################################################################ 
### Representing the phylogenetic community structure of bats across the ###
### gradient of sampling pool (or species pool) restriction              ###
############################################################################

# Representing NTI wrapping across Realms

NTI.Realm.sampling.boxplot <- ggplot(MPD.MNTD.LatLong.AllScales_phyloSampling_summary %>%
                                       filter(!is.na(ID_Realm)),
                                     aes(x = SamplingPool, y = nti.mean.samp, fill = SamplingPool)) +
  geom_boxplot(aes(x = SamplingPool, y = nti)) +
  facet_wrap(~ID_Realm,
             nrow = 1,
             strip.position = "bottom") +
  ggsci::scale_fill_uchicago() +
  # scale_fill_nord("baie_mouton", reverse = TRUE) +
  # scale_fill_okabeito(reverse = TRUE) +
  # scale_fill_brewer(palette = "Dark2") +
  # scale_fill_viridis(discrete = TRUE,
  #                    name="Sampling Pool") +
  scale_y_continuous(
    breaks =
      #   pretty(MPD.MNTD.LatLong.AllScales$nti, n = 5),
      round(
        c(min(MPD.MNTD.LatLong.AllScales_phyloSampling_summary$nti.mean.samp, na.rm = T),
          -2,
          0,
          2,
          4,
          max(MPD.MNTD.LatLong.AllScales_phyloSampling_summary$nti.mean.samp, na.rm = T) + 0.1)
      ),
    limits = c(min(MPD.MNTD.LatLong.AllScales_phyloSampling_summary$nti.mean.samp, na.rm = T) - 0.1,
               max(MPD.MNTD.LatLong.AllScales_phyloSampling_summary$nti.mean.samp, na.rm = T) + 0.1)
    # n.breaks = 4
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
        legend.title = element_text(size = 16),
        legend.title.align = 0.5,
        legend.margin=margin(0, 0, 0, 0),
        legend.box.margin=margin(-5, -5, 2, -5)
  ) +   
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             title.position = "bottom",
                             title="Lower  ⬅  Sampling Pool Restriction  ➡  Higher",)
  )



### Combining plots into a single figure ###

(fig.NRI.NTI.Realm.boxplot.phyloSampling <- ggarrange(
  NRI.Realm.sampling.boxplot +
    theme(axis.title.x = element_blank()), 
  NTI.Realm.sampling.boxplot,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  widths = c(1, 1),
  heights = c(1, 1.1),
  align = "v",
  common.legend = TRUE,
  legend = "bottom"
)
)




NRI.Realm.boxplot.phyloSampling.plot <- function(){
  (NRI.Realm.boxplot.phyloSampling <- ggplot(MPD.MNTD.LatLong.AllScales_phyloSampling %>%
                                               filter(phylo.tree.pos %in% c(138)) %>%
                                               filter(!is.na(ID_Realm)),
                                             aes(x = SamplingPool, y = nri)) +
     geom_boxplot(aes(x = SamplingPool, y = nri, fill = SamplingPool), 
                  alpha = 0.1,
                  outlier.alpha = 0.1, 
                  outlier.colour = NA,
                  lwd = 0.1) +
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
       #pretty(MPD.LatLong.Env.AllScales$NRI, n = 7)
       c(min(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T),
         0,
         5,
         10,
         15,
         max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T))),
       limits = c(min(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T), 
                  max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T))
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
           legend.margin = margin(0, 0, 0, 0),
           legend.box.margin = margin(-5, -5, 2, -5)
     ) +   
     guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                title.position = "bottom",
                                title="Lower  ⬅  Sampling Pool Restriction  ➡  Higher",)
     )
  )
  for(i in unique(MPD.MNTD.LatLong.AllScales_phyloSampling %>% 
                  dplyr::select(phylo.tree.pos))[-1, 1]){ 
    NRI.Realm.boxplot.phyloSampling <- NRI.Realm.boxplot.phyloSampling + 
      geom_boxplot(data = MPD.MNTD.LatLong.AllScales_phyloSampling %>%
                     filter(phylo.tree.pos %in% c(i)) %>%
                     filter(!is.na(ID_Realm)),
                   aes(x = SamplingPool, y = nri, fill = SamplingPool),
                   alpha = 0.1, outlier.colour = NA, lwd = 0.1)
  }
  return(NRI.Realm.boxplot.phyloSampling)
}


NRI.Realm.boxplot.phyloSampling <- NRI.Realm.boxplot.phyloSampling.plot()
NRI.Realm.boxplot.phyloSampling

NTI.Realm.boxplot.phyloSampling.plot <- function() {
  (NTI.Realm.boxplot.phyloSampling <- ggplot(MPD.MNTD.LatLong.AllScales_phyloSampling %>%
                                               filter(phylo.tree.pos %in% c(138)) %>%
                                               filter(!is.na(ID_Realm)),
                                             aes(x = SamplingPool, y = nti)) +
     geom_boxplot(aes(x = SamplingPool, y = nti, fill = SamplingPool), 
                  alpha = 0.1,
                  outlier.alpha = 0.1, 
                  outlier.colour = NA) +
     facet_wrap(~ID_Realm,
                nrow = 1,
                strip.position = "bottom") +
     ggsci::scale_fill_uchicago() +
     # scale_fill_nord("baie_mouton", reverse = TRUE) +
     # scale_fill_okabeito(reverse = TRUE) +
     # scale_fill_brewer(palette = "Dark2") +
     # scale_fill_viridis(discrete = TRUE,
     #                    name="Sampling Pool") +
     scale_y_continuous(
       breaks =
         #   pretty(MPD.MNTD.LatLong.AllScales$nti, n = 5),
         round(
           c(min(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nti, na.rm = T),
             -2,
             0,
             2,
             4,
             max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nti, na.rm = T) + 0.1)
         ),
       limits = c(min(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nti, na.rm = T) - 0.1,
                  max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nti, na.rm = T) + 0.1)
       # n.breaks = 4
     ) +
     geom_hline(yintercept = 0,
                alpha = 0.4) +
     labs(y="NTI", x = "") +
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
           legend.margin = margin(0, 0, 0, 0),
           legend.box.margin = margin(-5, -5, 2, -5)
     ) +   
     guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                title.position = "bottom",
                                title="Lower  ⬅  Sampling Pool Restriction  ➡  Higher",)
     )
  )
  for(i in unique(MPD.MNTD.LatLong.AllScales_phyloSampling %>% 
                  select(phylo.tree.pos))[-1, 1]){ 
    NTI.Realm.boxplot.phyloSampling <- NTI.Realm.boxplot.phyloSampling + 
      geom_boxplot(data = MPD.MNTD.LatLong.AllScales_phyloSampling %>%
                     filter(phylo.tree.pos %in% c(i)) %>%
                     filter(!is.na(ID_Realm)),
                   aes(x = SamplingPool, y = nti, fill = SamplingPool),
                   alpha = 0.1, outlier.colour = NA)
  }
  return(NTI.Realm.boxplot.phyloSampling)
}


NTI.Realm.boxplot.phyloSampling <- NTI.Realm.boxplot.phyloSampling.plot()




### Combining plots into a single figure ###

(fig.NRI.NTI.Realm.boxplot.phyloSampling <- ggarrange(
  NRI.Realm.boxplot.phyloSampling +
    theme(axis.title.x = element_blank()), 
  NTI.Realm.boxplot.phyloSampling,
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

# This is the Figure SX in the manuscript

ggsave(filename = "figures/fig.NRI.NTI.Realm.boxplot.phyloSampling.png", 
       dpi = 300, 
       width = 18.5, height = 10, 
       units = "in")

# Cleaning from memory

rm(NTI.Realm.boxplot,
   NRI.Realm.boxplot,
   NTI.Realm.boxplot.phyloSampling,
   NRI.Realm.boxplot.phyloSampling,
   fig.NRI.NTI.Realm.boxplot.phyloSampling)


