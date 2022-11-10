##################################################################################
## Combine p-values from the null models estimating standardized effect sizes of #
## community phylogenetic relatedness.
##################################################################################

# Create new p-value, where the p-value for negative NRI and NTI values becomes
# the complement of the original p-value, i.e. we take its reciprocal as in 1 - p

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  mutate(mntd.obs.p.complement.neg = ifelse(nti.rarefac.mean < 0, 1 - mntd.obs.p, mntd.obs.p),
         mpd.obs.p.complement.neg = ifelse(nri.rarefac.mean < 0, 1 - mpd.obs.p, mpd.obs.p)) 

# 
# MPD.MNTD.LatLong.AllScales <-  MPD.MNTD.LatLong.AllScales %>%
#   mutate(mntd.obs.p.complement.neg = ifelse(nti < 0, 1 - mntd.obs.p, mntd.obs.p),
#          mpd.obs.p.complement.neg = ifelse(nri < 0, 1 - mpd.obs.p, mpd.obs.p)) 


# Use Fisher's method for combining p-values

# Using my own function

pValueCombFisher <- function(p.vector){
  
  # not used
  # find p-value for the z-score under a two-tailed test
  # p.vector <- 2*pnorm(-abs(z.vector), lower.tail = TRUE)
  # end of not used
  
  k.length.p.vector <- length(p.vector)
  degrees.of.freedom <- k.length.p.vector*2
  y <- -2*sum(log(p.vector))
  
  fisher.comb.p.val <- pchisq(y, df = degrees.of.freedom, lower.tail = FALSE)
  
  # should be the same as
  # fisher.comb.p.val <- 1 - pchisq(y, df = degrees.of.freedom); 
  
  
  fisher.comb.p.val <- as.numeric(fisher.comb.p.val)
  
  return(fisher.comb.p.val)
}

zValueCombStouffer <- function(z.vector){
  
  k.length.z.vector <- length(z.vector)
  
  stouffer.comb.z.val <- sum(z.vector)/sqrt(k.length.z.vector)
  
  return(stouffer.comb.z.val)
}



# Testing

# MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
#   filter(SamplingPool == "Ecoregion sampling") %>%
#   filter(ID_Realm == "Afrotropical") %>%
#   drop_na(mntd.obs.p.complement.neg) %>%
#   pull(mntd.obs.p.complement.neg) %>%
#   pValueCombFisher()


# pValueCombFisher() for NRI ----

(fisher.comb.p.val.nri <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
   group_by(ID_Realm, SamplingPool) %>%
   drop_na(mntd.obs.p.complement.neg,
           mpd.obs.p.complement.neg,
           ID_Realm) %>%
   group_modify(~ {
     pValueCombFisher(.x$mpd.obs.p.complement.neg) %>%
       tibble::enframe(name = "mpd.name",
                       value = "mpd.fisher.comb.p.val")
   }) %>%
   as.data.frame()
)

# pValueCombFisher() for NTI ----

(fisher.comb.p.val.nti <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    group_by(ID_Realm, SamplingPool) %>%
    drop_na(mntd.obs.p.complement.neg,
            ID_Realm) %>%
    group_modify(~ {
      pValueCombFisher(.x$mntd.obs.p.complement.neg) %>%
        tibble::enframe(name = "mntd.name",
                        value = "mntd.fisher.comb.p.val")
    }) %>%
    as.data.frame()
)
# poolr::fisher()$p for NRI -----

(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  group_by(ID_Realm, SamplingPool) %>%
  drop_na(mntd.obs.p.complement.neg,
          mpd.obs.p.complement.neg,
          ID_Realm) %>%
  group_modify(~ {
    poolr::fisher(.x$mpd.obs.p.complement.neg)$p %>%
      tibble::enframe(name = "mpd.name",
                      value = "mpd.poolr.fisher.comb.p.val")
  }) %>%
  as.data.frame()
 )

# poolr::fisher()$p for NTI ----

(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  group_by(ID_Realm, SamplingPool) %>%
  drop_na(mntd.obs.p.complement.neg,
          mpd.obs.p.complement.neg,
          ID_Realm) %>%
  group_modify(~ {
    poolr::fisher(.x$mntd.obs.p.complement.neg)$p %>%
      tibble::enframe(name = "mntd.name",
                      value = "mntd.poolr.fisher.comb.p.val")
  }) %>%
  as.data.frame()
 )

# Using poolr::stouffer()$p for NRI ----

(stouffer.comb.p.val.nri <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  group_by(ID_Realm, SamplingPool) %>%
  drop_na(mntd.obs.p.complement.neg,
          mpd.obs.p.complement.neg,
          ID_Realm) %>%
  group_modify(~ {
    poolr::stouffer(.x$mpd.obs.p.complement.neg)$p %>%
      tibble::enframe(name = "mpd.name",
                      value = "mpd.poolr.stouffer.comb.p.val")
  }) %>%
  as.data.frame()
)


(stouffer.comb.p.val.nri <- MPD.MNTD.LatLong.AllScales %>%
    group_by(ID_Realm, SamplingPool) %>%
    drop_na(mntd.obs.p.complement.neg,
            mpd.obs.p.complement.neg,
            ID_Realm) %>%
    group_modify(~ {
      poolr::stouffer(.x$mpd.obs.p.complement.neg)$p %>%
        tibble::enframe(name = "mpd.name",
                        value = "mpd.poolr.stouffer.comb.p.val")
    }) %>%
    as.data.frame()
)

# Using poolr::stouffer()$p for NTI ----

(stouffer.comb.p.val.nti <- MPD.MNTD.LatLong.AllScales %>%
  group_by(ID_Realm, SamplingPool) %>%
  drop_na(mntd.obs.p.complement.neg,
          mpd.obs.p.complement.neg,
          ID_Realm) %>%
  group_modify(~ {
    poolr::stouffer(.x$mntd.obs.p.complement.neg)$p %>%
      tibble::enframe(name = "mntd.name",
                      value = "mntd.poolr.stouffer.comb.p.val")
  }) %>%
  as.data.frame()
)

# Tabling results ----

left_join(fisher.comb.p.val.nri,
          fisher.comb.p.val.nti) %>%
  dplyr::select(-mpd.name, -mntd.name) %>%
  kable("markdown")

left_join(stouffer.comb.p.val.nri,
          stouffer.comb.p.val.nti) %>%
  dplyr::select(-mpd.name, -mntd.name) %>%
  kable("markdown")

left_join(fisher.comb.p.val.nri,
          fisher.comb.p.val.nti) %>%
  dplyr::select(-mpd.name, -mntd.name) 



# pValueCombFisher() for NRI ----

(stouffer.comb.z.val.nri <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
   group_by(ID_Realm, SamplingPool) %>%
   drop_na(nri,
           nti,
           ID_Realm) %>%
   group_modify(~ {
     zValueCombStouffer(.x$nri) %>%
       tibble::enframe(name = "z.nri.name",
                       value = "nri.stoufer.comb.z.val")
   }) %>%
   as.data.frame()
)


####################
##### Figure 2 #####
####################

# Representing NRI wrapping across Realms

(NRI.Realm.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div,
                                    is.na(ID_Realm) == FALSE),
                             aes(x = SamplingPool, y = nri)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    facet_wrap(~ID_Realm,
               nrow = 1,
               strip.position = "bottom") +
    ggsci::scale_fill_uchicago() +
    scale_y_continuous(breaks = round(
      c(min(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T),
        0,
        5,
        10,
        15,
        max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T)))
    ) +
    geom_text(data = stouffer.comb.p.val.nri, 
              aes(y = max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nri, na.rm = T)*1.05, 
                  label = stouffer.comb.p.val.nri %>%
                    mutate(annot = ifelse(mpd.poolr.stouffer.comb.p.val < 0.05, "*", NA)) %>%
                    pull(annot)
                  ), 
              position = position_dodge(width = .75),
              size = 7) +
    geom_hline(yintercept = 0,
               alpha = 0.4) +
    labs(x = "",
         y = "NRI") +
    theme_minimal(base_size = 20) +
    theme(strip.background = element_rect(fill = "white",
                                          linetype = NULL,
                                          color = "white"),
          strip.text = element_text(color = "black",
                                    face = "bold",
                                    size = 15),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 16, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 12,
                                     margin = margin(r = 0.8, unit = 'cm')),
          legend.title = element_text(size = 13),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.key.width = unit(1.4, "line"), 
          legend.title.align = 0.5,
          legend.margin = margin(t = 0.5, r = 0, b = 1, l = 0),
          legend.box.margin = margin(-2, -2, -2, -2)
    ) +
    guides(fill = guide_legend(nrow = 1, 
                               byrow = TRUE,
                               title.position = "bottom",
                               title = "Lower \U2190  Geographical Extent Restriction  \U2192  Higher")
           
    )
  
           
)


# If necessary, save figure
# ggsave(filename = "NRI.Realm.boxplot.png", 
#        dpi = 300, 
#        width = 18.5, height = 5, 
#        units = "in")

# Representing NTI wrapping across Realms

(NTI.Realm.boxplot <- ggplot(filter(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div,
                                    is.na(ID_Realm) == FALSE),
                             aes(x = SamplingPool, y = nti)) +
    geom_boxplot(aes(fill = SamplingPool)) +
    facet_wrap(~ID_Realm,
               nrow = 1,
               strip.position = "bottom") +
    ggsci::scale_fill_uchicago() +
    geom_text(data = stouffer.comb.p.val.nti, 
              aes(y = max(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$nti, na.rm = T)*1.05, 
                  label = stouffer.comb.p.val.nri %>%
                    mutate(annot = ifelse(mntd.poolr.stouffer.comb.p.val < 0.05, "*", NA)) %>%
                    pull(annot)), 
              position = position_dodge(width = .75),
              size = 7) +
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
    labs(x = "",
         y = "NTI") +
    
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
          legend.text = element_text(size = 12,
                                     margin = margin(r = 0.8, unit = 'cm')),
          legend.title = element_text(size = 13),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.key.width = unit(1.4, "line"), 
          legend.title.align = 0.5,
          legend.margin = margin(t = 0.5, r = 0, b = 1, l = 0),
          legend.box.margin = margin(-2, -2, -2, -2)
    ) +
    guides(fill = guide_legend(nrow = 1, 
                               byrow = TRUE,
                               title.position = "bottom",
                               title = "Lower \U2190  Geographical Extent Restriction  \U2192  Higher")
           
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


ragg::agg_png("figures/fig.NRI.NTI.Realm.boxplot.png", 
        width = 1548*8, height = 900*8, res = 650)

fig.NRI.NTI.Realm.boxplot

invisible(dev.off())


ggsave(filename = "figures/fig.NRI.NTI.Realm.boxplot.png", 
       dpi = 600, 
       width = 22, height = 10, 
       units = "in")

