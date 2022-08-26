##################################################################################
## Combine p-values from the null models estimating standardized effect sizes of #
## community phylogenetic relatedness.
##################################################################################

# Create new p-value, where the p-value for negative NRI and NTI values becomes
# the complement of the original p-value, i.e. we take its reciprocal as in 1 - p

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  mutate(mntd.obs.p.complement.neg = ifelse(nti < 0, 1 - mntd.obs.p, mntd.obs.p),
         mpd.obs.p.complement.neg = ifelse(nri < 0, 1 - mpd.obs.p, mpd.obs.p)) 

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

# Using poolr::stouffer()$p for NTI ----

(stouffer.comb.p.val.nti <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
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


# Option B: Test whether NRI or NTI significantly differs from zero ----

lm(nti ~ 1,
   data = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
     dplyr::filter(SamplingPool == "Ecoregion sampling") %>%
     dplyr::filter(ID_Realm == "Afrotropical") %>%
     drop_na(nti)) %>%
  summary()

# Garbage ----

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ID_Realm == "Indomalay") %>%
  drop_na(ses.mpd.z.query) %>%
  pull(ses.mpd.z.query) %>%
  hist()



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

