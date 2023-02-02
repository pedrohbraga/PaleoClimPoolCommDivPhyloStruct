##################################################################################
##################################################################################

# To describe how phylogenetic relatedness of bat communities (for both NRI and
# NTI) changed as a function of historical change in temperature, historical
# change in precipitation and in situ net diversification rates, we started by
# estimating their percentiles (i.e., 100-quantiles). We plotted the mean
# phylogenetic relatedness of bat communities (for both NRI and NTI) for each of
# 100-quantiles (percentiles) of the predictor variables of interest (i.e.,
# historical change in temperature, historical change in precipitation and in
# situ net diversification rates; see Figures 3 and 4). These representations
# allowed us to describe how phylogenetic structure varies as a response to the
# predictors of interests.

# To test hypotheses H2 and H3 inferentially, we started by calculating for each
# local community whether its phylogenetic structure (either NRI or NTI) was
# greater than its upper unconditional proportional 90th percentile (we also
# fitted models based on its 75th percentile) across all communities; if
# greater, a value of one was assigned to that community and, if smaller, a
# value of zero was assigned instead. We then applied conditionally unbiased
# bounded influence robust logistic regressions (i.e., robust to reduce the
# influence of potential outliers; see Kunsch et al., 1989) in which the
# response variable was the vector of binary outcomes (ones and zeros)
# representing relatively high or low phylogenetic structure (separately for NRI
# and NTI) and the predictors were z-score standardized values (i.e., to have
# mean zero and variance one) of historical change in temperature, historical
# change in precipitation and in situ net diversification rates. With this, we
# could estimate the relative importance of each predictor within a single
# model. To estimate confidence intervals for each predictor, we used a
# bootstrap approach based on 1000 resamples of 2500 random communities each.
# This approach allowed us to explicitly test how changes in historical climatic
# stability and in situ diversification rates independently increased (or
# decreased) the log-odds (the logistic response) of a community being composed
# of highly phylogenetically related species.

# We repeated these analyses for each sampling frame extent to assess whether
# the influence of historical changes in temperature and precipitation and in
# situ net diversification rates on the phylogenetic structure of bat
# communities was consistent across spatial scales.

##################################################################################

ragg::agg_png("figures/fig.quantile.SamplingPool.NRI.diffTemp.diffPrec.netDiv.png", 
              width = 8*3, 
              height = 7*6, 
              units = "in", 
              res = 250,
              scaling = 1.1)
fig.quantile.SamplingPool.NRI.diffTemp.diffPrec.netDiv
dev.off()

# Bootstrap approach ----

## NRI ----

glm.boot.Global.NRI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nri",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    dplyr::select(c("nri", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Hemispheric.NRI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nri",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    dplyr::select(c("nri", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Realm.NRI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nri",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    dplyr::select(c("nri", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Plate.NRI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nri",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    dplyr::select(c("nri", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Biome.NRI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nri",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    dplyr::select(c("nri", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)


glm.boot.Ecoregion.NRI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nri",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    dplyr::select(c("nri", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

## NTI ----

glm.boot.Global.NTI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nti",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    dplyr::select(c("nti", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Hemispheric.NTI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nti",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    dplyr::select(c("nti", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Realm.NTI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nti",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    dplyr::select(c("nti", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Plate.NTI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nti",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    dplyr::select(c("nti", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

glm.boot.Biome.NTI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nti",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    dplyr::select(c("nti", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)


glm.boot.Ecoregion.NTI.ClimStab.Div.quantiles <- logistic.Phylo.Env(
  community.phylo.variable = "nti",
  X = c("diff.AnnTemp.LGM_cur", 
        "diff.log.AnnPrec.log.LGM_cur",
        "netDiv_CWM_std_tw_rs_1"),
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    dplyr::select(c("nti", 
                    "diff.AnnTemp.LGM_cur", 
                    "diff.log.AnnPrec.log.LGM_cur", 
                    "netDiv_CWM_std_tw_rs_1")) %>%
    drop_na(),
  quantiles = c(0.75, 0.90), 
  n.boot = 1000,
  boot.size = 2500 
)

# Summarise results -----

## Create Figure 5 ----

glm.boot.AllScales.NRI.ClimStab.Div.quantiles.09 <- bind_rows(
  data.frame(glm.boot.Global.NRI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Global sampling"),
  data.frame(glm.boot.Hemispheric.NRI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Hemispheric sampling"),
  data.frame(glm.boot.Realm.NRI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Realm sampling"),
  data.frame(glm.boot.Plate.NRI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Plate sampling"),
  data.frame(glm.boot.Biome.NRI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Biome sampling"),
  data.frame(glm.boot.Ecoregion.NRI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Ecoregion sampling")
) 

glm.boot.AllScales.NTI.ClimStab.Div.quantiles.09 <- bind_rows(
  data.frame(glm.boot.Global.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Global sampling"),
  data.frame(glm.boot.Hemispheric.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Hemispheric sampling"),
  data.frame(glm.boot.Realm.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Realm sampling"),
  data.frame(glm.boot.Plate.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Plate sampling"),
  data.frame(glm.boot.Biome.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Biome sampling"),
  data.frame(glm.boot.Ecoregion.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Ecoregion sampling")
) 


glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09 <- bind_rows(
  data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.09[, 2:5],
             PhyloStructure = "NRI") %>%
    reshape2::melt( 
      id.vars = c("SamplingPool", "PhyloStructure"),
      variable.name = "variable",
      value.name = "estimate"),
  
  data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.09[, 2:5],
             PhyloStructure = "NTI") %>%
    reshape2::melt( 
      id.vars = c("SamplingPool", "PhyloStructure"),
      variable.name = "variable",
      value.name = "estimate")
)


glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09$SamplingPool <- factor(
  glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09$SamplingPool,
  levels = c(
    "Global sampling",
    "Hemispheric sampling",
    "Realm sampling",
    "Plate sampling",
    "Biome sampling",
    "Ecoregion sampling"),
  labels = c(
    "Global~sampling",
    "Hemispheric~sampling",
    "Realm~sampling",
    "Plate~sampling",
    "Biome~sampling",
    "Ecoregion~sampling"
  )
)


glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09$PhyloStructure <- as.factor(
  glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09$PhyloStructure
)

levels(glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09$PhyloStructure)

levels(glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09$PhyloStructure) <- c(
  TeX("$Pr(NRI_{Q90}=1)$"),
  TeX("$Pr(NTI_{Q90}=1)$")
)

# TeX('$Pr(NRI_{Q_{90}}=1)$')

## Create Figure 5 ----

(fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.09 <- ggplot(
  data = glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09,
  aes(x = estimate,
      y = variable)
) +
    geom_tile(data = glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09,
              aes(fill = SamplingPool),
              alpha = 1) +
    geom_boxplot(
      lwd = 0.2,
      outlier.size = 1,
      outlier.stroke = 0.2,
      aes(fill = SamplingPool, 
          colour = SamplingPool
      ))  +
    ggsci::scale_fill_uchicago(alpha = 0.8) +
    ggsci::scale_colour_uchicago() +
    labs(y = "",
         x = "Partial bootstrapped coefficients") +
    scale_y_discrete(
      labels = c("diff.AnnTemp.LGM_cur" = "Historical change in temperature", 
                 "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
                 "netDiv_CWM_std_tw_rs_1" = "_In situ_ diversification rates")
    ) +
    
    # Add white line on top (Inf) of the plot (ie, betweem plot and facet)
    # geom_vline(xintercept = Inf, 
    #            color = "white", 
    #            size = 10) +
    geom_vline(xintercept = 0,
               alpha = 0.4) +
    facet_grid(PhyloStructure ~ SamplingPool,
               # scales = "free",
               switch = "y",
               labeller = label_parsed
    ) +
    # facet_wrap(SamplingPool~PhyloStructure,
    #            nrow = 6,
    #            ncol = 2,
    #            strip.position = c("top",
    #                               "right")) +
    
    theme_boot_quant() +
    theme(
      plot.margin = unit(
        c(1.2, # top
          1, 
          1.2,
          0.5), "lines"),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1.5, "lines"),
      strip.placement = "outside",
      #      strip.background = element_rect(size = 1),
      strip.text.y = element_text(margin = margin(t = 0, 
                                                  r = 2, 
                                                  b = 0, 
                                                  l = 2, "mm")),
      strip.text.x = element_text(margin = margin(t = 0, 
                                                  r = 0, 
                                                  b = 2, 
                                                  l = 0, "mm"),
                                  size = 17),
      axis.text.x = ggtext::element_markdown( 
        size = 15,
        # angle = 90, 
        vjust = 1, 
        hjust = 0.5,
        face = "bold"),
      axis.text.y = ggtext::element_markdown( 
        size = 15,
        hjust = 1),
      axis.title = element_text(size = 16),
      legend.position = "none"
      
    )
)

# Exporting 

as.character(TeX('$Pr(NRI_{Q_{90}}=1)$')) %>%
  names(c("NRI"))

# Export

ragg::agg_png("figures/fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.09.h.png", 
              width = 4*6, 
              height = 4, 
              units = "in", 
              res = 250,
              scaling = 1)

fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.09

dev.off()

# Exporting results for the upper quartile, which is used by the Supplementary
# Information S7

glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075 <- bind_rows(
  data.frame(glm.boot.Global.NRI.ClimStab.Div.quantiles$boot.coefs[[1]],
             SamplingPool = "Global sampling"),
  data.frame(glm.boot.Hemispheric.NRI.ClimStab.Div.quantiles$boot.coefs[[1]],
             SamplingPool = "Hemispheric sampling"),
  data.frame(glm.boot.Realm.NRI.ClimStab.Div.quantiles$boot.coefs[[1]],
             SamplingPool = "Realm sampling"),
  data.frame(glm.boot.Plate.NRI.ClimStab.Div.quantiles$boot.coefs[[1]],
             SamplingPool = "Plate sampling"),
  data.frame(glm.boot.Biome.NRI.ClimStab.Div.quantiles$boot.coefs[[1]],
             SamplingPool = "Biome sampling"),
  data.frame(glm.boot.Ecoregion.NRI.ClimStab.Div.quantiles$boot.coefs[[1]],
             SamplingPool = "Ecoregion sampling")
) 

saveRDS(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075, 
        "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075.RDS")

glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075 <- bind_rows(
  data.frame(glm.boot.Global.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Global sampling"),
  data.frame(glm.boot.Hemispheric.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Hemispheric sampling"),
  data.frame(glm.boot.Realm.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Realm sampling"),
  data.frame(glm.boot.Plate.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Plate sampling"),
  data.frame(glm.boot.Biome.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Biome sampling"),
  data.frame(glm.boot.Ecoregion.NTI.ClimStab.Div.quantiles$boot.coefs[[2]],
             SamplingPool = "Ecoregion sampling")
) 

saveRDS(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075,
        "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075.RDS")


## Create Table 2 -----

glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.09

library(gtsummary)

# bind_rows(
#   data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.09,
#              PhyloStructure = "NRI") %>%
#     reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
#                     variable.name = "variable",
#                     value.name = "estimate"),
#   
#   data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.09,
#              PhyloStructure = "NTI") %>%
#     reshape2::melt( 
#       id.vars = c("SamplingPool", 
#                   "PhyloStructure"),
#       variable.name = "variable",
#       value.name = "estimate")
# ) %>%
#   as.data.frame() %>%
#   select(-PhyloStructure) %>%
#   rename(variables = variable) %>%
#   group_by(variables, SamplingPool) %>%
#   ungroup() %>%
#   gtsummary::tbl_summary(by = SamplingPool,
#                          include = estimate)

# 
# tbl_stack(
#   list(
#     data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.09,
#                PhyloStructure = "NRI") %>%
#       reshape2::melt( 
#         id.vars = c("SamplingPool", "PhyloStructure"),
#         variable.name = "variable",
#         value.name = "estimate") %>%
#       as.matrix() %>%
#       as.data.frame() %>%
#       mutate(estimate, estimate = as.numeric(estimate),
#              variable, Variables = plyr::revalue(
#                variable, 
#                c("X.Intercept." = "Intercept",
#                  "diff.AnnTemp.LGM_cur" = "Historical change in temperature",
#                  "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
#                  "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
#              SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
#                "Global sampling",
#                "Hemispheric sampling",
#                "Realm sampling",
#                "Plate sampling",
#                "Biome sampling",
#                "Ecoregion sampling"
#              ))) %>%
#       mutate(estimate, 
#              estimate = as.numeric(estimate),
#              Variables, 
#              Variables = factor( Variables, 
#                                  levels = c("Intercept",
#                                             "Historical change in temperature",
#                                             "Historical change in precipitation",
#                                             "In situ diversification rates"))) %>%
#       select(-variable, -PhyloStructure) %>%
#       group_by(Variables, SamplingPool) %>%
#       ungroup() %>%
#       gtsummary::tbl_continuous(variable = estimate, 
#                                 by = SamplingPool,
#                                 statistic = NULL) %>%
#       modify_header(all_stat_cols() ~ "**{level}**") %>%
#       bold_labels(), 
#     
#     data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.09,
#                PhyloStructure = "NTI") %>%
#       reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
#                       variable.name = "variable",
#                       value.name = "estimate") %>%
#       as.matrix() %>%
#       as.data.frame() %>%
#       mutate(estimate, estimate = as.numeric(estimate),
#              variable, Variables = plyr::revalue(
#                variable, 
#                c("X.Intercept." = "Intercept",
#                  "diff.AnnTemp.LGM_cur" = "Historical change in temperature",
#                  "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
#                  "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
#              SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
#                "Global sampling",
#                "Hemispheric sampling",
#                "Realm sampling",
#                "Plate sampling",
#                "Biome sampling",
#                "Ecoregion sampling"
#              )
#              )
#       ) %>%
#       mutate(estimate, estimate = as.numeric(estimate),
#              Variables, Variables = factor( 
#                Variables, 
#                levels = c("Intercept",
#                           "Historical change in temperature",
#                           "Historical change in precipitation",
#                           "In situ diversification rates")
#              )
#       ) %>%
#       select(-variable, -PhyloStructure) %>%
#       group_by(Variables, SamplingPool) %>%
#       ungroup() %>%
#       gtsummary::tbl_continuous(variable = estimate, 
#                                 by = SamplingPool,
#                                 statistic = estimate ~ "{mean} ({sd})") %>%
#       modify_header(all_stat_cols() ~ "**{level}**") %>%
#       bold_labels()
#   ), 
#   group_header = c("NRI", "NTI")) # %>%
# as_flex_table() %>%
#   flextable::save_as_docx(path = "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/manuscript/supp_info/table_S3.docx")

# # Mean exp
# 
# mean_exp <- function(data, ...){
#   dplyr::tibble(
#     mean_exp. = exp(mean(data$estimate, 
#                          na.rm = TRUE)),
#     mean = mean(data$estimate, 
#                 na.rm = TRUE),
#     ci.lower = confit(data$estimate, 
#                       level = 0.95)[1],
#     ci.upper = confit(data$estimate, 
#                       level = 0.95)[2]
#     
#   )
# }
# 
# 
# tbl_stack(
#   list(
#     data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.09,
#                PhyloStructure = "NRI") %>%
#       reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
#                       variable.name = "variable",
#                       value.name = "estimate") %>%
#       as.matrix() %>%
#       as.data.frame() %>%
#       mutate(estimate, estimate = as.numeric(estimate),
#              variable, Variables = plyr::revalue(
#                variable, 
#                c("X.Intercept." = "Intercept",
#                  "diff.AnnTemp.LGM_cur" = "Historical change in temperature",
#                  "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
#                  "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")
#              ),
#              SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
#                "Global sampling",
#                "Hemispheric sampling",
#                "Realm sampling",
#                "Plate sampling",
#                "Biome sampling",
#                "Ecoregion sampling"
#              ))) %>%
#       mutate(estimate, 
#              estimate = as.numeric(estimate),
#              Variables, 
#              Variables = factor( Variables, 
#                                  levels = c("Intercept",
#                                             "Historical change in temperature",
#                                             "Historical change in precipitation",
#                                             "In situ diversification rates"))) %>%
#       select(-variable, -PhyloStructure) %>%
#       group_by(Variables, SamplingPool) %>%
#       ungroup() %>%
#       group_nest(Variables) %>%
#       mutate(
#         tbl = map2(
#           Variables, data, 
#           ~tbl_custom_summary(.y, 
#                               by = SamplingPool, 
#                               type = list(estimate ~ 'continuous'),
#                               label = list(estimate = paste("", .x)), 
#                               missing = "no",
#                               statistic = "estimate" ~ "{mean} ({mean_exp.})",
#                               stat_fns = everything() ~ mean_exp
#           )
#         )
#       ) %>%
#       pull(tbl) %>%
#       tbl_stack() %>%
#       modify_footnote(all_stat_cols() ~ "Mean (Odds Ratio)") %>%
#       bold_labels(),
#     
#     #
#     #
#     
#     data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.09,
#                PhyloStructure = "NTI") %>%
#       reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
#                       variable.name = "variable",
#                       value.name = "estimate") %>%
#       as.matrix() %>%
#       as.data.frame() %>%
#       mutate(estimate, estimate = as.numeric(estimate),
#              variable, Variables = plyr::revalue(
#                variable, 
#                c("X.Intercept." = "Intercept",
#                  "diff.AnnTemp.LGM_cur" = "Historical change in temperature",
#                  "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
#                  "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
#              SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
#                "Global sampling",
#                "Hemispheric sampling",
#                "Realm sampling",
#                "Plate sampling",
#                "Biome sampling",
#                "Ecoregion sampling"
#              ))) %>%
#       mutate(estimate, 
#              estimate = as.numeric(estimate),
#              Variables, 
#              Variables = factor( Variables, 
#                                  levels = c("Intercept",
#                                             "Historical change in temperature",
#                                             "Historical change in precipitation",
#                                             "In situ diversification rates"))) %>%
#       select(-variable, -PhyloStructure) %>%
#       group_by(Variables, SamplingPool) %>%
#       ungroup() %>%
#       group_nest(Variables) %>%
#       mutate(
#         tbl = map2(
#           Variables, data, 
#           ~tbl_custom_summary(.y, 
#                               by = SamplingPool, 
#                               type = list(estimate ~ 'continuous'),
#                               label = list(estimate = paste("", .x)), 
#                               missing = "no",
#                               statistic = "estimate" ~ "{mean} ({mean_exp.})",
#                               stat_fns = everything() ~ mean_exp
#           )
#         )
#       ) %>%
#       pull(tbl) %>%
#       tbl_stack() %>%
#       modify_footnote(all_stat_cols() ~ "Mean (Odds Ratio)") %>%
#       bold_labels() 
#   ),
#   group_header = c("NRI", 
#                    "NTI")
# )  %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path = "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/manuscript/supp_info/supp_info_S2_table_S2.3.docx")
# 


mean_exp <- function(data, ...){
  dplyr::tibble(
    mean_exp. = exp(mean(data$estimate, 
                         na.rm = TRUE)),
    mean = mean(data$estimate, na.rm = TRUE),
    ci.lower = confint(data$estimate,
                       level = 0.95,
                       method = "stderr")[,1],
    ci.upper = confint(data$estimate, 
                       level = 0.95,
                       method = "stderr")[,2]
    
  )
}

# Table 2

tbl_stack(
  list(
    data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.09,
               PhyloStructure = "NRI") %>%
      reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
                      variable.name = "variable",
                      value.name = "estimate") %>%
      as.matrix() %>%
      as.data.frame() %>%
      mutate(estimate, estimate = as.numeric(estimate),
             variable, Variables = plyr::revalue(variable, 
                                                 c("X.Intercept." = "Intercept",
                                                   "diff.AnnTemp.LGM_cur" = "Historical change in temperature", 
                                                   "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
                                                   "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates"
                                                 )
             ),
             SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
               "Global sampling",
               "Hemispheric sampling",
               "Realm sampling",
               "Plate sampling",
               "Biome sampling",
               "Ecoregion sampling"
             ))) %>%
      mutate(estimate, 
             estimate = as.numeric(estimate),
             Variables, 
             Variables = factor(Variables, 
                                levels = c("Intercept",
                                           "Historical change in temperature",
                                           "Historical change in precipitation",
                                           "In situ diversification rates"
                                ))) %>%
      select(-variable, -PhyloStructure) %>%
      group_by(Variables, SamplingPool) %>%
      ungroup() %>%
      group_nest(Variables) %>%
      mutate(
        tbl = map2(
          Variables, data, 
          ~tbl_custom_summary(.y, 
                              by = SamplingPool, 
                              type = list(estimate ~ 'continuous'),
                              label = list(estimate = paste("", .x)), 
                              missing = "no",
                              statistic = "estimate" ~ "{mean} ({ci.lower}; {ci.upper})",
                              stat_fns = everything() ~ mean_exp
          )
        )
      ) %>%
      pull(tbl) %>%
      tbl_stack() %>%
      modify_footnote(all_stat_cols() ~ "Mean (Lower CI; Upper CI)") %>%
      bold_labels(),
    
    #
    data.frame(glm.boot.AllScales.mean.NTI.rarefaction.relative.ClimStab.Div.quantiles.09,
               PhyloStructure = "NTI") %>%
      reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
                      variable.name = "variable",
                      value.name = "estimate") %>%
      as.matrix() %>%
      as.data.frame() %>%
      mutate(estimate, estimate = as.numeric(estimate),
             variable, Variables = plyr::revalue(
               variable, 
               c("X.Intercept." = "Intercept",
                 "diff.AnnTemp.LGM_cur" = "Historical change in temperature", 
                 "diff.log.AnnPrec.log.LGM_cur" = "Historical change in precipitation",
                 "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates"
                 
               )),
             SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
               "Global sampling",
               "Hemispheric sampling",
               "Realm sampling",
               "Plate sampling",
               "Biome sampling",
               "Ecoregion sampling"
             ))) %>%
      mutate(estimate, 
             estimate = as.numeric(estimate),
             Variables, 
             Variables = factor(
               Variables, 
               levels = c("Intercept",
                          "Historical change in temperature",
                          "Historical change in precipitation",
                          "In situ diversification rates"
               ))) %>%
      select(-variable, -PhyloStructure) %>%
      group_by(Variables, SamplingPool) %>%
      ungroup() %>%
      group_nest(Variables) %>%
      mutate(
        tbl = map2(
          Variables, data, 
          ~tbl_custom_summary(
            .y, 
            by = SamplingPool, 
            type = list(estimate ~ 'continuous'),
            label = list(estimate = paste("", .x)), 
            missing = "no",
            statistic = "estimate" ~ "{mean} ({ci.lower}; {ci.upper})",
            stat_fns = everything() ~ mean_exp
          )
        )
      ) %>%
      pull(tbl) %>%
      tbl_stack() %>%
      modify_footnote(all_stat_cols() ~ "Mean (Lower CI; Upper CI)") %>%
      modify_caption(caption = "") %>%
      bold_labels() 
  ),
  group_header = c("NRI", 
                   "NTI")
) %>% 
  bold_labels() %>%
  # as_hux_table() %>%      
  # set_width(1.3) %>% 
  # set_font_size(8) %>% 
  # huxtable::theme_article(header_rows = TRUE, header_cols = TRUE) %>%
  # set_caption(NA) %>%
  # set_caption_pos("topleft") %>%
  #  
  as_flex_table() %>%
  flextable::fontsize(size = 7, part = "all") %>%
  # flextable::autofit() %>%
  # flextable::set_table_properties(layout = "autofit") %>%
  # flextable::fit_to_width(7.5) %>%
  flextable::save_as_docx(path = "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/manuscript/Table-2.docx")


