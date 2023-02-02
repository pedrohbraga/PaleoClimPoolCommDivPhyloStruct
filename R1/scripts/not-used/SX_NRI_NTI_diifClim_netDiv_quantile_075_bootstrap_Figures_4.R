
##### Figure creation #####
# q = 0.9 ######

glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075 <- bind_rows(data.frame(glm.boot.Global.NRI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Global sampling"),
                                                              data.frame(glm.boot.Hemispheric.NRI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Hemispheric sampling"),
                                                              data.frame(glm.boot.Realm.NRI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Realm sampling"),
                                                              data.frame(glm.boot.Plate.NRI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Plate sampling"),
                                                              data.frame(glm.boot.Biome.NRI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Biome sampling"),
                                                              data.frame(glm.boot.Ecoregion.NRI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Ecoregion sampling")
) 

glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075 <- bind_rows(data.frame(glm.boot.Global.NTI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Global sampling"),
                                                              data.frame(glm.boot.Hemispheric.NTI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Hemispheric sampling"),
                                                              data.frame(glm.boot.Realm.NTI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Realm sampling"),
                                                              data.frame(glm.boot.Plate.NTI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Plate sampling"),
                                                              data.frame(glm.boot.Biome.NTI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Biome sampling"),
                                                              data.frame(glm.boot.Ecoregion.NTI.ClimStab.Div.quantiles$boot.coefs[[5]],
                                                                         SamplingPool = "Ecoregion sampling")
) 


glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075 <- bind_rows(
  data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075[, 2:5],
             PhyloStructure = "NRI") %>%
    reshape2::melt( 
      id.vars = c("SamplingPool", "PhyloStructure"),
      variable.name = "variable",
      value.name = "estimate"),
  
  data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075[, 2:5],
             PhyloStructure = "NTI") %>%
    reshape2::melt( 
      id.vars = c("SamplingPool", "PhyloStructure"),
      variable.name = "variable",
      value.name = "estimate")
)


glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075$SamplingPool <- factor(glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075$SamplingPool,
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


glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075$PhyloStructure <- as.factor(glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075$PhyloStructure)

levels(glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075$PhyloStructure)

levels(glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075$PhyloStructure) <- c(TeX("$Pr(NRI_{Q90}=1)$"), TeX("$Pr(NTI_{Q90}=1)$"))

saveRDS(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075, 
        "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075.RDS")
saveRDS(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075, 
        "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075.RDS")

# TeX('$Pr(NRI_{Q_{90}}=1)$')

(fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.075 <- ggplot(data = glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075,
                                                                    aes(x = estimate,
                                                                        y = variable)
) +
    geom_tile(data = glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075,
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



as.character(TeX('$Pr(NRI_{Q_{90}}=1)$')) %>%
  names(c("NRI"))


ggsave(
  file = "figures/fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.09.h.svg",
  fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.075,
  width = 4*6, 
  height = 4, 
  units = c("in"),
  dpi = 250,
  limitsize = FALSE
)

ragg::agg_png("figures/fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.075.h.png", 
              width = 4*6, 
              height = 4, 
              units = "in", 
              res = 250,
              scaling = 1)
fig.wrap.boot.AllScales.NRI.NTI.ClimStab.Div.quantile.075
dev.off()

## Table generation ####

glm.boot.AllScales.NRI.NTI.ClimStab.Div.quantiles.075


library(gtsummary)

bind_rows(
  data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075,
             PhyloStructure = "NRI") %>%
    reshape2::melt( id.vars = c("SamplingPool", "PhyloStructure"),
                    variable.name = "variable",
                    value.name = "estimate"),
  
  data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075,
             PhyloStructure = "NTI") %>%
    reshape2::melt( 
      id.vars = c("SamplingPool", 
                  "PhyloStructure"),
      variable.name = "variable",
      value.name = "estimate")
) %>%
  as.data.frame() %>%
  select(-PhyloStructure) %>%
  rename(variables = variable) %>%
  group_by(variables, SamplingPool) %>%
  ungroup() %>%
  gtsummary::tbl_summary(by = SamplingPool,
                         include = estimate)


tbl_stack(
  list(
    data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075,
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
                                                   "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
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
             Variables = factor( Variables, 
                                 levels = c("Intercept",
                                            "Historical change in temperature",
                                            "Historical change in precipitation",
                                            "In situ diversification rates"))) %>%
      select(-variable, -PhyloStructure) %>%
      group_by(Variables, SamplingPool) %>%
      ungroup() %>%
      gtsummary::tbl_continuous(variable = estimate, 
                                by = SamplingPool,
                                statistic = NULL) %>%
      modify_header(all_stat_cols() ~ "**{level}**") %>%
      bold_labels(), 
    
    data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075,
               PhyloStructure = "NTI") %>%
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
                                                   "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
             SamplingPool, SamplingPool = factor(SamplingPool, levels = c(
               "Global sampling",
               "Hemispheric sampling",
               "Realm sampling",
               "Plate sampling",
               "Biome sampling",
               "Ecoregion sampling"
             )
             )
      ) %>%
      mutate(estimate, estimate = as.numeric(estimate),
             Variables, Variables = factor( Variables, 
                                            levels = c("Intercept",
                                                       "Historical change in temperature",
                                                       "Historical change in precipitation",
                                                       "In situ diversification rates"))) %>%
      select(-variable, -PhyloStructure) %>%
      group_by(Variables, SamplingPool) %>%
      ungroup() %>%
      gtsummary::tbl_continuous(variable = estimate, 
                                by = SamplingPool,
                                statistic = estimate ~ "{mean} ({sd})") %>%
      modify_header(all_stat_cols() ~ "**{level}**") %>%
      bold_labels()
  ), 
  group_header = c("NRI", "NTI")) # %>%
as_flex_table() %>%
  flextable::save_as_docx(path = "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/manuscript/supp_info/table_S3.docx")





mean_exp <- function(data, ...){
  dplyr::tibble(
    mean_exp. = exp(mean(data$estimate, na.rm = TRUE)),
    mean = mean(data$estimate, na.rm = TRUE),
    ci.lower = confit(data$estimate, 
                      level = 0.95)[1],
    ci.upper = confit(data$estimate, 
                      level = 0.95)[2]
    
  )
}


tbl_stack(
  list(
    data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075,
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
                                                   "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
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
             Variables = factor( Variables, 
                                 levels = c("Intercept",
                                            "Historical change in temperature",
                                            "Historical change in precipitation",
                                            "In situ diversification rates"))) %>%
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
                              statistic = "estimate" ~ "{mean} ({mean_exp.})",
                              stat_fns = everything() ~ mean_exp
          )
        )
      ) %>%
      pull(tbl) %>%
      tbl_stack() %>%
      modify_footnote(all_stat_cols() ~ "Mean (Odds Ratio)") %>%
      bold_labels(),
    
    #
    #
    
    data.frame(glm.boot.AllScales.NTI.ClimStab.Div.quantiles.075,
               PhyloStructure = "NTI") %>%
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
                                                   "netDiv_CWM_std_tw_rs_1" = "In situ diversification rates")),
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
             Variables = factor( Variables, 
                                 levels = c("Intercept",
                                            "Historical change in temperature",
                                            "Historical change in precipitation",
                                            "In situ diversification rates"))) %>%
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
                              statistic = "estimate" ~ "{mean} ({mean_exp.})",
                              stat_fns = everything() ~ mean_exp
          )
        )
      ) %>%
      pull(tbl) %>%
      tbl_stack() %>%
      modify_footnote(all_stat_cols() ~ "Mean (Odds Ratio)") %>%
      bold_labels() 
  ),
  group_header = c("NRI", 
                   "NTI")
)  %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/manuscript/supp_info/supp_info_S1_table_S2.3.docx")



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

tbl_stack(
  list(
    data.frame(glm.boot.AllScales.NRI.ClimStab.Div.quantiles.075,
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
    #
    
    data.frame(glm.boot.AllScales.mean.NTI.rarefaction.relative.ClimStab.Div.quantiles.075,
               PhyloStructure = "NTI") %>%
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
             Variables = factor( Variables, 
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
  flextable::save_as_docx(path = "~/chapter-PaleoClimPoolCommDivPhyloStruct/R1/manuscript/supp_info/supp_info_S7_table_S2.3.docx")


