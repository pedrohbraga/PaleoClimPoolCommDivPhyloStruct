# quantile.XY for NTI

quantile.nti.XY.all <- data.frame()

for(SamplingPool_i in levels(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$SamplingPool)){
  for(X.var_i in c("diff.AnnTemp.LGM_cur", "diff.log.AnnPrec.log.LGM_cur", "netDiv_CWM_std_tw_rs_1")){
    quantile.nti.XY.i <- quantile.XY(
      n.classes = 100,
      Y.var = "nti",
      X.var = X.var_i,
      SamplingPool.i = SamplingPool_i,
      dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div
    )
    
    quantile.nti.XY.all <- rbind(quantile.nti.XY.all, 
                             quantile.nti.XY.i)
  }
}

# as.factor(quantile.nti.XY.all$SamplingPool)

quantile.nti.XY.all$SamplingPool <- factor(quantile.nti.XY.all$SamplingPool, 
                                       levels(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div$SamplingPool))

# quantile.nti.XY.all$x.lab <- factor(quantile.nti.XY.all$x.lab , 
#                                 levels = c("diff.AnnTemp.LGM_cur", "diff.log.AnnPrec.log.LGM_cur", "netDiv_CWM_std_tw_rs_1"),
#                                 )
quantile.nti.XY.all$x.lab. <- quantile.nti.XY.all$x.lab 

quantile.nti.XY.all$x.lab <- factor(quantile.nti.XY.all$x.lab., 
                                levels = c("diff.AnnTemp.LGM_cur", "diff.log.AnnPrec.log.LGM_cur", "netDiv_CWM_std_tw_rs_1"),
                                labels = c(
                                  c(expression(atop(
                                    ~italic(Q)*"(Historical Change in Temperature) (°C)",
                                    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                  ))),
                                  c(expression(atop(
                                    ~italic(Q)*"(Historical Change in Precipitation)",
                                    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))))),
                                  c(expression(atop(
                                    ~italic(Q)*"(Local Net Diversification Rate)",
                                    lambda[CWM[tw]] - mu[CWM[tw]])))
                                )
)
# 


quantile.nti.XY.all$x.lab <- factor(quantile.nti.XY.all$x.lab., 
                                levels = c("diff.AnnTemp.LGM_cur", "diff.log.AnnPrec.log.LGM_cur", "netDiv_CWM_std_tw_rs_1"),
                                labels = c(
                                  c(expression(
                                    atop(
                                      "Historical Change in Temperature °C",
                                      scriptstyle("MAT"[Contemporary] ~ "\U2212" ~ "MAT"[LGM])
                                    )
                                  )
                                  ),
                                  c(expression(atop(
                                    "Historical Change in Precipitation",
                                    scriptstyle(log("MAP"[Contemporary]) ~ "\U2212" ~ log("MAP"[LGM]))))),
                                  c(expression(atop(
                                    "Local Net Diversification Rate",
                                    "\U03BB"[CWM[tw]] ~ "\U2212" ~ "\U03BC"[CWM[tw]])))
                                )
)
 

# quantile.nti.XY.all$x.lab <- factor(quantile.nti.XY.all$x.lab., 
#                                 levels = c("diff.AnnTemp.LGM_cur", "diff.log.AnnPrec.log.LGM_cur", "netDiv_CWM_std_tw_rs_1"),
#                                 labels = c(
#                                   c(expression(atop(
#                                     "Historical Change\n in Temperature (°C)",
#                                     scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
#                                   ))),
#                                   c(expression(atop(
#                                     "Historical Change\n in Precipitation",
#                                     scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))))),
#                                   c(expression(atop(
#                                     "Local Net Diversification Rate",
#                                     lambda[CWM[tw]] - mu[CWM[tw]])))
#                                 )
# )


labels_x_vars <- c(
  `diff.AnnTemp.LGM_cur` = c(expression(atop(
    "Historical Change\n in Temperature (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  `diff.log.AnnPrec.log.LGM_cur` = c(expression(atop(
    "Historical Change\n in Precipitation",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))))),
  `netDiv_CWM_std_tw_rs_1` = c(expression(atop(
    "Local Net Diversification Rate",
    ""[CWM[tw]] - mu[CWM[tw]])))
)

# factor(quantile.nti.XY.all$x.lab)
# 
# levels(quantile.nti.XY.all$x.lab) <- 
# 
# lab_x = c(expression(atop(
#   ~italic(Q)*"(Historical Change in Temperature) (°C)",
#   scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
# )))

facet_labels <- quantile.nti.XY.all %>% 
  distinct(x.lab, y.lab, SamplingPool, .keep_all = F)

facet_labels <- facet_labels %>%
  mutate(facet_labels = LETTERS[1:nrow(facet_labels)])

(fig.quantile.SamplingPool.NTI.diffTemp.log.diffPrec.netDiv <- ggplot(
  quantile.nti.XY.all,
  aes(
    x = mean.X,
    y = mean.Y
  )
) +
    geom_errorbar(aes(ymin = mean.Y - 1.96*se.Y,
                      ymax = mean.Y + 1.96*se.Y),
                  cex = 0.1) +
    # geom_errorbarh(aes(xmin = min.X,
    #                    xmax = max.X),
    #                cex = 1,
    #                alpha = 0.75) +
    # scico::scale_fill_scico_d(
    #   palette = "acton",
    #   direction = -1,
    #   drop = FALSE
    # ) +
    geom_point(cex = 2) +
    labs(
      y = "NTI",
      x = NULL
      # x = c(expression(atop(
      #   ~italic(Q)*"(Local Net Diversification Rate)",
      #   lambda[CWM[tw]] - mu[CWM[tw]])))
    ) +
    facet_grid(
      SamplingPool ~ x.lab,
      # SamplingPool ~  .,
      scales = "free",
      switch = "both",
      labeller = labeller(x.lab = as_labeller(labels_x_vars, 
                                              label_parsed))
    ) +
    geom_text(data = facet_labels,
              aes(x = -Inf, y = Inf, 
                  label = facet_labels, 
                  # group = facet_labels
              ),
              size = 6,
              hjust = -0.6,
              vjust = 1.5,
              inherit.aes = FALSE) +
    geom_hline(yintercept = 0, 
               alpha = 0.25) +
    geom_vline(xintercept = 0, 
               alpha = 0.25) +
    # scale_y_continuous(
    #   breaks = pretty(c(quantile.nti.XY.all$mean.Y,
    #                     0), n = 6),
    #   limits = c(-0.05, 0.05) + range(pretty(c(quantile.nti.XY.all$mean.Y, 0), n = 7)),
    #   #   expand = expansion(mult = c(0, 0)),
    #   position = "left"
    # ) +
    # scale_x_continuous(
    #   breaks = pretty(c(X, 0), n = 5),
    #   limits = c(-0.05, 0.05) + range(pretty(c(quantile.nti.XY.all$mean.X, 0),
    #                                          n = 6)),
  #   #        expand = expansion(mult = c(0, 0)),
  # ) +
  theme_classic(base_size = 15* 1.6) +
    theme(
      #text = element_text(family = "Karla", color = "#22211d"),
      legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_blank(),
      legend.title = element_blank(),
      legend.box = "horizontal",
      # axis.title.x = element_blank(),
      axis.text.x = element_text( #size = 15 * 1.6
        angle = 30,
        hjust = 1),
      axis.text.y = element_text( #size = 15 * 1.6
      ),
      axis.title.y = element_text(
        # size = 17 * 2,
        face = "bold"
      ),
      axis.title.x = element_text(
        # size = 14 * 1.4,
        face = "bold"
      ),
      plot.margin = unit(c(0.5, 0.6, 0, 0.6), "cm"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.switch.pad.grid = unit(0.5, "cm"),
      strip.text.x = element_text(size = 14),
      panel.border = element_rect(colour = "black", 
                                  fill = NA)
    ) +
    guides(colour = guide_legend(nrow = 1))
)


ragg::agg_png("fig.quantile.SamplingPool.NTI.diffTemp.log.diffPrec.netDiv.png",
              width = 10.5*1.3, 
              height = 18*1.3, 
              units = "in",
              res = 750)

fig.quantile.SamplingPool.NTI.diffTemp.log.diffPrec.netDiv 

dev.off()

# ggsave("fig.quantile.SamplingPool.NTI.diffTemp.log.diffPrec.netDiv.png", 
#        plot = fig.quantile.SamplingPool.NTI.diffTemp.log.diffPrec.netDiv , 
#        device = "png", 
#        width = 10.5*1.2, 
#        height = 18*1.2, 
#        units = "in",
#        dpi = 450) # units = "in"

