#################################################################################
### Code to compute independent binned bivariate kernel density estimates to  ###
### obtain the frequency probability in the association of community          ###
### phylogenetic structure – NRI and NTI – and historical climatic stability  ###
### (as measured by the changes in temperature and in precipitation since the ###
### last glacial maximum).                                                    ###
#                                                                               # 
#                                                                               # 
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2021-07-02"                                                     #
#                                                                               # 
#################################################################################

CWM.Div.MPD.MNTD.diff.worldClimate.Global.data <- CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
  filter(SamplingPool == "Global sampling") %>%
  left_join(MPD.MNTD.LatLong.diff.Env.Global %>%
              select(ID, 
                     diff.AnnTemp.mP2_cur:diff.AnnPrec.LGM_cur), 
            by = "ID") %>%
  mutate(diff.AnnTemp.LGM_cur = diff.AnnTemp.LGM_cur/10,
         # diff.log.AnnPrec.LGM_cur = diff.AnnPrec.LGM_cur
  )

# NRI and Annual Temperature

Global.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks::Hpi(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                 select(NRI, diff.AnnTemp.LGM_cur) %>%
                                                 drop_na(),
                                               pilot = "samse")

Global.NRI.diff.AnnTemp.LGM_cur.Hpi

# Compute kde 
kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks::kde(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                     select(diff.AnnTemp.LGM_cur, NRI) %>%
                                                     drop_na(), 
                                                   H = Global.NRI.diff.AnnTemp.LGM_cur.Hpi,
                                                   binned = TRUE,
                                                   compute.cont = TRUE)

# Represent with ks::plot.kde()

# "cont" specifies the density contours, which are upper percentages of highest
# density regions. The default contours are at 25%, 50%, and 75%

percentiles <- percent(c(0.01, seq(0.1, 0.9, length.out = 5), 0.99), 
                       accuracy = 1)

plot(kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi, 
     display = "filled.contour2", 
     cont = as.numeric(sub("%", "", percentiles, 
                           fixed = TRUE)))

# plot(kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi, display = "persp")

## Subset and prepare information to represent

fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi <- kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi

dimnames(fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi[['estimate']]) <- list(fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi[["eval.points"]][[1]], 
                                                                             fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi[["eval.points"]][[2]])

molten.fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi <- reshape2::melt(fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi[['estimate']],
                                                                      value.name = "estimate") %>%
  rename(diff.AnnTemp.LGM_cur = Var1,
         NRI = Var2)

# Represent contours with ggplot() and geom_contour

(fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi <- ggplot(molten.fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi, 
                                                       aes(x = diff.AnnTemp.LGM_cur, 
                                                           y = NRI)) +
    geom_contour_filled(aes(z = estimate), 
                        breaks = fhat.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi[["cont"]][percentiles]) +
    # scale_color_continuous_sequential(palette = "Reds", l1 = 20, c2 = 70, p1 = 1) +
    # scale_fill_brewer(palette = "BuPu") +
    # scale_fill_viridis_d(option = "magma", direction = -1) +
    # paletteer::scale_fill_paletteer_d("viridis::plasma") +
    # scale_fill_gradient(low="blue", high="red") +
    # scale_fill_manual(values = muted(scales::seq_gradient_pal("blue", "red", "Lab")(seq(0, 1, length.out = 7)))) +
    scico::scale_fill_scico_d(palette = "acton", direction = -1) + 
    labs(x = c(expression(atop("Historical Change in Temperature (°C)", 
                               scriptstyle("MAT"[Contemporary]-"MAT"[LGM])))),
         y = c(expression("NRI"["Global"]))) +
    geom_hline(yintercept = 0, alpha = 0.25) +  
    scale_y_continuous(#breaks = round(seq(min(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NRI, na.rm = TRUE), 
      #                    max(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NRI, na.rm = TRUE), length.out = 6), 1),
      breaks = pretty_breaks(7, 
                             high.u.bias = 0.1,
                             eps.correct = 0, 
                             min.n = 6),
      limits = c(-0.5,0.5) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NRI, 
                                   na.rm = TRUE)
    ) +
    scale_x_continuous(limits = c(-0.1,0.1) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$diff.AnnTemp.LGM_cur, 
                                                    na.rm = TRUE)) +
    theme_classic(base_size = 17) +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal",
          # axis.title.x = element_blank(),
          # axis.text.x = element_blank(),
          # axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16,
                                    face="bold")) + 
    guides(colour = guide_legend(nrow = 1)
    )
)

##


# NRI and Annual Precipitation

Global.NRI.diff.AnnPrec.LGM_cur.Hpi <- ks::Hpi(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                 select(NRI,
                                                        diff.AnnPrec.LGM_cur) %>%
                                                 drop_na(),
                                               pilot = "samse")

Global.NRI.diff.AnnPrec.LGM_cur.Hpi

CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
  select(NRI) %>%
  drop_na() %>%
  range()

# Compute kde 
kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi <- ks::kde(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                     select(diff.AnnPrec.LGM_cur, 
                                                            NRI) %>%
                                                     drop_na(), 
                                                   H = Global.NRI.diff.AnnPrec.LGM_cur.Hpi,
                                                   # xmin = c(-4978.519, -1.894606),
                                                   # xmax = c(1000.534, 18.352332),
                                                   binned = TRUE,
                                                   compute.cont = TRUE)


# (evpts <- do.call(expand.grid,  lapply(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
#                                         select(diff.AnnPrec.LGM_cur, NRI) %>%
#                                         drop_na(), 
#                                       quantile, 
#                                       prob=c(0.1, .25, .5, .75, .99)) ))



# Represent with ks::plot.kde()

# "cont" specifies the density contours, which are upper percentages of highest
# density regions. The default contours are at 25%, 50%, and 75%

percentiles <- percent(c(0.01, seq(0.1, 0.9, length.out = 5), 0.99), 
                       accuracy = 1)

plot(kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi, 
     display = "filled.contour2", 
     cont = as.numeric(sub("%", "", percentiles, 
                           fixed = TRUE)))

# plot(kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi, display = "persp")

## Subset and prepare information to represent

fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi <- kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi

dimnames(fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi[['estimate']]) <- list(fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi[["eval.points"]][[1]], 
                                                                             fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi[["eval.points"]][[2]])

molten.fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi <- reshape2::melt(fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi[['estimate']],
                                                                      value.name = "estimate") %>%
  rename(diff.AnnPrec.LGM_cur = Var1,
         NRI = Var2)

# Represent contours with ggplot() and geom_contour

(fig.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi <- ggplot(molten.fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi,
                                                       aes(x = diff.AnnPrec.LGM_cur, 
                                                           y = NRI)) +
    geom_contour_filled(aes(z = estimate), 
                        breaks = fhat.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi[["cont"]][percentiles]) + 
    # scale_fill_brewer(palette = "OrRd") +
    scico::scale_fill_scico_d(palette = "acton", direction = -1) + 
    labs(x = c(expression(atop("Historical Change in Precipitation (mm)", 
                               scriptstyle("MAP"[Contemporary]-"MAP"[LGM])))),
         y = c(expression("NRI"["Global"]))) +
    geom_hline(yintercept = 0, alpha = 0.25) +  
    scale_y_continuous(#breaks = round(seq(min(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NRI, na.rm = TRUE), 
      #                    max(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NRI, na.rm = TRUE), length.out = 6), 1),
      breaks = pretty_breaks(7, 
                             high.u.bias = 0.1,
                             eps.correct = 0, 
                             min.n = 6),
      limits = c(-0.5,0.5) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NRI, na.rm = TRUE)
    ) +
    scale_x_continuous(limits = c(-0.01,0.01) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$diff.AnnPrec.LGM_cur, 
                                                      na.rm = TRUE)
    ) +
    theme_classic(base_size = 17) +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal",
          #        axis.title.x = element_blank(),
          #        axis.text.x = element_blank(),
          # axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16, face="bold")) + 
    guides(colour = guide_legend(nrow = 1)
    )
)

myBreaks <- function(x){
  breaks <- c(min(x), median(x), max(x))
  names(breaks) <- attr(breaks,"labels")
  breaks
}

# NTI and Annual Temperature

Global.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks::Hpi(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                 select(NTI, diff.AnnTemp.LGM_cur) %>%
                                                 drop_na(),
                                               pilot = "samse")

Global.NTI.diff.AnnTemp.LGM_cur.Hpi

# Compute kde 
kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks::kde(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                     select(diff.AnnTemp.LGM_cur, NTI) %>%
                                                     drop_na(), 
                                                   H = Global.NTI.diff.AnnTemp.LGM_cur.Hpi,
                                                   binned = TRUE,
                                                   compute.cont = TRUE)

# Represent with ks::plot.kde()

# "cont" specifies the density contours, which are upper percentages of highest
# density regions. The default contours are at 25%, 50%, and 75%

percentiles <- percent(c(0.01, seq(0.1, 0.9, length.out = 5), 0.99), 
                       accuracy = 1)

plot(kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi, 
     display = "filled.contour2", 
     cont = as.numeric(sub("%", "", percentiles, 
                           fixed = TRUE)))


# plot(kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi, display = "persp")

## Subset and prepare information to represent

fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi <- kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi

dimnames(fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi[['estimate']]) <- list(fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi[["eval.points"]][[1]], 
                                                                             fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi[["eval.points"]][[2]])

molten.fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi <- reshape2::melt(fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi[['estimate']],
                                                                      value.name = "estimate") %>%
  rename(diff.AnnTemp.LGM_cur = Var1,
         NTI = Var2)

# Represent contours with ggplot() and geom_contour

(fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi <- ggplot(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data, 
                                                       aes(x = diff.AnnTemp.LGM_cur, 
                                                           y = NTI)) +
    # geom_point() +
    geom_contour_filled(data = molten.fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi, 
                        aes(z = estimate), 
                        breaks = fhat.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi[["cont"]][percentiles]) + 
    # scale_fill_brewer(palette = "OrRd") +
    scico::scale_fill_scico_d(palette = "acton", direction = -1) + 
    labs(x = c(expression(atop("Historical Change in Temperature (°C)", 
                               scriptstyle("MAT"[Contemporary]-"MAT"[LGM])))),
         y = c(expression("NTI"["Global"]))) +
    geom_hline(yintercept = 0, alpha = 0.25) +  
    scale_y_continuous(breaks = round(seq(min(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NTI, 
                                              na.rm = TRUE), 
                                          max(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NTI, 
                                              na.rm = TRUE), length.out = 5),1),
                       limits = c(-0.5,0.5) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NTI, 
                                                    na.rm = TRUE)
    ) +
    scale_x_continuous(limits = c(-0.1,0.1) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$diff.AnnTemp.LGM_cur,
                                                    na.rm = TRUE)
    ) +
    theme_classic(base_size = 17) +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal",
          #        axis.title.x = element_blank(),
          #        axis.text.x = element_blank(),
          # axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16,
                                    face="bold")) + 
    guides(colour = guide_legend(nrow = 1)
    )
)

### NTI and Annual Precipitation ####

Global.NTI.diff.AnnPrec.LGM_cur.Hpi <- ks::Hpi(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                 select(NTI, diff.AnnPrec.LGM_cur) %>%
                                                 drop_na(),
                                               pilot = "samse")

Global.NTI.diff.AnnPrec.LGM_cur.Hpi

# Compute kde 
kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi <- ks::kde(x = CWM.Div.MPD.MNTD.diff.worldClimate.Global.data %>%
                                                     select(diff.AnnPrec.LGM_cur, NTI) %>%
                                                     drop_na(), 
                                                   H = Global.NTI.diff.AnnPrec.LGM_cur.Hpi,
                                                   binned = TRUE,
                                                   compute.cont = TRUE)

# Represent with ks::plot.kde()

# "cont" specifies the density contours, which are upper percentages of highest
# density regions. The default contours are at 25%, 50%, and 75%

percentiles <- percent(c(0.01, seq(0.1, 0.9, length.out = 5), 0.99), 
                       accuracy = 1)

plot(kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi,
     display = "filled.contour2", 
     drawpoints=FALSE,
     cont = as.numeric(sub("%", "", percentiles, 
                           fixed = TRUE)))

# plot(kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi, display = "persp")

## Subset and prepare information to represent

fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi <- kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi

dimnames(fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi[['estimate']]) <- list(fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi[["eval.points"]][[1]], 
                                                                             fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi[["eval.points"]][[2]])

molten.fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi <- reshape2::melt(fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi[['estimate']],
                                                                      value.name = "estimate") %>%
  rename(diff.AnnPrec.LGM_cur = Var1,
         NTI = Var2)

# Represent contours with ggplot() and geom_contour

(fig.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi <- ggplot(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data, 
                                                       aes(x = diff.AnnPrec.LGM_cur, 
                                                           y = NTI)) +
    # geom_point() +
    geom_contour_filled(data = molten.fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi,
                        aes(z = estimate), 
                        breaks = fhat.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi[["cont"]][percentiles]) + 
    # scale_fill_brewer(palette = "OrRd") +
    scico::scale_fill_scico_d(palette = "acton", direction = -1) + 
    labs(x = c(expression(atop("Historical Change in Precipitation", 
                               scriptstyle(log("MAP"[Contemporary])-log("MAP"[LGM]))))),
         y = c(expression("NTI"["Global"]))) +
    geom_hline(yintercept = 0, alpha = 0.25) +  
    scale_y_continuous(breaks = round(seq(min(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NTI, na.rm = TRUE), 
                                          max(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NTI, na.rm = TRUE), length.out = 5),1),
                       limits = c(-0.5,0.5) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$NTI, na.rm = TRUE)) +
    scale_x_continuous(limits = c(-0.01,0.01) + range(CWM.Div.MPD.MNTD.diff.worldClimate.Global.data$diff.AnnPrec.LGM_cur, na.rm = TRUE)) +
    theme_classic(base_size = 17) +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal",
          #        axis.title.x = element_blank(),
          #        axis.text.x = element_blank(),
          # axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16, face="bold")) + 
    guides(colour = guide_legend(nrow = 1)
    )
)


(fig.kde.Global.NRI.NTI.diffTemp.diffPrec <- grid.arrange(
  cbind(
    rbind(ggplotGrob(fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi +
                       ggtitle("a") +
                       theme(axis.title.x = element_blank())), 
          ggplotGrob(fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi), 
          size = "last"),
    rbind(ggplotGrob(fig.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi + theme(axis.title.x = element_blank())), 
          ggplotGrob(fig.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi),
          size = "last")
  ),
  # common.legend, 
  nrow = 2,
  heights = c(11, 1)
)
)

(fig.kde.Global.NRI.NTI.diffTemp.diffPrec <- ggarrange(
  ggplotGrob(fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi +
               theme(axis.title.x = element_blank())), 
  ggplotGrob(fig.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi + 
               theme(axis.title.x = element_blank())), 
  ggplotGrob(fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi),
  ggplotGrob(fig.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi),
  labels = c("A", "B", "C", "D"),
  nrow = 2,
  ncol = 2,
  #  heights = c(11, 1),
  widths = c(1, 1),
  align = "hv"
)
)

(fig.kde.Global.NRI.NTI.diffTemp.diffPrec <- ggarrange(
  arrangeGrob(fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi +
                theme(axis.title.x = element_blank())), 
  arrangeGrob(fig.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi + 
                theme(axis.title.x = element_blank())),
  arrangeGrob(fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi),
  arrangeGrob(fig.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi),
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  widths = c(1, 1),
  align = "v"
)
)

(fig.kde.Global.NRI.NTI.diffTemp.diffPrec <- ggarrange(
  fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()), 
  fig.kde.Global.NRI.diff.AnnPrec.LGM_cur.Hpi + 
    theme(axis.title.x = element_blank()),
  fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi,
  fig.kde.Global.NTI.diff.AnnPrec.LGM_cur.Hpi,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  widths = c(1, 1),
  heights = c(0.9, 1),
  align = "v"
)
)

ggsave(file="fig.kde.Global.NRI.NTI.diffTemp.diffPrec.png", 
       fig.kde.Global.NRI.NTI.diffTemp.diffPrec,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 200)

