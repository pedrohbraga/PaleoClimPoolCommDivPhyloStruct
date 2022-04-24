#################################################################################
### Code to test and represent the relationship between indices for communtiy   #
### phylogenetic structure (NRI and NTI) and community weighted diversification #
### rate means                                                                  #
#                                                                               #
# The code begins by fitting ordinary least square models and assesing model    #
# assumptions. Assumptions were violated in a series of models.                 #
#
# When representing the data, we noted that there was a non-stationary          #
# relationship between NRi and net diversification rates. We proceeded to apply #
# a quantile regression framework.                                              #
#                                                                               #
# The final model is a forth-order polynomial quantile regression.              #
#                                                                               #
# Model selection between second-order, third-order and fourth-order polynomial #
# was done using BIC.                                                           #
#                                                                               #
# This code also produces a table and coefficient estimates and goodness-of-fit #
# for the quantile regressions on the effects of the proportional quantiles.    #
#                                                                               #
#                                                                               #
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2022-02-28"                                                     #
#                                                                               # 
#################################################################################

### Temporary environment preparation for PPN

# Check for needed packages, and install the missing ones

required.libraries <- c("ggplot2",  
                        "viridis", 
                        "broom", 
                        "tidyr", 
                        "scales",
                        "kableExtra",
                        "tidyverse", 
                        "reshape2",
                        "quantreg",
                        "WRTDStidal",
                        "ggrepel")

needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]

if(length(needed.libraries)) install.packages(needed.libraries, Ncpus = 12)

# Load all required libraries at once

lapply(required.libraries, 
       require, 
       character.only = TRUE)

### Loading data

MPD.MNTD.LatLong.AllScales.rarefaction.relative.CWM.netDiv <- MPD.MNTD.LatLong.AllScales.rarefaction.relative %>%
  filter(SamplingPool == "Global sampling") %>%
  left_join(CWM.Div.MPD.Chiroptera.Comm %>%
              filter(SamplingPool == "Global sampling") %>%
              select(ID,
                     lambda_CWM:netDiv_CWM_std_w),
            by = c("ID" = "ID")) %>%
  drop_na(nri.sec) %>%
  drop_na(ID_Realm)

phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff <- MPD.MNTD.LatLong.AllScales.rarefaction.relative.CWM.netDiv %>%
  left_join(MPD.MNTD.LatLong.AllScales.rarefaction.relative.worldClimate.diff %>%
              filter(SamplingPool == "Global sampling") %>%
              select(ID,
                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur),
            by = c("ID" = "ID")) %>%
  select(ID,
         ID_SamplingPool,
         Longitude, Latitude,
         LONB:LATDD,
         SamplingPool,
         ID_Realm,
         ID_PlateName,
         ID_Biome,
         ID_Biome_Realm,
         ID_Ecoregion,
         ntaxa,
         mpd.obs.query:nri,
         ses.mpd.z.query.rarefac.mean:nri.rarefac.mean,
         mntd.obs.query:nti,
         ses.mntd.z.query.rarefac.mean:nti.rarefac.mean,
         ntaxa.pool.size,
         rarefac.pool.size,
         lambda_CWM:diff.log.AnnPrec.log.LGM_cur)

# write.csv(phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff,
#           "data/matrices/phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff.csv")

phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff <- read.csv("data/matrices/phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff.csv",
                                                                          h = T,
                                                                          row.names = 1)

phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$ID_Realm <- factor(phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$ID_Realm, 
                                                                                 levels = c('Neotropical',
                                                                                            'Nearctic',
                                                                                            'Afrotropical',
                                                                                            'Palearctic', 
                                                                                            'Indomalay', 
                                                                                            'Australasian'))
###################################
#####  Visualization: Part 1  #####
###################################

# Standard Scatterplot ----------------------------------------------------

# NRI and Net Diversification Rate ###

ggplot(data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
       aes(x = netDiv_CWM_std_tw, 
           y = nri.sec)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NRI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

# NTI and Net Diversification Rate ###

ggplot(data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
       aes(x = netDiv_CWM_std_tw, 
           y = nti.sec)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NTI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

#############################################################################################
### Representation and analyses for net diversification rates and phylogenetic structure ####
#############################################################################################

# Linear model ------------------------------------------------------------

lm.NRI.netDiv_CWM_std_tw <- lm(nri.sec ~ netDiv_CWM_std_tw, 
                               data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff
)


summary(lm.NRI.netDiv_CWM_std_tw)

# Verify assumptions

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(lm.NRI.netDiv_CWM_std_tw, which = 1)
plot(lm.NRI.netDiv_CWM_std_tw, which = 2)
plot(lm.NRI.netDiv_CWM_std_tw, which = 3)
plot(lm.NRI.netDiv_CWM_std_tw, which = 4)

dev.off()

# Diagnostic plots show a clear violation of linear assumptions

# Quadratic model ---------------------------------------------------------

lm.quad.NRI.netDiv_CWM_std_tw <- lm(nri.sec ~ netDiv_CWM_std_tw + I(netDiv_CWM_std_tw^2), 
                                    data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

summary(lm.quad.NRI.netDiv_CWM_std_tw)

# Verify assumptions

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(lm.quad.NRI.netDiv_CWM_std_tw, which=1)
plot(lm.quad.NRI.netDiv_CWM_std_tw, which=2)
plot(lm.quad.NRI.netDiv_CWM_std_tw, which=3)
plot(lm.quad.NRI.netDiv_CWM_std_tw, which=5)

dev.off()

# 4th-degree polynomial model ---------------------------------------------------------

lm.poly.4.NRI.netDiv_CWM_std_tw <- lm(nri.sec ~ poly(netDiv_CWM_std_tw, 4), 
                                      data =  phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

summary(lm.poly.4.NRI.netDiv_CWM_std_tw)

# Verify assumptions

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(lm.poly.4.NRI.netDiv_CWM_std_tw, which=1)
plot(lm.poly.4.NRI.netDiv_CWM_std_tw, which=2)
plot(lm.poly.4.NRI.netDiv_CWM_std_tw, which=3)
plot(lm.poly.4.NRI.netDiv_CWM_std_tw, which=5)

dev.off()

# LOESS -------------------------------------------------------------------

## Attempting with a locally estimated scatterplot smoothing (loess)

loess.NRI.netDiv_CWM_std_tw <- loess(nri.sec ~ netDiv_CWM_std_tw, 
                                     data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
                                     span = 0.50)

loess.NRI.netDiv_CWM_std_tw.fit <- augment(loess.NRI.netDiv_CWM_std_tw)

# get smoothed output

smoothed.loess.NRI.netDiv_CWM_std_tw <- predict(loess.NRI.netDiv_CWM_std_tw) 

# Representing LOESS ------------------------------------------------------

ggplot(data = loess.NRI.netDiv_CWM_std_tw.fit, 
       aes(x = netDiv_CWM_std_tw, 
           y = nri.sec)) +
  geom_point(alpha = 0.7, 
             size = 1.4) +
  geom_line(aes(y = .fitted), 
            color = "red") +
  # facet_grid(. ~ ID_Realm) +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NRI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

# Quantile regression -----------------------------------------------------

## Defining tau

# (q20 <- seq(0.05, 0.95, by = 0.05))

(q10 <- seq(0.05, 0.95, by = 0.1))

(q100 <- seq(0.01, 0.99, by = 0.01))


## Standard quantile regression fit

rq.NRI.netDiv_CWM_std_tw <- rq(nri.sec ~ netDiv_CWM_std_tw, 
                               tau = q10,
                               method = "br",
                               data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

rq.NRI.netDiv_CWM_std_tw

summary(rq.NRI.netDiv_CWM_std_tw)

# Plot distributions

plot(summary(rq.NRI.netDiv_CWM_std_tw,
             se = "nid"),
     main = c("Intercept", 
              expression("Net Diversification Rate"[CWM[STD[tw]]])))

# Represent in ggplot

# NRI and Net Diversification Rate ####

ggplot(data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
       aes(x = netDiv_CWM_std_tw, 
           y = nri.sec)) +
  geom_point(alpha = 0.15,
             size = 2,
             stroke = 0.01 #,
             #             aes(colour = factor(ID_Realm))
  ) +
  # scale_colour_viridis(option="D", 
  #                     begin = 0,
  #                     end = 1,
  #                     direction = 1, 
  #                     discrete = TRUE,
  #                     name = "Realm") +
  stat_quantile(quantiles = q10, 
                formula = y ~ x, 
                colour = "black"
  ) +
  # xlim(min(prebinning_test$netDiv_CWM_std_tw),
  #      0.35) +
  # # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NRI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

# NTI and Net Diversification Rate ####

ggplot(data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
       aes(x = netDiv_CWM_std_tw, 
           y = nti.sec)) +
  geom_point(alpha = 0.15,
             size = 2,
             stroke = 0.01 #,
             #             aes(colour = factor(ID_Realm))
  ) +
  # scale_colour_viridis(option="D", 
  #                     begin = 0,
  #                     end = 1,
  #                     direction = 1, 
  #                     discrete = TRUE,
  #                     name = "Realm") +
  stat_quantile(quantiles = q100, 
                formula = y ~ x, 
                colour = "black"
  ) +
  # xlim(min(prebinning_test$netDiv_CWM_std_tw),
  #      0.35) +
  # # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NTI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

## Polynomial quantile regression fit ----------------------

rq.poly.4.NRI.netDiv_CWM_std_tw <- rq(nri.sec ~ poly(netDiv_CWM_std_tw, 4), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

rq.poly.3.NRI.netDiv_CWM_std_tw <- rq(nri.sec ~ poly(netDiv_CWM_std_tw, 3), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

rq.poly.2.NRI.netDiv_CWM_std_tw <- rq(nri.sec ~ poly(netDiv_CWM_std_tw, 2), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)


## Calculating the BIC

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw$fitted.values)))

AIC(rq.poly.3.NRI.netDiv_CWM_std_tw, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw$fitted.values)))

AIC(rq.poly.4.NRI.netDiv_CWM_std_tw, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw$fitted.values)))

anova(lm.poly.4.NRI.netDiv_CWM_std_tw,
      rq.poly.4.NRI.netDiv_CWM_std_tw[[1]])

anova.rq(rq.poly.4.NRI.netDiv_CWM_std_tw,
         rq.poly.3.NRI.netDiv_CWM_std_tw,
         rq.poly.2.NRI.netDiv_CWM_std_tw, 
         joint = FALSE)

class(rq.poly.4.NRI.netDiv_CWM_std_tw)

KhmaladzeTest.rq.poly.4.NRI.netDiv_CWM_std_tw <- KhmaladzeTest(nri.sec ~ poly(netDiv_CWM_std_tw, 4), 
                                                               tau = q10,
                                                               method = "br",
                                                               # ci = TRUE,
                                                               data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

# Create Part 1 of Table S3

# Create a formatted table containing the model coefficients, which will go in
# the supplementary material.

# In the original submission, this table corresponds to Table S3.

# This table is generated with stargazer::stargazer(), and because it estimates
# the SE using bootstrapping, it may take up to 30 min to run.

# I was able to reduce time by using trace(".stargazer.wrap", edit = TRUE), and
# modifying the summary part for the rq.sq to change the number of replicates in
# the bootstrapping method from 200 (default) to 100.

# To do this, look for the code below, and then add "R = 100" within summary()

# if (.model.identify(object.name) == "rq") {
#   .summary.object <<- suppressMessages(summary(object.name, se=.format.rq.se))
# }

rq.poly.4.NRI.netDiv_CWM_std_tw.mapped <- map(q10, 
                                              ~rq(nri.sec ~ poly(netDiv_CWM_std_tw, 4), 
                                                  tau = .x,
                                                  method = "br",
                                                  data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)
)

rq.poly.4.NRI.netDiv_CWM_std_tw.formatted.coef.table <- stargazer::stargazer(lm.poly.4.NRI.netDiv_CWM_std_tw,
                                                                             rq.poly.4.NRI.netDiv_CWM_std_tw.mapped, 
                                                                             rq.se = "boot",
                                                                             column.labels = c("OLS", paste("tau = ", q10)),
                                                                             covariate.labels = NULL,
                                                                             dep.var.labels = "NRI",
                                                                             omit = c("Constant"),
                                                                             model.numbers = TRUE,
                                                                             model.names =  FALSE,
                                                                             keep.stat = c('n'),
                                                                             type ='html')


# Table S3 - Part A

write(rq.poly.4.NRI.netDiv_CWM_std_tw.formatted.coef.table, 
      "rq.poly.4.NRI.netDiv_CWM_std_tw.formatted.coef.table.html")

# pandoc_convert("rq.formatted.table.tex", to = "markdown")

# Calculate quantile regression goodness of fit using residuals and non-conditional residuals

# The goodness of fit measure for quantile regression is estimated as 1 minus
# the ratio between the sum of absolute deviations in the fully parameterized
# models and the sum of absolute deviations in the null (non-conditional)
# quantile model. The values are useful for comparisons between quantile models,
# but they are not comparable to standard coefficients of determination. The
# latter is based on the variance of squared deviations, whereas goodness of fit
# values for quantile regression are based on absolute deviations. Goodness of
# fit values will always be smaller than R2 values.

# Koenker, R., Machado, J.A.F. 1999. Goodness of fit and related inference
# processes for quantile regression. Journal of the American Statistical
# Association. 94(448):1296-1310.

# ?WRTDStidal::goodfit()

goodfit.rq <- vector(mode = "numeric", length(q10))

names(goodfit.rq) <- q10

for(i in 1:length(q10)){
  
  ## conditional quantile model
  rq.poly.4.NRI.netDiv_CWM_std_tw.i <- rq(nri.sec ~ poly(netDiv_CWM_std_tw, 4), 
                                          tau = q10[i],
                                          #  method = "conquer",
                                          #  ci = TRUE,
                                          data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)
  
  resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i <- resid(rq.poly.4.NRI.netDiv_CWM_std_tw.i)
  
  ## non-conditional quantile model
  
  null.rq.poly.4.NRI.netDiv_CWM_std_tw.i <- rq(nri.sec ~ 1,
                                               tau = q10[i],
                                               # method = "conquer",
                                               # ci = TRUE,
                                               data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)
  
  null.resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i <- resid(null.rq.poly.4.NRI.netDiv_CWM_std_tw.i)
  
  ## quantile regression goodness of fit
  
  goodfit.rq.i <- goodfit(resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i, 
                          null.resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i, 
                          q10[i])
  
  (goodfit.rq[i] <- goodfit.rq.i)
}

goodfit.rq

# Create Part B of Table S3

kable(goodfit.rq, format = "markdown")

### Representing in ggplot2 ####

# Obtain predicted values for the same x

pred.rq.poly.4.NRI.netDiv_CWM_std_tw <- as.data.frame(
  predict(rq.poly.4.NRI.netDiv_CWM_std_tw, 
          newdata = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff %>% 
            select(netDiv_CWM_std_tw)))

pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw <- phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$netDiv_CWM_std_tw

# Reshape predictions to long format

pred.rq.poly.4.NRI.netDiv_CWM_std_tw<- reshape2::melt(pred.rq.poly.4.NRI.netDiv_CWM_std_tw, 
                                                      id.var = "netDiv_CWM_std_tw",
                                                      variable.name = "tau",
                                                      value.name = "nri.sec")

pred.rq.poly.4.NRI.netDiv_CWM_std_tw$tau <- as.numeric(gsub("tau= (.*)", 
                                                            "\\1", 
                                                            pred.rq.poly.4.NRI.netDiv_CWM_std_tw$tau))
pred.rq.poly.4.NRI.netDiv_CWM_std_tw$label <- NA
pred.rq.poly.4.NRI.netDiv_CWM_std_tw$label[which(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw == max(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw))] <- pred.rq.poly.4.NRI.netDiv_CWM_std_tw$tau[which(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw == max(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw))]

# Represent figure

(NRI.poly.4.netDiv_CWM_std_tw.rq.plot <- ggplot(data = pred.rq.poly.4.NRI.netDiv_CWM_std_tw, 
                                                aes(x = netDiv_CWM_std_tw, 
                                                    y = nri.sec,
                                                    colour = tau, 
                                                    group = factor(tau)))  +   
    geom_point(data= phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff,
               aes(x = netDiv_CWM_std_tw, 
                   y = nri.sec),
               alpha = 0.15) + 
    geom_line(size = 1) + 
    geom_text_repel(
      data = pred.rq.poly.4.NRI.netDiv_CWM_std_tw %>%
        drop_na() %>%
        group_by(tau) %>%
        summarise(tau = mean(tau),
                  nri.sec = mean(nri.sec),
                  netDiv_CWM_std_tw = mean(netDiv_CWM_std_tw),
                  label = mean(label)
        ),
      aes(label = paste(tau*100, "%")),
      size = 4.5,
      nudge_x = 0.05,
      segment.color = "grey",
      min.segment.length = 0.1,
      colour = "black"
    )  +
    geom_hline(yintercept = 0, alpha = 0.25) +  
    scale_x_continuous( limits = c(c(-0.025, 0.05) + 
                                     range(phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$netDiv_CWM_std_tw))) +
    scico::scale_colour_scico(palette = "acton", direction = -1) + 
    # facet_grid(. ~ ID_Realm,
    #  scales = "free") +
    labs(y = c(expression("NRI"["Global"])),
         x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
    theme_classic(base_size = 20) +
    theme(legend.position = "none", 
          legend.direction = "horizontal",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 9),
          legend.title = element_text(face = "bold", size = 9),
          legend.box = "horizontal") + 
    guides(colour = guide_legend(nrow = 1))
)



## Ignore ###

library(quantreg)
library(splines)

set.seed(142857) # Reproducibility
x <- phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$netDiv_CWM_std_tw
d <- data.frame(x = x, y = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$nri)

# Plots
plot(y~x, data = d, type = "n", las = 1)
points(d$x, d$y, cex = 0.5, col = "#0000004C")
# Fit the quantile regressions

for (tau in 1:4/5) {
  rqfit <- rq(y ~ ns(x, 
                     df = 5), 
              tau = tau, data = d)
  xx <- phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$netDiv_CWM_std_tw
  fit <- predict(rqfit, 
                 newdata = 
                   data.frame(x = xx))
  lines(fit ~ xx, 
        col = "red", 
        lwd = 1.5)
}


# Represent in ggplot

# NRI and Net Diversification Rate ####

ggplot(data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
       aes(x = netDiv_CWM_std_tw, 
           y = nri.sec)) +
  geom_point(alpha = 0.15,
             size = 2,
             stroke = 0.01 #,
             #             aes(colour = factor(ID_Realm))
  ) +
  # scale_colour_viridis(option="D", 
  #                     begin = 0,
  #                     end = 1,
  #                     direction = 1, 
  #                     discrete = TRUE,
  #                     name = "Realm") +
  stat_quantile(quantiles = q10, 
                formula = y ~ ns(x, 
                                 df = 5), 
                colour = "black"
  ) +
  # xlim(min(prebinning_test$netDiv_CWM_std_tw),
  #      0.35) +
  # # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NRI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))



# NRI and Net Diversification Rate ####

ggplot(data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff, 
       aes(x = diff.AnnTemp.LGM_cur, 
           y = netDiv_CWM_std_tw)) +
  geom_point(alpha = 0.15,
             size = 2,
             stroke = 0.01 #,
             #             aes(colour = factor(ID_Realm))
  ) +
  # scale_colour_viridis(option="D", 
  #                     begin = 0,
  #                     end = 1,
  #                     direction = 1, 
  #                     discrete = TRUE,
  #                     name = "Realm") +
  stat_quantile(quantiles = q10, 
                formula = y ~ x, 
                colour = "black"
  ) +
  # xlim(min(prebinning_test$netDiv_CWM_std_tw),
  #      0.35) +
  # # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("NRI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

###


rq.poly.4.NRI.netDiv_CWM_std_tw <- rq(nri.sec ~ poly(netDiv_CWM_std_tw, 4), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)

rq.poly.4.NRI.diff.AnnTemp.LGM_cur <- rq(nri.sec ~ poly(diff.AnnTemp.LGM_cur, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff)


resid(rq.poly.4.NRI.diff.AnnTemp.LGM_cur)


nri.sec.resid.rq.poly.4.netDiv.diff.AnnTemp.LGM_cur <- data.frame(
  nri.sec = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff$nri.sec, 
  resid.netDiv.diff.AnnTemp = resid(rq(netDiv_CWM_std_tw ~ poly(diff.AnnTemp.LGM_cur, 4), 
                                       tau = q10,
                                       method = "br",
                                       # ci = TRUE,
                                       data = phyloStr.LatLong.AllScales.raref.rel.CWM.netDiv.ClimStab.diff))[, 9]
  )

#Now generate your plot using this NEWDF as your data

ggplot(data = nri.sec.resid.rq.poly.4.netDiv.diff.AnnTemp.LGM_cur, 
       aes(x = resid.netDiv.diff.AnnTemp, 
           y = nri.sec)) + 
  geom_point(na.rm = T) +
  labs(x = c(expression("E[Net Diversification Rate | Historical Climatic Stability]")),
       y = c(expression("NRI"[Global]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

