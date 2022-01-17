#################################################################################
### Code to obtain local diversification rates for each community based on    ### 
### rates estimated for each tip in the phylogeny                             ###
#                                                                               #
# Current measures calculated:                                                  #
#                                                                               #
# 1. arithmeticMeans.lambda, arithmeticMeans.mu and arithmeticMeans.netdiv:     #
#                                                                               #
# For the arithmetic mean, the rate for the kth grid cell is computed as        #
# λk = (Σwiλi)/Σwi in which Σ denotes a summation over all N species (i = 1 to  #
# i = N) present in the kth cell, and wi is the weight assigned to the ith      #
# species. We computed weights for each species as the inverse of the number of #
# grid cells in which the species was found                                     #
#                                                                               #
# 2. lambda.CWM.Chiroptera.Comm, mu.CWM.Chiroptera.Comm,                  
# netdiv.CWM.Chiroptera.Comm:
# 
# 
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2021-05-26"                                                     #
#                                                                               # 
#################################################################################

# Check lengths
length(tip.rates$lambda.avg); length(tip.rates$mu.avg)

### Computing arithmetic weighted means for diversification rates

# For the arithmetic mean, the rate for the kth grid cell is computed as
# λk=(Σwiλi)/Σwi in which Σ denotes a summation over all N species (i=1 to
# i=N) present in the kth cell, and wi is the weight assigned to the ith
# species. We computed weights for each species as the inverse of the number of
# grid cells in which the species was found21.

### For speciation rates (lambda)

dim(Chiroptera.FaurSven.comm)
lambda.Chiroptera.Comm <- Chiroptera.FaurSven.comm

weightSpecies <- colSums(Chiroptera.FaurSven.comm)
weightSpecies <- 1/colSums(Chiroptera.FaurSven.comm)

for(i in names(tip.rates$lambda.avg)){
  lambda.Chiroptera.Comm[, i] <- Chiroptera.FaurSven.comm[, i] * tip.rates$lambda.avg[i]
}

arithmeticMeans.lambda <- matrix(NA, nrow(Chiroptera.FaurSven.comm), 1)

rownames(arithmeticMeans.lambda) <- rownames(Chiroptera.FaurSven.comm)

for(i in 1:nrow(Chiroptera.FaurSven.comm)){
  subset.lambda.Chiroptera.Comm <- lambda.Chiroptera.Comm[i, lambda.Chiroptera.Comm[i, ] != 0, drop = FALSE]
  if(ncol(subset.lambda.Chiroptera.Comm) != 0) {
    arithmeticMeans.lambda[i, ] <- rowSums(subset.lambda.Chiroptera.Comm * weightSpecies[colnames(subset.lambda.Chiroptera.Comm)])/sum(weightSpecies[colnames(subset.lambda.Chiroptera.Comm)])
  } else
    arithmeticMeans.lambda[i, ] <- 0
}

# arithmeticMeans.lambda

colnames(arithmeticMeans.lambda) <- "arithmeticMeans.lambda"

mean(lambda.Chiroptera.Comm[lambda.Chiroptera.Comm != 0, ], 
     na.rm = TRUE)

### For extinction rates (mu)

dim(Chiroptera.FaurSven.comm)
mu.Chiroptera.Comm <- Chiroptera.FaurSven.comm

for(i in names(tip.rates$mu.avg)){
  mu.Chiroptera.Comm[, i] <- Chiroptera.FaurSven.comm[, i] * tip.rates$mu.avg[i]
}

arithmeticMeans.mu <- matrix(NA, nrow(Chiroptera.FaurSven.comm), 1)
rownames(arithmeticMeans.mu) <- rownames(Chiroptera.FaurSven.comm)

for(i in 1:nrow(Chiroptera.FaurSven.comm)){
  subset.mu.Chiroptera.Comm <- mu.Chiroptera.Comm[i, mu.Chiroptera.Comm[i, ] != 0, 
                                                  drop = FALSE]
  if(ncol(subset.mu.Chiroptera.Comm) != 0) {
    arithmeticMeans.mu[i, ] <- rowSums(subset.mu.Chiroptera.Comm * weightSpecies[colnames(subset.mu.Chiroptera.Comm)])/sum(weightSpecies[colnames(subset.mu.Chiroptera.Comm)])
  } else
    arithmeticMeans.mu[i, ] <- 0
}

# arithmeticMeans.mu

colnames(arithmeticMeans.mu) <- "arithmeticMeans.mu"

# arithmeticMeans.mu

## For local net diversification rates

netDiv.arithmeticMeans.Chiroptera.Comm <- arithmeticMeans.lambda - arithmeticMeans.mu

colnames(netDiv.arithmeticMeans.Chiroptera.Comm) <- "arithmeticMeans.netDiv"

#### Calculating CWM for Net Diversification Rates ####

### CWM for Speciation Rates

lambda.CWM.Chiroptera.Comm <- CWM_Std_TW(Trait = tip.rates$lambda.avg, 
                                         Distrib = Chiroptera.FaurSven.comm)

colnames(lambda.CWM.Chiroptera.Comm) <- paste("lambda", 
                                              colnames(lambda.CWM.Chiroptera.Comm), 
                                              sep = "_")

mu.CWM.Chiroptera.Comm <- CWM_Std_TW(Trait = tip.rates$mu.avg, 
                                     Distrib = Chiroptera.FaurSven.comm)

colnames(mu.CWM.Chiroptera.Comm) <- paste("mu", 
                                          colnames(mu.CWM.Chiroptera.Comm), 
                                          sep = "_")

## Subtract mu from lambda to obtain per community net diversification rates
# (note that this is different then subtracting mu from lambda directly at the
# tip-rates)

netDiv.CWM.Chiroptera.Comm <- lambda.CWM.Chiroptera.Comm - mu.CWM.Chiroptera.Comm

colnames(netDiv.CWM.Chiroptera.Comm) <- c("netDiv_CWM", 
                                          "netDiv_CWM_std_tw", 
                                          "netDiv_CWM_std_w" )

netDiv.CWM.Chiroptera <- data.frame(lambda.CWM.Chiroptera.Comm,
                                    mu.CWM.Chiroptera.Comm, 
                                    netDiv.CWM.Chiroptera.Comm) %>% 
  rownames_to_column("ID")

netDiv.CWM.Chiroptera$ID <- as.numeric(netDiv.CWM.Chiroptera$ID)

CWM.Div.MPD.Chiroptera.Comm <- merge(x = MPD.LatLong.Env.AllScales %>%
                                       select(-c(bio_1:bio_19)), 
                                     y =  netDiv.CWM.Chiroptera, 
                                     by = "ID", 
                                     all.x = TRUE)

CWM.Div.MPD.MNTD.Chiroptera.Comm <- merge(x = MPD.MNTD.LatLong.AllScales, 
                                          y =  netDiv.CWM.Chiroptera, 
                                          by = "ID", 
                                          all.x = TRUE)

Div.MPD.Chiroptera.Comm <- merge(x = MPD.LatLong.Env.AllScales[MPD.LatLong.Env.AllScales$SamplingPool == "Global sampling", ], 
                                 y =  netDiv.arithmeticMeans.Chiroptera.Comm, 
                                 by = 0, 
                                 all.x = TRUE)

head(CWM.Div.MPD.Chiroptera.Comm)

dim(CWM.Div.MPD.Chiroptera.Comm)

# write.csv(CWM.Div.MPD.Chiroptera.Comm, "data/matrices/CWM.Div.MPD.Chiroptera.Comm.SamplingCorrected_2.csv")

# Explore units that have more than 2 co-occurring taxa.

CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2)

# CWM.Div.MPD.Chiroptera.Comm <- read.csv("data/matrices/CWM.Div.MPD.Chiroptera.Comm.csv", h = T)

CWM.Div.MPD.Chiroptera.Comm$ID_Realm <- factor(CWM.Div.MPD.Chiroptera.Comm$ID_Realm, 
                                               levels=c('Neotropical',
                                                        'Nearctic',
                                                        'Afrotropical',
                                                        'Palearctic', 
                                                        'Indomalay', 
                                                        'Australasian'))

CWM.Div.MPD.Chiroptera.Comm$SamplingPool <-  factor(CWM.Div.MPD.Chiroptera.Comm$SamplingPool,  
                                                    levels = c("Global sampling",
                                                               "Hemispheric sampling",
                                                               "Realm sampling",
                                                               "Plate sampling",
                                                               "Biome sampling",
                                                               "Ecoregion sampling"))


#############################################################################################
### Representation and analyses for net diversification rates and phylogenetic structure ####
#############################################################################################


# Linear model ------------------------------------------------------------

lm.NRI.netDiv_CWM_std_tw <- lm(NRI ~ netDiv_CWM_std_tw, 
                               data = CWM.Div.MPD.Chiroptera.Comm %>%
                                 filter(SamplingPool == "Global sampling")
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

lm.quad.NRI.netDiv_CWM_std_tw <- lm(NRI ~ netDiv_CWM_std_tw + I(netDiv_CWM_std_tw^2), 
                                    data = CWM.Div.MPD.Chiroptera.Comm %>%
                                      filter(SamplingPool == "Global sampling"))

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

# Polynomial model ---------------------------------------------------------

lm.poly.4.NRI.netDiv_CWM_std_tw <- lm(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                      data = CWM.Div.MPD.Chiroptera.Comm %>%
                                        filter(SamplingPool == "Global sampling") %>%
                                        drop_na(NRI))

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

###################################
#####  Visualization: Part 1  #####
###################################


# Standard Scatterplot ----------------------------------------------------

# NRI
ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


# NTI

ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_tw, 
           y = NTI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

# Scatterplot with smooth lines -------------------------------------------

# NRI

ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  stat_smooth(method = "gam",
              formula = y ~ poly(x, 4), 
              se = TRUE,
              colour = "black") +
  #  stat_smooth(method = "lm", formula= (y ~ poly(x, 2)), se = TRUE, colour = "black") +
  #  stat_smooth(method = "loess", 
  #              formula = y ~ poly(x, 2), size = 1)
  facet_grid(. ~ ID_Realm) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


# NTI

ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_tw, 
           y = NTI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  stat_smooth(method = "gam",
              formula = y ~ poly(x, 4), 
              se = TRUE,
              colour = "black") +
  #  stat_smooth(method = "lm", formula= (y ~ poly(x, 2)), se = TRUE, colour = "black") +
  #  stat_smooth(method = "loess", 
  #              formula = y ~ poly(x, 2), size = 1)
  facet_grid(. ~ ID_Realm) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


## Attempting with a locally estimated scatterplot smoothing (loess)

loess.NRI.netDiv_CWM_std_tw <- loess(NRI ~ netDiv_CWM_std_tw, 
                                     data = CWM.Div.MPD.Chiroptera.Comm %>%
                                       filter(SamplingPool == "Global sampling"), 
                                     span = 0.50)

loess.NRI.netDiv_CWM_std_tw.fit <- augment(loess.NRI.netDiv_CWM_std_tw)

# get smoothed output

smoothed.loess.NRI.netDiv_CWM_std_tw <- predict(loess.NRI.netDiv_CWM_std_tw) 


# Representing LOESS ------------------------------------------------------

ggplot(loess.NRI.netDiv_CWM_std_tw.fit, 
       aes(netDiv_CWM_std_tw, NRI)) +
  geom_point() +
  geom_line(aes(y = .fitted), 
            color = "red")

ggplot(data = loess.NRI.netDiv_CWM_std_tw.fit, 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4) +
  geom_line(aes(y = .fitted), 
            color = "red") +
  # facet_grid(. ~ ID_Realm) +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


# Quantile regression -----------------------------------------------------

# Preparing dataset

prebinning_test <- CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
  filter(SamplingPool == "Global sampling") %>%
  dplyr::select(c(ID, ID_Realm, netDiv_CWM_std_tw, NRI, NTI)) %>%
  drop_na() %>%
  as_tibble()

(q20 <- seq(0.05, 0.95, by = 0.05))

(q10 <- seq(0.05, 0.95, by = 0.1))

(q7 <- seq(0.05, 0.95, by = 0.15))

## Standard quantile regression fit

# NRI

rq.NRI.netDiv_CWM_std_tw <- rq(NRI ~ netDiv_CWM_std_tw, 
                               tau = q10,
                               method = "br",
                               data = prebinning_test)

rq.NRI.netDiv_CWM_std_tw

summary(rq.NRI.netDiv_CWM_std_tw)

# Represent 

plot(rq.NRI.netDiv_CWM_std_tw)

# Plot distributions

plot(summary(rq.NRI.netDiv_CWM_std_tw,
             se = "nid"),
     main = c("Intercept", 
              expression("Net Diversification Rate"[CWM[STD[tw]]])))

## Polynomial quantile regression fit ----------------------

# Worldwide model

rq.poly.5.NRI.netDiv_CWM_std_tw <- rq(NRI ~ poly(netDiv_CWM_std_tw, 5), 
                                      tau = q20,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

rq.poly.4.NRI.netDiv_CWM_std_tw <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

rq.poly.3.NRI.netDiv_CWM_std_tw <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

rq.poly.2.NRI.netDiv_CWM_std_tw <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw, k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw$fitted.values)))
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw, k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw$fitted.values)))
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw, k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw$fitted.values)))

anova(rq.poly.3.NRI.netDiv_CWM_std_tw, joint = FALSE)

anova(rq.poly.4.NRI.netDiv_CWM_std_tw,
      rq.poly.3.NRI.netDiv_CWM_std_tw,
      rq.poly.2.NRI.netDiv_CWM_std_tw[[1]], joint = FALSE)

## Obtaining the formatted table of coefficients for the preferred model --------

# Using map()

rq.poly.4.NRI.netDiv_CWM_std_tw.mapped <- map(q10, ~rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                                       tau = .x,
                                                       method = "br",
                                                       data = prebinning_test)
)

coef.formatted.table.rq.poly.4.NRI <- stargazer::stargazer(rq.poly.4.NRI.netDiv_CWM_std_tw.mapped, 
                                                           rq.se = "boot",
                                                           column.labels = c(paste("tau = ", q10)),
                                                           covariate.labels = NULL,
                                                           dep.var.labels = "NRI",
                                                           omit = c("Constant"),
                                                           model.numbers = TRUE,
                                                           model.names =  FALSE,
                                                           keep.stat = c('n'),
                                                           type ='html')

write(coef.formatted.table.rq.poly.4.NRI, "coef.formatted.table.rq.poly.4.NRI.html")


pandoc_convert("coef.formatted.table.rq.poly.4.NRI.tex", to = "markdown")

# Overlay quantile model estimates on the data ------

summary.boot.rq.poly.4.NRI.netDiv_CWM_std_tw <- summary(rq.poly.4.NRI.netDiv_CWM_std_tw,
                                                        se = "boot",
                                                        R = 100)

# Performing the goodness-of-fit test -----
# ?WRTDStidal::goodfit()

goodfit.rq.poly.4.NRI <- vector(mode = "numeric", length(q10))

names(goodfit.rq.poly.4.NRI) <- q10

for(i in 1:length(q10)){
  ## conditional quantile model
  rq.poly.4.NRI.netDiv_CWM_std_tw.i <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                          tau = q10[i],
                                          #  method = "conquer",
                                          #  ci = TRUE,
                                          data = prebinning_test)
  
  resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i <- resid(rq.poly.4.NRI.netDiv_CWM_std_tw.i)
  
  ## non-conditional quantile model
  null.rq.poly.4.NRI.netDiv_CWM_std_tw.i <- rq(NRI ~ 1,
                                               tau = q10[i],
                                               # method = "conquer",
                                               # ci = TRUE,
                                               data = prebinning_test)
  
  null.resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i <- resid(null.rq.poly.4.NRI.netDiv_CWM_std_tw.i)
  
  ## Quantile regression goodness of fit
  
  goodfit.rq.poly.4.NRI.i <- goodfit(resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i, 
                                     null.resid.rq.poly.4.NRI.netDiv_CWM_std_tw.i, 
                                     q10[i])
  
  (goodfit.rq.poly.4.NRI[i] <- goodfit.rq.poly.4.NRI.i)
}

goodfit.rq.poly.4.NRI

kable(goodfit.rq.poly.4.NRI, format = "markdown")

# Qtools::GOFTest does not work

# GOFTest.rq.poly.4.NRI.netDiv_CWM_std_tw <- GOFTest(tau.rq.poly.4.NRI.netDiv_CWM_std_tw,
#                                                   B = 10, 
#                                                   seed = 12312) # does not work
# 

dev.off()

plot.summary.boot.rq.poly.4.NRI.netDiv_CWM_std_tw <- plot(summary.boot.rq.poly.4.NRI.netDiv_CWM_std_tw #,
                                                          # main = c("Intercept", 
                                                          # expression("Net Diversification Rate"[CWM[STD[tw]]]))
)

plot(rq.poly.4.NRI.netDiv_CWM_std_tw)

###

stargazer::stargazer(rq.poly.4.NRI.netDiv_CWM_std_tw[[1]], 
                     rq.se = "boot",
                     column.labels = c(paste("tau = ", q10)),
                     covariate.labels = c("In-situ diversification rates"),
                     dep.var.labels = "NRI",
                     omit = c("Constant"),
                     model.numbers = TRUE,
                     model.names =  TRUE,
                     keep.stat = c('n'),
                     type ='text')

visreg::visreg(rq.poly.4.NRI.netDiv_CWM_std_tw, 
               "netDiv_CWM_std_tw", 
               overlay = TRUE, 
               collapse = TRUE)

# Using geom_quantile()

ggplot(data = CWM.Div.MPD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling")%>%
         drop_na(ID_Realm), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
   geom_point(alpha = 0.25, 
              size = 1.4,
              stroke = 0.01,
              aes(colour = factor(ID_Realm))) +
  # scale_colour_viridis(option="D", 
  #                     begin = 0,
  #                     end = 1,
  #                     direction = 1, 
  #                     discrete = TRUE,
  #                     name = "Realm") +
  stat_quantile(quantiles = q10, 
                formula = y ~ poly(x, 4), 
                colour = "black") +
  #  xlim(min(prebinning_test$netDiv_CWM_std_tw),
  #       0.35) +
  # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


ggsave(filename = "figures/Fig.X.rq.poly.4.NRI.netDiv_CWM_std_tw_clean.png", 
       dpi = 300, 
       width = 8.5 , height = 8.5, 
       units = "in")


rq.poly.4.NRI.netDiv_CWM_std_tw.broom <- rq.poly.4.NRI.netDiv_CWM_std_tw %>% 
  broom::tidy(se.type = "boot",
              R = 100)

rq.poly.4.NRI.netDiv_CWM_std_tw.broom.f <- rq.poly.4.NRI.netDiv_CWM_std_tw.broom %>%
  filter(!grepl("factor", term)) %>% 
  mutate(term = as.factor(term))

rq.poly.4.NRI.netDiv_CWM_std_tw.broom.f$term <- recode(rq.poly.4.NRI.netDiv_CWM_std_tw.broom.f$term, 
                                                       "poly(netDiv_CWM_std_tw, 4)1" = "Linear",
                                                       "poly(netDiv_CWM_std_tw, 4)2" = "Quadratic",
                                                       "poly(netDiv_CWM_std_tw, 4)3" = "Cubic",
                                                       "poly(netDiv_CWM_std_tw, 4)4" = "Quartic")

rq.poly.4.NRI.netDiv_CWM_std_tw.broom.f %>%
  ggplot(aes(x = tau, 
             y = estimate))+
  geom_point(color = "#27408b", 
             size = 1)+ 
  geom_line(color = "#27408b", 
            size = 0.2) + 
  scale_x_continuous(breaks =  seq(0.05, 0.95, by = 0.15)) +
  # geom_smooth(method=  "lm", colour = "red", se = T) +  
  facet_wrap(~ term, 
             scales = "free", 
             ncol = 5) + 
  geom_ribbon(aes(
    ymin = estimate - std.error, 
    ymax = estimate + std.error),
    alpha = 0.25, 
    fill = "#27408b") +
  labs(y = "Estimate",
       x = c(expression("Quantile"~(tau)))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal",
        strip.background = element_rect(colour = NA,
                                        fill = "white")) + 
  guides(colour = guide_legend(nrow = 1))

# per-Realm models

# Neotropical

rq.poly.4.NRI.netDiv_CWM_std_tw.NT <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Neotropical"))

rq.poly.3.NRI.netDiv_CWM_std_tw.NT <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Neotropical"))

rq.poly.2.NRI.netDiv_CWM_std_tw.NT <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Neotropical"))

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw.NT, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw.NT$fitted.values)))
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw.NT, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw.NT$fitted.values)))
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw.NT, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw.NT$fitted.values)))

plot(rq.poly.4.NRI.netDiv_CWM_std_tw.NT)

# Nearctic

rq.poly.4.NRI.netDiv_CWM_std_tw.NA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Nearctic"))

rq.poly.3.NRI.netDiv_CWM_std_tw.NA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Nearctic"))

rq.poly.2.NRI.netDiv_CWM_std_tw.NA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Nearctic"))

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw.NA, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw.NA$fitted.values))) %>%
  mean()
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw.NA, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw.NA$fitted.values))) %>%
  mean()
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw.NA, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw.NA$fitted.values))) %>%
  mean()

plot(rq.poly.4.NRI.netDiv_CWM_std_tw.NA)

# Palearctic

rq.poly.4.NRI.netDiv_CWM_std_tw.PA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Palearctic"))

rq.poly.3.NRI.netDiv_CWM_std_tw.PA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Palearctic"))

rq.poly.2.NRI.netDiv_CWM_std_tw.PA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Palearctic"))

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw.PA, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw.PA$fitted.values)))  %>%
  mean()
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw.PA, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw.PA$fitted.values)))  %>%
  mean()
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw.PA, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw.PA$fitted.values))) %>%
  mean()

plot(rq.poly.4.NRI.netDiv_CWM_std_tw.PA)


# Afrotropical

rq.poly.4.NRI.netDiv_CWM_std_tw.AT <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Afrotropical"))

rq.poly.3.NRI.netDiv_CWM_std_tw.AT <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Afrotropical"))

rq.poly.2.NRI.netDiv_CWM_std_tw.AT <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Afrotropical"))

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw.AT, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw.AT$fitted.values)))  %>%
  mean()
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw.AT, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw.AT$fitted.values))) %>%
  mean()
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw.AT, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw.AT$fitted.values))) %>%
  mean()

plot(rq.poly.4.NRI.netDiv_CWM_std_tw.AT)

# Indomalay

rq.poly.4.NRI.netDiv_CWM_std_tw.IM <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Indomalay"))

rq.poly.3.NRI.netDiv_CWM_std_tw.IM <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Indomalay"))

rq.poly.2.NRI.netDiv_CWM_std_tw.IM <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Indomalay"))

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw.IM, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw.IM$fitted.values))) %>% 
  mean()
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw.IM, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw.IM$fitted.values))) %>% 
  mean()
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw.IM, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw.IM$fitted.values))) %>% 
  mean()

plot(rq.poly.4.NRI.netDiv_CWM_std_tw.IM)

# Australasian

rq.poly.4.NRI.netDiv_CWM_std_tw.AA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 4), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Australasian"))

rq.poly.3.NRI.netDiv_CWM_std_tw.AA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 3), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Australasian"))

rq.poly.2.NRI.netDiv_CWM_std_tw.AA <- rq(NRI ~ poly(netDiv_CWM_std_tw, 2), 
                                         tau = q10,
                                         method = "br",
                                         # ci = TRUE,
                                         data = prebinning_test %>%
                                           filter(ID_Realm == "Australasian"))

AIC(rq.poly.2.NRI.netDiv_CWM_std_tw.AA, 
    k = log(nrow(rq.poly.2.NRI.netDiv_CWM_std_tw.AA$fitted.values))) %>% 
  mean()
AIC(rq.poly.3.NRI.netDiv_CWM_std_tw.AA, 
    k = log(nrow(rq.poly.3.NRI.netDiv_CWM_std_tw.AA$fitted.values))) %>% 
  mean()
AIC(rq.poly.4.NRI.netDiv_CWM_std_tw.AA, 
    k = log(nrow(rq.poly.4.NRI.netDiv_CWM_std_tw.AA$fitted.values))) %>% 
  mean()

plot(rq.poly.4.NRI.netDiv_CWM_std_tw.AA)


# Using geom_quantile()

ggplot(data = CWM.Div.MPD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling")%>%
         drop_na(ID_Realm), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  geom_point(alpha = 0.25,
             size = 1.4,
             stroke = 0.01,
             aes(colour = factor(ID_Realm))) +
  scale_colour_viridis(option = "D",
                       begin = 0,
                       end = 1,
                       direction = 1,
                       discrete = TRUE,
                       name = "Realm") +
  stat_quantile(quantiles = q10, 
                formula = y ~ poly(x, 4), 
                colour = "black",
                method = "rq") +
  # ylim(min(prebinning_test$NRI) - 10,
  #      max(prebinning_test$NRI) + 5) +
  facet_grid(. ~ ID_Realm,
             scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", 
                                    size = 9),
        legend.box = "horizontal", 
        strip.background = element_rect(colour = NA,
                                        fill = "white")) + 
  guides(colour = guide_legend(nrow = 1))


# Intercept-only model (not used)

rq.poly.4.NRI.1 <- rq(NRI ~ 1, 
                      tau = q10,
                      # method = "conquer",
                      data = prebinning_test)

plot(rq.poly.4.NRI.1)


# Quantile GAM framework ----------------------------------------

qus <- seq(0.05, 0.95, 
           length.out = 10)

fit <- mqgam(NRI ~ s(netDiv_CWM_std_tw,  
                     k = 20, 
                     bs = "tp"), 
             data = prebinning_test, 
             qu = q10)

fitCycleM <- mgcViz::getViz(fit)
plot(fitCycleM)

#

plot(NRI ~ netDiv_CWM_std_tw, 
     data = prebinning_test, 
     col = "grey", 
     ylab = "y")

lines(x, f + qgamma(0.95, 4, 1), lty = 2)
lines(x, f + qgamma(0.5, 4, 1), lty = 2)
lines(x, f + qgamma(0.05, 4, 1), lty = 2)

lines(prebinning_test$netDiv_CWM_std_tw, qdo(fit, qus[1], predict), col = 2)
lines(prebinning_test$netDiv_CWM_std_tw, qdo(fit, qus[2], predict), col = 2)
lines(prebinning_test$netDiv_CWM_std_tw, qdo(fit, qus[3], predict), col = 2)

#

form <- mqgam(NRI ~ s(netDiv_CWM_std_tw,  
                      k = 10, 
                      bs = "ad"), 
              data = prebinning_test, 
              qu = qus)

set.seed(41241)

sigSeq <- seq(4, 8, length.out = 16)

closs <- tuneLearnFast(form = NRI ~ s(netDiv_CWM_std_tw,  
                                      k = 10, 
                                      bs = "ad"), 
                       data = prebinning_test, 
                       # lsig = sigSeq,  control = list("K" = 10), 
                       qu = qus,
                       multicore = TRUE, 
                       ncores = 10)

cqcheckI(obj = fit, 
         v = c("netDiv_CWM_std_tw"),
         X = prebinning_test,
         y = prebinning_test$NRI)

# Fit using estimated sigmas

fit <- mqgam(NRI ~ s(netDiv_CWM_std_tw,  k = 10, bs = "ad"), 
             data = prebinning_test, 
             qu = qus, 
             lsig = closs$lsig)

n <- nlrq(NRI ~ netDiv_CWM_std_tw, 
          data = prebinning_test, 
          start = list(k = 1, e = 2))

tidy(n)
augment(n)
glance(n)

library(ggplot2)
ggplot(augment(n), aes(wt, mpg)) +
  geom_point() +
  geom_line(aes(y = .fitted))

newdata <- head(mtcars)
newdata$wt <- newdata$wt + 1
augment(n, newdata = newdata)



# Plot fitted quantiles

plot(NRI ~ netDiv_CWM_std_tw, data = as.data.frame(prebinning_test), col = "grey", ylab = "NRI")

for(iq in qus){
  pred <- qdo(fit, iq, predict)
  lines(as.vector(prebinning_test$netDiv_CWM_std_tw), 
        pred, col = 2)
}  

fitCycleM <- mgcViz::getViz(form)
plot(fitCycleM)

single.fit <- qgam(NRI ~ s(netDiv_CWM_std_tw,  k = 20, bs = "ad") + ID_Realm, 
                   data = prebinning_test, 
                   qu = qus[5])

single.fit.getViz <- getViz(single.fit)

check(single.fit.getViz)

cqcheck(obj = single.fit, 
        v = c("netDiv_CWM_std_tw"), 
        X = prebinning_test, 
        y = prebinning_test$NRI)

print(plot(single.fit.getViz, 
           allTerms = TRUE, 
           #     select = c(1:3, 5)
), 
pages = 1)


xseq <- with(prebinning_test, 
             seq(min(netDiv_CWM_std_tw), max(netDiv_CWM_std_tw), 
                 length.out = 100))
preds <- sapply(fitCycleM, predict, newdata = data.frame(netDiv_CWM_std_tw = xseq))
plot(prebinning_test[,3:4], 
     ylim = range(preds),
     col = "grey", 
     xlab = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))
)

for(ii in 1:5) {
  lines(xseq, preds[ , ii], 
        col = 2)
}

tmp <- qdo(fit, iq)
aaa <- plot(tmp)

plot(NRI ~ netDiv_CWM_std_tw, data = prebinning_test, col = "grey", ylab = "y")

for(iq in q10){ 
  new.pred.data <- data.frame(pred.val = qdo(fit, iq, predict), 
                              netDiv_CWM_std_tw = prebinning_test$netDiv_CWM_std_tw,
                              row.names = prebinning_test$ID)
  lines(pred.val ~ netDiv_CWM_std_tw, 
        data = new.pred.data)
}

legend("top", c("truth", "fitted"), col = 2:1, lty = rep(1, 2))

library(splines)

plot(NRI ~ netDiv_CWM_std_tw, data = prebinning_test,
     # type="n", 
     cex = .75,
     col = "grey", ylab = "y")

X <- model.matrix(NRI ~ bs(netDiv_CWM_std_tw, df = 15), data= prebinning_test)

for(tau in 1:3/4){
  fit <- rq(NRI ~ bs(netDiv_CWM_std_tw, df = 15), tau = tau, data= prebinning_test)
  accel.fit <- X %*% fit$coef
  lines(accel.fit ~ netDiv_CWM_std_tw, data = prebinning_test)
}



###################################################
###################################################
###################################################

q10 <- seq(0.05, 0.95, 
           by = 0.1)

qs <- c(0.025, 0.25, 0.50, 0.75, 0.975)

# Using geom_quantile()

ggplot(data = CWM.Div.MPD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling") %>%
         drop_na(ID_Realm), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  geom_quantile(quantiles = q10,
                method = "rqss",
                colour="red") +
  facet_grid(. ~ ID_Realm) +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

#

ggsave(filename = "Fig.X.rq.NRI.netDiv_CWM_std_tw.png", 
       dpi = 300, 
       width = 8.5 , height = 8.5, 
       units = "in")


#

ggplot(data = CWM.Div.MPD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete = TRUE,
                      name="Realm") +
  geom_quantile(method = "rqss",
                lambda = 0.90, 
                quantiles = q10, colour="red") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

###############################

plot(prebinning_test$netDiv_CWM_std_tw, prebinning_test$NRI, 
     xlab = c(expression("Net Diversification Rate"[CWM[STD[tw]]])), 
     ylab = "acceleration (in g)",
     col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1),
     pch = 16)

###############################

hs <- c(1, 2, 3)

for(i in hs){
  h = hs[i]
  fit <- lprq(prebinning_test$netDiv_CWM_std_tw,
              prebinning_test$NRI,
              h = h,
              tau= q10)
  lines(fit$xx,
        fit$fv,
        lty=i)
}

legend(50,-70,c("h=1","h=2","h=3","h=4"),lty=1:length(hs))

###############################



###############################
########### Binning ###########
###############################
# 

dummy_double <- data.frame(x = rnorm(n = 500, mean = 0, sd = 1),
                           y = rnorm(n = 500, mean = 0, sd = 1)*2)

rbin::rbin_equal_freq(dummy_double,
                      y,
                      x,
                      bins = 10)

library(rbin)

prebinning_test <- CWM.Div.MPD.Chiroptera.Comm %>%
  filter(SamplingPool == "Global sampling") %>%
  dplyr::select(c(ID, ID_Realm, netDiv_CWM_std_tw, NRI)) %>%
  drop_na() %>%
  as_tibble()

post_binning <- prebinning_test %>% 
  mutate(NRI_bins = cut(NRI, breaks = 10),
         NRI_bins_seq = cut(NRI, seq(min(prebinning_test$NRI), max(prebinning_test$NRI), 2)),
         netDiv_CWM_std_tw_bins = cut(netDiv_CWM_std_tw, breaks = 10),
         netDiv_CWM_std_tw_bins_seq = cut(netDiv_CWM_std_tw, seq(min(prebinning_test$netDiv_CWM_std_tw), max(prebinning_test$netDiv_CWM_std_tw), 0.12))
  )

ggplot(data = post_binning, 
       aes(y = netDiv_CWM_std_tw, 
           x = NRI_bins_seq)) +
  geom_boxplot(aes(fill = factor(ID_Realm)),
               position = position_dodge(preserve = "single")) + 
  #  geom_jitter(position=position_jitter(0.1)) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete=TRUE,
                      name="Realm") +
  labs(y = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       x = c(expression("Net Relatedness Index"))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))



ggplot(data = post_binning, 
       aes(y = NRI, 
           x = netDiv_CWM_std_tw_bins)) +
  geom_boxplot(aes(fill = factor(ID_Realm)),
               position = position_dodge(preserve = "single")) + 
  # geom_jitter(position=position_jitter(0.1)) +
  scale_fill_viridis(option="D", 
                     begin = 0,
                     end = 1,
                     direction = 1, 
                     discrete=TRUE,
                     name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]])),
       y = c(expression("Net Relatedness Index"))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


ggsave(filename = "Fig.X.binned.boxplot.perRealm.NRI.netDiv_CWM_std_tw.png", 
       dpi = 300, 
       width = 9.5 , height = 6.5, 
       units = "in")

# Arithmetic means ######

ggplot(data = Div.MPD.Chiroptera.Comm, 
       aes(x = arithmeticMeans.netDiv, y = NRI)) +
  geom_point(alpha = 0.7, size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  stat_smooth(method = "loess") +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete=TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[arithmeticMean]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


#################################
########### Visualize ###########
#################################

# We will create a map using sp and maptools, because an issue prevented me from doing it with "sf"
# This implies in an slower procedure, but it will work.

# Load the polygon grid
polygonGrid_world_50km = readOGR(dsn="data/grids/", 
                                 layer="world_grid_prev_50km")

# Define its projection, if it isn't defined
proj4string(polygonGrid_world_50km) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

# Attribute an "id" variable
polygonGrid_world_50km@data$id = 1:length(polygonGrid_world_50km)


# Make it "ggplot2-friendly" by converting it to data.frame and using broom::tidy()
polygonGrid_world_50km.points <- broom::tidy(polygonGrid_world_50km)

# Explore to understand the result of fortifying: 
# head(polygonGrid_world_50km.points)

CWM.Div.MPD.Chiroptera.Comm$QuadratID <-  (CWM.Div.MPD.Chiroptera.Comm$ID - 1)

# Join CRI data to the fortified data
polygonGrid_world_50km.df <- polygonGrid_world_50km.points %>% 
  mutate_each(list(as.numeric), id) %>%
  left_join(CWM.Div.MPD.Chiroptera.Comm %>%
              filter(SamplingPool == "Global sampling"), 
            by = c("id" = "QuadratID"))

head(polygonGrid_world_50km.df)
head(CWM.Div.MPD.Chiroptera.Comm)

# Inspect: head(polygonGrid_world_50km.df)

#### Straightfoward visualisation ####

#pal <- wes_palette("Zissou1", 100, type = "continuous")
#pal <- wes_palette("Moonrise3", 100, type = "continuous")

raw.colPalette <- brewer.pal(9, 'YlOrBr')
pie(rep(1,9), col = raw.colPalette)

newcol <- colorRampPalette(raw.colPalette)
ncols <- 100

extended.colPalette <- newcol(ncols) #apply the function to get 100 colours

pie(rep(1, ncols), 
    col = extended.colPalette, 
    border = NA, 
    labels = NA)

netDiv_CWM_std_w.map <- ggplot() +
  # municipality polygons
  geom_polygon(data = polygonGrid_world_50km.df, aes(fill = netDiv_CWM_std_w, 
                                                     x = long, 
                                                     y = lat, 
                                                     group = group)) +
  # no outline
  #geom_path(data = map_data, aes(x = long, 
  #                               y = lat, 
  #                               group = group), 
  #          color = "white", size = 0.1) +
  
  coord_equal() +
  
  # add a basic theme
  theme_map() +
  
  theme(legend.position = "bottom") +
  
  scale_fill_gradientn(colours = extended.colPalette,
                       name = "Net Diversification Rates", 
                       guide = guide_colorbar(
                         direction = "horizontal",
                         barheight = unit(2, units = "mm"),
                         barwidth = unit(75, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         # some shifting around
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       )) + 
  # insert a North arrow
  ggsn::north(polygonGrid_world_50km.df, 
              symbol = 12, scale = 0.04,
              location = "bottomright",
              anchor = c(x = +14500000, y = -6850000)) +
  
  labs(x = NULL, 
       y = NULL, 
       title = "Net Diversification Rates", 
       subtitle = "Community Weighted Mean",
       caption = "Bats")

netDiv_CWM_std_w.map




# Looking a normal linear model

lm.NRI.netDiv_CWM_std_w <- lm(NRI ~ netDiv_CWM_std_w, 
                              data = CWM.Div.MPD.Chiroptera.Comm %>%
                                filter(SamplingPool == "Global sampling"))


summary(lm.NRI.netDiv_CWM_std_w)

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(lm.NRI.netDiv_CWM_std_w, which=1)
plot(lm.NRI.netDiv_CWM_std_w, which=2)
plot(lm.NRI.netDiv_CWM_std_w, which=3)
plot(lm.NRI.netDiv_CWM_std_w, which=5)

dev.off()

# Diagnostic plots show a clear violation of linear assumptions

# Looking a normal linear model

quad.lm.NRI.netDiv_CWM_std_w <- lm(NRI ~ netDiv_CWM_std_w + I(netDiv_CWM_std_w^2), 
                                   data = CWM.Div.MPD.Chiroptera.Comm %>%
                                     filter(SamplingPool == "Global sampling"))

summary(quad.lm.NRI.netDiv_CWM_std_w)

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(quad.lm.NRI.netDiv_CWM_std_w, which=1)
plot(quad.lm.NRI.netDiv_CWM_std_w, which=2)
plot(quad.lm.NRI.netDiv_CWM_std_w, which=3)
plot(quad.lm.NRI.netDiv_CWM_std_w, which=5)

dev.off()

# Diagnostic plots show a clear violation of linear assumptions

# Stepping up to a double-generalized linear model appraoch that consists of
# three parts: a model for the variance, a GLM for the mean and a GLM for the
# disperson

dglm.NRI.netDiv<- dglm::dglm(NRI ~ netDiv_CWM_std_w + I(netDiv_CWM_std_w^2), 
                             data = CWM.Div.MPD.Chiroptera.Comm %>%
                               filter(SamplingPool == "Global sampling"),
                             dformula = ~netDiv_CWM_std_w,
                             # dlink = "identity",
                             method = "reml")
# Summarize the mean model as for a glm
summary.glm(dglm.NRI.netDiv)

# Summarize the dispersion model as for a glm
summary(dglm.NRI.netDiv$dispersion.fit)

# Examine goodness of fit of dispersion model by plotting residuals
plot(fitted(dglm.NRI.netDiv$dispersion.fit),
     residuals(dglm.NRI.netDiv$dispersion.fit))


dglm.Pvalues <- function(dglm.fit){
  P.disp = anova.dglm(dglm.fit)$Adj.P[2]
  P.mean = summary(dglm.fit)$coef[2,4]
  list(P.mean=P.mean, P.disp=P.disp)
}

dglm.Pvalues(dglm.NRI.netDiv) # Not working

anova(dglm.NRI.netDiv)

summary.glm(dglm.NRI.netDiv)

summary(dglm.NRI.netDiv$dispersion.fit)

#

ggplot(data = CWM.Div.MPD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_w, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  stat_smooth(method = "gam",
              formula = y ~ poly(x, 2), 
              se = TRUE,
              colour = "black") +
  #  stat_smooth(method = "lm", formula= (y ~ poly(x, 2)), se = TRUE, colour = "black") +
  #  stat_smooth(method = "loess", 
  #              formula = y ~ poly(x, 2), size = 1) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete=TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[w]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


## Attempting with a locally estimated scatterplot smoothing (loess)

loess.NRI.netDiv_CWM_std_w <- loess(NRI ~ netDiv_CWM_std_w, 
                                    data = CWM.Div.MPD.Chiroptera.Comm %>%
                                      filter(SamplingPool == "Global sampling"), 
                                    span = 0.50)

loess.NRI.netDiv_CWM_std_w.fit <- augment(loess.NRI.netDiv_CWM_std_w)

# get smoothed output
smoothed.loess.NRI.netDiv_CWM_std_w <- predict(loess.NRI.netDiv_CWM_std_w) 

ggplot(loess.NRI.netDiv_CWM_std_w.fit, 
       aes(netDiv_CWM_std_w, NRI)) +
  geom_point() +
  geom_line(aes(y = .fitted), 
            color = "red")

ggplot(data = loess.NRI.netDiv_CWM_std_w.fit, 
       aes(x = netDiv_CWM_std_w, 
           y = NRI)) +
  geom_point(alpha = 0.7, 
             size = 1.4) +
  geom_line(aes(y = .fitted), 
            color = "red") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[w]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))

# 

ggplot(data = Div.MPD.Chiroptera.Comm, 
       aes(x = arithmeticMeans.netDiv, y = NRI)) +
  geom_point(alpha = 0.7, size = 1.4, 
             aes(colour = factor(ID_Realm))) +
  stat_smooth(method = "loess") +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete=TRUE,
                      name="Realm") +
  labs(x = c(expression("Net Diversification Rate"[CWM[arithmeticMean]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


#################################
########### Visualize ###########
#################################

# We will create a map using sp and maptools, because an issue prevented me from doing it with "sf"
# This implies in an slower procedure, but it will work.

# Load the polygon grid
polygonGrid_world_50km = readOGR(dsn="data/grids/", 
                                 layer="world_grid_prev_50km")

# Define its projection, if it isn't defined
proj4string(polygonGrid_world_50km) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

# Attribute an "id" variable
polygonGrid_world_50km@data$id = 1:length(polygonGrid_world_50km)


##

polygonGrid_world_50km.sf <- st_as_sf(polygonGrid_world_50km)

# we can fix this using the lwgeom library...
polygonGrid_world_50km.sf <- st_make_valid(polygonGrid_world_50km.sf)

# Make it "ggplot2-friendly" by converting it to data.frame and using broom::tidy()
polygonGrid_world_50km.points <- broom::tidy(polygonGrid_world_50km)

# Explore to understand the result of fortifying: 
# head(polygonGrid_world_50km.points)

CWM.Div.MPD.Chiroptera.Comm$QuadratID <-  (CWM.Div.MPD.Chiroptera.Comm$ID - 1)

# Join CRI data to the fortified data
polygonGrid_world_50km.df <- polygonGrid_world_50km.points %>% 
  mutate_each(list(as.numeric), id) %>%
  left_join(CWM.Div.MPD.Chiroptera.Comm %>%
              filter(SamplingPool == "Global sampling"), 
            by = c("id" = "QuadratID"))

head(polygonGrid_world_50km.df)
head(CWM.Div.MPD.Chiroptera.Comm)

# Inspect: head(polygonGrid_world_50km.df)

#### Straightfoward visualisation ####

#pal <- wes_palette("Zissou1", 100, type = "continuous")
#pal <- wes_palette("Moonrise3", 100, type = "continuous")

raw.colPalette <- brewer.pal(9, 'YlOrBr')
pie(rep(1,9), col = raw.colPalette)

newcol <- colorRampPalette(raw.colPalette)
ncols <- 100

extended.colPalette <- newcol(ncols) #apply the function to get 100 colours

pie(rep(1, ncols), 
    col = extended.colPalette, 
    border = NA, 
    labels = NA)

p <- ggplot() +
  # municipality polygons
  geom_polygon(data = polygonGrid_world_50km.df, aes(fill = netDiv_CWM_std_w, 
                                                     x = long, 
                                                     y = lat, 
                                                     group = group)) +
  # no outline
  #geom_path(data = map_data, aes(x = long, 
  #                               y = lat, 
  #                               group = group), 
  #          color = "white", size = 0.1) +
  
  coord_equal() +
  
  # add a basic theme
  theme_map() +
  
  theme(legend.position = "bottom") +
  
  scale_fill_gradientn(colours = extended.colPalette,
                       name = "Net Diversification Rates", 
                       guide = guide_colorbar(
                         direction = "horizontal",
                         barheight = unit(2, units = "mm"),
                         barwidth = unit(75, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         # some shifting around
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       )) + 
  # insert a North arrow
  ggsn::north(polygonGrid_world_50km.df, 
              symbol = 12, scale = 0.04,
              location = "bottomright",
              anchor = c(x = +14500000, y = -6850000)) +
  
  labs(x = NULL, 
       y = NULL, 
       title = "Net Diversification Rates", 
       subtitle = "Community Weighted Mean",
       caption = "Bats")
p


# NTI

rq.NTI.netDiv_CWM_std_tw <- rq(NTI ~ netDiv_CWM_std_tw, 
                               tau = q10,
                               method = "br",
                               data = prebinning_test)

rq.NTI.netDiv_CWM_std_tw

summary(rq.NTI.netDiv_CWM_std_tw)

# Represent 

plot(rq.NTI.netDiv_CWM_std_tw)

# Plot distributions

plot(summary(rq.NTI.netDiv_CWM_std_tw,
             se = "nid"),
     main = c("Intercept", 
              expression("Net Diversification Rate"[CWM[STD[tw]]])))

## Polynomial quantile regression fit ----------------------

# Worldwide model

rq.poly.5.NTI.netDiv_CWM_std_tw <- rq(NTI ~ poly(netDiv_CWM_std_tw, 5), 
                                      tau = q20,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

rq.poly.4.NTI.netDiv_CWM_std_tw <- rq(NTI ~ poly(netDiv_CWM_std_tw, 4), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

rq.poly.3.NTI.netDiv_CWM_std_tw <- rq(NTI ~ poly(netDiv_CWM_std_tw, 3), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

rq.poly.2.NTI.netDiv_CWM_std_tw <- rq(NTI ~ poly(netDiv_CWM_std_tw, 2), 
                                      tau = q10,
                                      method = "br",
                                      # ci = TRUE,
                                      data = prebinning_test)

AIC(rq.poly.2.NTI.netDiv_CWM_std_tw, k = log(nrow(rq.poly.2.NTI.netDiv_CWM_std_tw$fitted.values)))
AIC(rq.poly.3.NTI.netDiv_CWM_std_tw, k = log(nrow(rq.poly.3.NTI.netDiv_CWM_std_tw$fitted.values)))
AIC(rq.poly.4.NTI.netDiv_CWM_std_tw, k = log(nrow(rq.poly.4.NTI.netDiv_CWM_std_tw$fitted.values)))

AIC(rq.poly.2.NTI.netDiv_CWM_std_tw)
AIC(rq.poly.3.NTI.netDiv_CWM_std_tw)
AIC(rq.poly.4.NTI.netDiv_CWM_std_tw)


anova(rq.poly.3.NTI.netDiv_CWM_std_tw, joint = FALSE)

anova(rq.poly.4.NTI.netDiv_CWM_std_tw,
      rq.poly.3.NTI.netDiv_CWM_std_tw,
      rq.poly.2.NTI.netDiv_CWM_std_tw[[1]], joint = FALSE)

## Obtaining the formatted table of coefficients for the preferred model --------

# Using map()

rq.poly.4.NTI.netDiv_CWM_std_tw.mapped <- map(q10, ~rq(NTI ~ poly(netDiv_CWM_std_tw, 4), 
                                                       tau = .x,
                                                       method = "br",
                                                       data = prebinning_test)
)

coef.formatted.table.rq.poly.4.NTI <- stargazer::stargazer(rq.poly.4.NTI.netDiv_CWM_std_tw.mapped, 
                                                           rq.se = "boot",
                                                           column.labels = c(paste("tau = ", q10)),
                                                           covariate.labels = NULL,
                                                           dep.var.labels = "NTI",
                                                           omit = c("Constant"),
                                                           model.numbers = TRUE,
                                                           model.names =  FALSE,
                                                           keep.stat = c('n'),
                                                           type ='html')

write(coef.formatted.table.rq.poly.4.NTI, "coef.formatted.table.rq.poly.4.NTI.html")


pandoc_convert("coef.formatted.table.rq.poly.4.NTI.tex", to = "markdown")

# Overlay quantile model estimates on the data ------

summary.boot.rq.poly.4.NTI.netDiv_CWM_std_tw <- summary(rq.poly.4.NTI.netDiv_CWM_std_tw,
                                                        se = "boot",
                                                        R = 100)

# Performing the goodness-of-fit test -----
# ?WRTDStidal::goodfit()

goodfit.rq.poly.4.NTI <- vector(mode = "numeric", length(q10))

names(goodfit.rq.poly.4.NTI) <- q10

for(i in 1:length(q10)){
  ## conditional quantile model
  rq.poly.4.NTI.netDiv_CWM_std_tw.i <- rq(NTI ~ poly(netDiv_CWM_std_tw, 4), 
                                          tau = q10[i],
                                          #  method = "conquer",
                                          #  ci = TRUE,
                                          data = prebinning_test)
  
  resid.rq.poly.4.NTI.netDiv_CWM_std_tw.i <- resid(rq.poly.4.NTI.netDiv_CWM_std_tw.i)
  
  ## non-conditional quantile model
  null.rq.poly.4.NTI.netDiv_CWM_std_tw.i <- rq(NTI ~ 1,
                                               tau = q10[i],
                                               # method = "conquer",
                                               # ci = TRUE,
                                               data = prebinning_test)
  
  null.resid.rq.poly.4.NTI.netDiv_CWM_std_tw.i <- resid(null.rq.poly.4.NTI.netDiv_CWM_std_tw.i)
  
  ## Quantile regression goodness of fit
  
  goodfit.rq.poly.4.NTI.i <- goodfit(resid.rq.poly.4.NTI.netDiv_CWM_std_tw.i, 
                                     null.resid.rq.poly.4.NTI.netDiv_CWM_std_tw.i, 
                                     q10[i])
  
  (goodfit.rq.poly.4.NTI[i] <- goodfit.rq.poly.4.NTI.i)
}

goodfit.rq.poly.4.NTI

kable(goodfit.rq.poly.4.NTI, format = "markdown")

# Qtools::GOFTest does not work

# GOFTest.rq.poly.4.NTI.netDiv_CWM_std_tw <- GOFTest(tau.rq.poly.4.NTI.netDiv_CWM_std_tw,
#                                                   B = 10, 
#                                                   seed = 12312) # does not work
# 

dev.off()

plot.summary.boot.rq.poly.4.NTI.netDiv_CWM_std_tw <- plot(summary.boot.rq.poly.4.NTI.netDiv_CWM_std_tw #,
                                                          # main = c("Intercept", 
                                                          # expression("Net Diversification Rate"[CWM[STD[tw]]]))
)

plot(rq.poly.4.NTI.netDiv_CWM_std_tw)

###

stargazer::stargazer(rq.poly.4.NTI.netDiv_CWM_std_tw[[1]], 
                     rq.se = "boot",
                     column.labels = c(paste("tau = ", q10)),
                     covariate.labels = c("In-situ diversification rates"),
                     dep.var.labels = "NTI",
                     omit = c("Constant"),
                     model.numbers = TRUE,
                     model.names =  TRUE,
                     keep.stat = c('n'),
                     type ='text')



visreg::visreg(rq.poly.4.NTI.netDiv_CWM_std_tw, 
               "netDiv_CWM_std_tw", 
               overlay = TRUE, 
               collapse = TRUE)


# Using geom_quantile()

ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling")%>%
         drop_na(ID_Realm), 
       aes(x = netDiv_CWM_std_tw, 
           y = NTI)) +
  # geom_point(alpha = 0.25, 
  #            size = 1.4,
  #            stroke = 0.01,
  #            aes(colour = factor(ID_Realm))) +
  # scale_colour_viridis(option="D", 
  #                     begin = 0,
  #                     end = 1,
  #                     direction = 1, 
  #                     discrete = TRUE,
  #                     name = "Realm") +
  stat_quantile(quantiles = q10, 
                formula = y ~ poly(x, 4), 
                colour = "black"
  ) +
  xlim(min(prebinning_test$netDiv_CWM_std_tw),
       0.35) +
  # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))


ggsave(filename = "figures/Fig.X.rq.poly.4.NTI.netDiv_CWM_std_tw_clean.png", 
       dpi = 300, 
       width = 8.5 , height = 8.5, 
       units = "in")


rq.poly.4.NTI.netDiv_CWM_std_tw.broom <- rq.poly.4.NTI.netDiv_CWM_std_tw %>% 
  broom::tidy(se.type = "boot",
              R = 100)

rq.poly.4.NTI.netDiv_CWM_std_tw.broom.f <- rq.poly.4.NTI.netDiv_CWM_std_tw.broom %>%
  filter(!grepl("factor", term)) %>% 
  mutate(term = as.factor(term))

rq.poly.4.NTI.netDiv_CWM_std_tw.broom.f$term <- recode(rq.poly.4.NTI.netDiv_CWM_std_tw.broom.f$term, 
                                                       "poly(netDiv_CWM_std_tw, 4)1" = "Linear",
                                                       "poly(netDiv_CWM_std_tw, 4)2" = "Quadratic",
                                                       "poly(netDiv_CWM_std_tw, 4)3" = "Cubic",
                                                       "poly(netDiv_CWM_std_tw, 4)4" = "Quartic")

rq.poly.4.NTI.netDiv_CWM_std_tw.broom.f %>%
  ggplot(aes(x = tau, 
             y = estimate))+
  geom_point(color = "#27408b", 
             size = 1)+ 
  geom_line(color = "#27408b", 
            size = 0.2) + 
  scale_x_continuous(breaks =  seq(0.05, 0.95, 
                                   by = 0.15)) +
  # geom_smooth(method=  "lm", colour = "red", se = T) +  
  facet_wrap(~ term, 
             scales = "free", 
             ncol = 5) + 
  geom_ribbon(aes(
    ymin = estimate - std.error, 
    ymax = estimate + std.error),
    alpha = 0.25, 
    fill = "#27408b") +
  labs(y = "Estimate",
       x = c(expression("Quantile"~(tau)))) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal",
        strip.background = element_rect(colour = NA,
                                        fill = "white")) + 
  guides(colour = guide_legend(nrow = 1))
