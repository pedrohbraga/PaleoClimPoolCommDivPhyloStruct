#################################################################################
### Code to obtain local diversification rates for each community based on    ### 
### rates estimated for each tip in the phylogeny                             ###
#                                                                               #
# Current measures calculated:                                                  #
#
# 1. arithmeticMeans.lambda, arithmeticMeans.mu and arithmeticMeans.netdiv:     #
#                                                                               #
# For the arithmetic mean, the rate for the kth grid cell is computed as        
# λk = (Σwiλi)/Σwi in which Σ denotes a summation over all N species (i = 1 to
# i = N) present in the kth cell, and wi is the weight assigned to the ith
# species. We computed weights for each species as the inverse of the number of
# grid cells in which the species was found
# 
# 2. lambda.CWM.Chiroptera.Comm, mu.CWM.Chiroptera.Comm,
# netdiv.CWM.Chiroptera.Comm:
# 
# 
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2020-04-02"                                                     #
#                                                                               # 
#################################################################################

length(tip.rates$lambda.avg)
length(tip.rates$mu.avg)


### Computing arithmetic weighted means for diversification rates

# For the arithmetic mean, the rate for the kth grid cell is computed as
# λk = (Σwiλi)/Σwi in which Σ denotes a summation over all N species (i = 1 to
# i = N) present in the kth cell, and wi is the weight assigned to the ith
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


arithmeticMeans.lambda
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


arithmeticMeans.mu

colnames(arithmeticMeans.mu) <- "arithmeticMeans.mu"

arithmeticMeans.mu

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
CWM.Div.MPD.Chiroptera.Comm %>%
  filter(SamplingPool == "Global sampling") %>%
  filter(ntaxa > 2)

# CWM.Div.MPD.Chiroptera.Comm <- read.csv(data/matrices/CWM.Div.MPD.Chiroptera.Comm.csv", h = T)

##################################################################
### Representation and analyses for net diversification rates ####
##################################################################

# Looking a normal linear model

lm.NRI.netDiv_CWM_std_tw <- lm(NRI ~ netDiv_CWM_std_tw, 
                               data = CWM.Div.MPD.Chiroptera.Comm %>%
                                 filter(SamplingPool == "Global sampling"))


summary(lm.NRI.netDiv_CWM_std_tw)

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(lm.NRI.netDiv_CWM_std_tw, which=1)
plot(lm.NRI.netDiv_CWM_std_tw, which=2)
plot(lm.NRI.netDiv_CWM_std_tw, which=3)
plot(lm.NRI.netDiv_CWM_std_tw, which=5)

dev.off()

# Diagnostic plots show a clear violation of linear assumptions

# Looking a normal linear model

quad.lm.NRI.netDiv_CWM_std_tw <- lm(NRI ~ netDiv_CWM_std_tw + I(netDiv_CWM_std_tw^2), 
                                    data = CWM.Div.MPD.Chiroptera.Comm %>%
                                      filter(SamplingPool == "Global sampling"))

summary(quad.lm.NRI.netDiv_CWM_std_tw)

layout(matrix(1:4, 
              nrow = 2, 
              ncol = 2, 
              byrow = TRUE))

plot(quad.lm.NRI.netDiv_CWM_std_tw, which=1)
plot(quad.lm.NRI.netDiv_CWM_std_tw, which=2)
plot(quad.lm.NRI.netDiv_CWM_std_tw, which=3)
plot(quad.lm.NRI.netDiv_CWM_std_tw, which=5)

dev.off()

# Diagnostic plots show a clear violation of linear assumptions

# Stepping up to a double-generalized linear model appraoch that consists of
# three parts: a model for the variance, a GLM for the mean and a GLM for the
# disperson

dglm.NRI.netDiv<- dglm::dglm(NRI ~ netDiv_CWM_std_tw + I(netDiv_CWM_std_tw^2), 
                             data = CWM.Div.MPD.Chiroptera.Comm %>%
                               filter(SamplingPool == "Global sampling"),
                             dformula = ~netDiv_CWM_std_tw,
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

#####
#####
#####

ggplot(data = CWM.Div.MPD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling"), 
       aes(x = netDiv_CWM_std_tw, 
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


##


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

#####
#####
#####

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
