# Variation Partitioning Framework

################################
#### Variation partitioning ####
################################


###########################################
### PCA on the environmental variables ####
###########################################

# Create a final list that will contain a list for the results for each biome being tested 

varPartStore.all <- list()


(Biomes <- names(sort(t(as.data.frame(lapply(SubRegion.Biome.Comm.Phylo, 
                                             function(x) nrow(x$Comm))))[, 1]))) 

# Selecting the world's biomes

for(SelectedBiome in Biomes){
  
  Chiroptera.Biome.MPD <- MPD.LatLong.Env.AllScales %>%
    filter(ID %in% row.names(SubRegion.Biome.Comm.Phylo[[SelectedBiome]]$Comm)) %>%
    filter(SamplingPool == "Global sampling") %>%
    filter(!is.na(NRI)) %>%
    column_to_rownames('ID') %>%
    select(NRI)
  
  if(nrow(Chiroptera.Biome.MPD) >= 5){
    # Subset environmental datasets
    biome.current.Env <- worldClimate %>%
      filter(row.names(worldClimate) %in% row.names(Chiroptera.Biome.MPD)) %>%
      select(starts_with("cur_"))
    
    biome.mid.holo.Env <-   worldClimate %>%
      filter(row.names(worldClimate) %in% row.names(Chiroptera.Biome.MPD)) %>%
      select(starts_with("MH_"))
    
    biome.lgm.Env <-   worldClimate %>%
      filter(row.names(worldClimate) %in% row.names(Chiroptera.Biome.MPD)) %>%
      select(starts_with("LGM_"))
    
    # biome.legacy.current.lgm.Env <- SubRegion.Biome.Comm.Phylo[[SelectedBiome]]$legacy.current.lgm.Env
    
    # Perform PCA for the current, LGM, Mid-Holocene periods
    
    pca.biome.current.Env <- dudi.pca(biome.current.Env,
                                      row.w = rep(1, nrow(biome.current.Env))/nrow(biome.current.Env), # uniform row weights
                                      col.w = rep(1, ncol(biome.current.Env)), # unit column weights 
                                      center = TRUE, # center by the mean
                                      scale = TRUE, # norm the vectors for the row.w weighting
                                      scannf = FALSE, nf = 18)
    
    pca.biome.mid.holo.Env <- dudi.pca(biome.mid.holo.Env,
                                       row.w = rep(1, nrow(biome.mid.holo.Env))/nrow(biome.mid.holo.Env),
                                       col.w = rep(1, ncol(biome.mid.holo.Env)), 
                                       center = TRUE, 
                                       scale = TRUE, 
                                       scannf = FALSE, nf = 18)
    
    pca.biome.lgm.Env <- dudi.pca(biome.lgm.Env,
                                  row.w = rep(1, nrow(biome.lgm.Env))/nrow(biome.lgm.Env),
                                  col.w = rep(1, ncol(biome.lgm.Env)), 
                                  center = TRUE, 
                                  scale = TRUE, 
                                  scannf = FALSE, nf = 18)
    
    # Explore the proportion of thevariance associated to each axis obtained in the PCAs
    
    biome.explained.variation <- cbind(100 * pca.biome.current.Env$eig/sum(pca.biome.current.Env$eig),
                                       100 * pca.biome.mid.holo.Env$eig/sum(pca.biome.mid.holo.Env$eig),
                                       100 * pca.biome.lgm.Env$eig/sum(pca.biome.lgm.Env$eig))
    
    rownames(biome.explained.variation) <- c(paste0("Axis ", as.character(1:nrow(biome.explained.variation))))
    
    colnames(biome.explained.variation) <- c("Current",
                                             "Mid-Holocene",
                                             "LGM")
    
    # For further analyses, apply eigenvector selection criteria
    
    # See brokenStick.selection() function on the Supplementary Material
    
    nAxes.PCA.current.Env <- brokenStick.selection(pca.biome.current.Env)$brokenStick.kept
    nAxes.PCA.mid.holo.Env <- brokenStick.selection(pca.biome.mid.holo.Env)$brokenStick.kept
    nAxes.PCA.lgm.Env <- brokenStick.selection(pca.biome.lgm.Env)$brokenStick.kept
    
    # Selection of Principal Components
    
    biome.current.Env.PC <- as.matrix(pca.biome.current.Env$li[ , 1:nAxes.PCA.current.Env])
    biome.mid.holo.Env.PC <- as.matrix(pca.biome.mid.holo.Env$li[ , 1:nAxes.PCA.mid.holo.Env])
    biome.lgm.Env.PC <- as.matrix(pca.biome.lgm.Env$li[ , 1:nAxes.PCA.lgm.Env])
    
    colnames(biome.current.Env.PC) <- c(paste0("Current.", "PC", 1:nAxes.PCA.current.Env))
    colnames(biome.mid.holo.Env.PC) <- c(paste0("MidHolocene.", "PC", 1:nAxes.PCA.mid.holo.Env))
    colnames(biome.lgm.Env.PC) <- c(paste0("LGM.", "PC", 1:nAxes.PCA.lgm.Env))
    
    biome.Env.PC <- cbind(pca.biome.current.Env$li[ , 1:nAxes.PCA.current.Env],
                          pca.biome.mid.holo.Env$li[ , 1:nAxes.PCA.mid.holo.Env],
                          pca.biome.lgm.Env$li[ , 1:nAxes.PCA.lgm.Env])
    
    colnames(biome.Env.PC) <- c(paste0("Current.", "PC", 1:nAxes.PCA.current.Env),
                                paste0("MidHolocene.", "PC", 1:nAxes.PCA.mid.holo.Env),
                                paste0("LGM.", "PC", 1:nAxes.PCA.lgm.Env))
    
    
    ##################################### 
    ###### Moran's Eigenvector Maps #####
    #####################################
    
    Biome.ID.LatLong <- regions.LatLong %>%
      filter(row.names(regions.LatLong) %in% row.names(Chiroptera.Biome.MPD)) %>%
      column_to_rownames('ID') %>%
      select(c("Longitude", "Latitude"))
    
    Biome.SpatialDistances <- dist(Biome.ID.LatLong)
    
    Biome.MEM <- adespatial::dbmem(Biome.SpatialDistances, 
                                   thresh = NULL, # create a connectivity matrix based on a threshold distance and minimum spanning tree algorithm
                                   MEM.autocor = "positive", # keep all eigenvalues, except the null ones
                                   silent = FALSE)
    dev.off()
    barplot(attr(Biome.MEM, "values"), 
            main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)
    
    ################################
    ##### Variation paritioning ####
    ################################
    
    # Chiroptera.Biome.MPD
    
    # Redundancy analysis between the Hellinger-transformed community data and the first spatial eigenvalues
    Chiroptera.Biome.MPD.MEM.lm <- lm(as.matrix(Chiroptera.Biome.MPD) ~ as.matrix(Biome.MEM))
    
    # Obtain 
    (Chiroptera.Biome.MPD.MEM.lm.Rsq <- RsquareAdj(Chiroptera.Biome.MPD.MEM.lm))
    
    #if(Chiroptera.Biome.MPD.MEM.lm.Rsq$adj.r.squared > 0)
    # Forward selection of the MEM that best explain the species variation
    ForwardSel.MEM.Biome <- try(adespatial::forward.sel(as.matrix(Chiroptera.Biome.MPD), 
                                                        as.matrix(Biome.MEM),
                                                        adjR2thresh = Chiroptera.Biome.MPD.MEM.lm.Rsq$adj.r.squared,
                                                        nperm = 999))
    if(class(ForwardSel.MEM.Biome) == "try-error") {
      cat("Caught an error during the foward selection step. Keeping the original MEM. \n")
      Selected.MEM <- as.matrix(Biome.MEM)
    } else {
      Selected.MEM <- as.matrix(Biome.MEM[, sort(ForwardSel.MEM.Biome$order)])
    }
    
    # You may visualize the MEMs by: 
    ## require("ade4", quietly = TRUE) & require("adegraphics", quietly = TRUE)
    ## s.value(Biome.ID.LatLong, Selected.MEM)
    
    # Peform variation partitioning on the basis of redundancy analyses
    
    Biome.varpart <- varpart(Chiroptera.Biome.MPD, # Response variable (Y)
                             biome.current.Env.PC, # X1
                             biome.mid.holo.Env.PC,# X2
                             biome.lgm.Env.PC,   # X3
                             as.matrix(Selected.MEM)) # X4
    
    #############################
    ##### Testing fractions #####
    #############################
    
    # [a] Pure Current Climate
    aFrac <- rda(Chiroptera.Biome.MPD ~ biome.current.Env.PC +
                   Condition(biome.mid.holo.Env.PC) + Condition(biome.lgm.Env.PC) + Condition(Selected.MEM))
    
    # [c] Pure Mid-Holocene Climate
    cFrac <-rda(Chiroptera.Biome.MPD ~ biome.mid.holo.Env.PC +
                  Condition(biome.current.Env.PC) + Condition(biome.lgm.Env.PC) + Condition(Selected.MEM))
    
    # [b] Pure LGM Climate
    bFrac <- rda(Chiroptera.Biome.MPD ~ biome.lgm.Env.PC +
                   Condition(biome.current.Env.PC) + Condition(biome.mid.holo.Env.PC) + Condition(Selected.MEM))
    
    # [d] Pure Spatial Climate
    dFrac <- rda(Chiroptera.Biome.MPD ~ Selected.MEM  +
                   Condition(biome.current.Env.PC) + Condition(biome.mid.holo.Env.PC) + Condition(biome.lgm.Env.PC) )
    
    # [aebglfc] Current and Past Climate (Condition: Space)
    aebglfcFrac <- rda(Chiroptera.Biome.MPD ~ biome.current.Env.PC + biome.mid.holo.Env.PC + biome.lgm.Env.PC + Condition(Selected.MEM))
    
    # Save most results directly on Biome.varpart
    
    Biome.varpart$part$indfract$P.value <- rep("NA", nrow(Biome.varpart$part$indfract))
    Biome.varpart$part$fract$P.value <- rep("NA", nrow(Biome.varpart$part$fract))
    Biome.varpart$part$contr1$P.value <- rep("NA", nrow(Biome.varpart$part$contr1))
    Biome.varpart$part$contr2$P.value <- rep("NA", nrow(Biome.varpart$part$contr2))
    
    # Test individual fractions 
    Biome.varpart$part$indfract$P.value[1:4] <- c(anova(aFrac, permutations = 99)$`Pr(>F)`[1],
                                                  anova(bFrac, permutations = 99)$`Pr(>F)`[1],
                                                  anova(cFrac, permutations = 99)$`Pr(>F)`[1],
                                                  anova(dFrac, permutations = 99)$`Pr(>F)`[1])
    
    # Test Current and Past Climate (Condition: Space)
    aebglfc.FracTest <- cbind(t(as.matrix(RsquareAdj(aebglfcFrac))), 
                              as.matrix(anova(aebglfcFrac, permutations = 99))[1, 4, drop = FALSE])
    
    Biome.varpart
    
    ####################################
    ##### Store results of interest ####
    ####################################
    
    varPartStore <- list(SelectedBiome = SelectedBiome,
                         biome.Env.PC = biome.Env.PC,
                         biome.explained.variation = biome.explained.variation,
                         Biome.MEM = Biome.MEM,
                         Biome.varpart = Biome.varpart,
                         ClimateOnlyFrac = aebglfc.FracTest)
    
    
    varPartStore.all[[paste0(SelectedBiome)]] <- varPartStore
  }
}
}


######################################
######### Export figures #############
######################################

lapply(varPartStore.all, 
       function(x) {
         png(paste0(x$SelectedBiome, ".VarPart", ".png"), width = 600, height = 600)
         plot(x$Biome.varpart,
              bg = c("lightgray","lightgray", "lightgray", "gray"), 
              Xnames = c("Current",
                         "Mid-Holocene",
                         "LGM",
                         "Space"),
              main = paste0(x$SelectedBiome),
              id.size = 1.1,
              cex = 1)
         dev.off()
       })

####################################################
##### Export tables and supplementary materials ####
####################################################

#### Export the results of the variation partitioning ####

# X1 - Contemporary Climate
# X2 - Mid-Holocene Climate
# X3 - Last Glaciation Maximum Climate
# X4 - Spatial Structure

# [a], [b], [c], [d], [e], [f], [l], [g]

data.frame.parts <- do.call(rbind,
                            lapply(varPartStore.all, 
                                   function(x) {
                                     rbind(x$Biome.varpart$part$indfract[c(1:15), ])
                                   }
                            )
)

name.vects <- do.call('rbind', 
                      strsplit(as.character(row.names(data.frame.parts)),
                               '.',
                               fixed=TRUE))

varpart.results <- data.frame(Biomes = name.vects[,1],
                              Fraction = name.vects[,2],
                              Explained = round(data.frame.parts$Adj.R.square, 3))
head(varpart.results)

varpart.results <- varpart.results %>%
  mutate(Explained = replace(Explained, 
                             which(Explained < 0), 0)) %>%
  spread(key = "Fraction", 
         value = "Explained") %>%
  mutate(efgl = `[e]` + `[f]` + `[g]` + `[l]`, # 
         jmno = `[j]` + `[m]` + `[n]` + `[o]`) %>%
  gather(key = "Fraction", 
         value = "Explained", -Biomes)

varpart.results$Fraction <- recode(varpart.results$Fraction, 
         "[a] = X1 | X2+X3+X4" = "[Current Climate]",
         "[b] = X2 | X1+X3+X4" = "[Mid-Holocene]", 
         "[c] = X3 | X1+X2+X4" = "[LGM]",
         "[d] = X4 | X1+X2+X3" = "[Space]",
         "efgl" = "[Contemporary:Paleoclimate]",
         "jmno" = "[Contemporary:Paleoclimate:Space]")

varpart.results$Fraction <- as.character(varpart.results$Fraction)


# 
# 
# 
# 
# colnames(varpart.results)[1:5] <- c("Biomes", 
#                                "[Contemporary]", 
#                                "[Mid-Holocene]", 
#                                "[LGM]", 
#                                "[Space]")
# 
# colnames(varpart.results)[17] <- c("[Contemporary and Paleoclimate]")
# colnames(varpart.results)[18] <- c("[Contemporary, Paleoclimate and Space]")

# X1 - Contemporary Climate
# X2 - Mid-Holocene Climate
# X3 - Last Glaciation Maximum Climate

# Common shared between Current and PaleoClimate
# [e] + [g] + [l] + [f]

# Common shared between Current, Paleo, and Space
# [n] + [o] + [m] + [j]

# Export table to the screen
write.table(varpart.results) 

#### Export loadings of the principal component analyses ####

ExplainedVar.PCA.results <- do.call(rbind,
                                    lapply(varPartStore.all, 
                                           function(x) {
                                             cbind.data.frame(x$biome.explained.variation[c(1:2), c(1:3)], 
                                                              c(x$SelectedBiome, x$SelectedBiome),
                                                              row.names(x$biome.explained.variation[c(1:2),]))
                                           }
                                    )
)


write.table(ExplainedVar.PCA.results)

#######################################################################################


ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "fill")

ggplot(data = varpart.results) + 
  geom_bar(mapping = aes(x = Biomes, fill = Explained), position = "fill")


ggplot(data = varpart.results, aes(x = Biomes, y = Explained)) +
  geom_bar(aes(fill = Fraction))

###

varpart.results$Fraction <- revalue(varpart.results$Fraction, c("[p] = Residuals" = "[Residuals]",
                                                                "[a] = X1 | X2+X3+X4" = "[Current Climate]",
                                                                "[b] = X2 | X1+X3+X4" = "[Mid-Holocene]", 
                                                                "[c] = X3 | X1+X2+X4" = "[LGM]",
                                                                "[d] = X4 | X1+X2+X3" = "[Space]"))


varpart.results.b <- varpart.results %>% 
  drop_na() %>%
  filter(Fraction %in% c("[a] = X1 | X2+X3+X4" = "[Current Climate]",
                         "[b] = X2 | X1+X3+X4" = "[Mid-Holocene]", 
                         "[c] = X3 | X1+X2+X4" = "[LGM]",
                         "[d] = X4 | X1+X2+X3" = "[Space]",
                         "efgl" = "[Contemporary:Paleoclimate]",
                         "jmno" = "[Contemporary:Paleoclimate:Space]"))

# varpart.results.b <- varpart.results %>% 
#  dplyr::mutate(Biomes = forcats::fct_reorder(Biomes, as.numeric(Explained), fun = mean))

varpart.results.b$Biomes <- droplevels(varpart.results.b$Biomes)

meanAnnTemp.biomes <- worldClimate %>% 
  select(cur_bio_1) %>%
  bind_cols(regions.LatLong[ , c("ID_Biome_Realm"), drop = FALSE]) %>%
  group_by(ID_Biome_Realm) %>% 
  summarise_all(funs(mean)) %>%
  arrange(desc(cur_bio_1)) %>%
  as.data.frame()

meanAnnTemp.biomes$ID_Biome_Realm <- as.factor(meanAnnTemp.biomes$ID_Biome_Realm)

varpart.results.b <- varpart.results.b %>%
  left_join(meanAnnTemp.biomes, c("Biomes" = "ID_Biome_Realm"))

varpart.results.b$Biomes <- as.factor(varpart.results.b$Biomes)


varpart.results.b$Biomes <- reorder(varpart.results.b$Biomes, varpart.results.b$cur_bio_1)

# nrow(meanAnnTemp.biomes[!grepl("NA__", meanAnnTemp.biomes$ID_Biome_Realm),])

varpart.results.b$Biomes <- factor(varpart.results.b$Biomes, varpart.results.b$ID_Biome_Realm)

varpart.results.b$ID_Biome_Realm <- varpart.results.b$Biomes

varpart.results.b$Biome_Acronym = factor(varpart.results.b$Biomes, 
                                         labels = abbreviate(gsub("_", 
                                                                  " ",
                                                                  levels(varpart.results.b$Biomes))))
ggplot(data = varpart.results.b, 
       aes(x = Biome_Acronym, 
           y = Explained, 
           fill = Fraction)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "") +
  scale_fill_viridis_d(option = "E", 
                       direction = 1) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(face = "bold", 
                                   size = 14),
        axis.title = element_text(size = 16, 
                                  face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
  )

# Bio_1=Annual Mean Temperature [°C*10]
# Bio_2=Mean Diurnal Range [°C]
# Bio_3=Isothermality [Bio_2/Bio_7]
# Bio_4=Temperature Seasonality [standard deviation*100]
# Bio_5=Max Temperature of Warmest Month [°C*10]
# Bio_6=Min Temperature of Coldest Month [°C*10]
# Bio_7=Temperature Annual Range [°C*10]
# Bio_8=Mean Temperature of Wettest Quarter [°C*10]
# Bio_9=Mean Temperature of Driest Quarter [°C*10]
# Bio_10=Mean Temperature of Warmest Quarter [°C*10]
# Bio_11=Mean Temperature of Coldest Quarter [°C*10]
# Bio_12=Annual Precipitation [mm/year]
# Bio_13=Precipitation of Wettest Month [mm/month]
# Bio_14=Precipitation of Driest Month [mm/month]
# Bio_15=Precipitation Seasonality [coefficient of variation]
# Bio_16=Precipitation of Wettest Quarter [mm/quarter]
# Bio_17=Precipitation of Driest Quarter [mm/quarter]
# Bio_18=Precipitation of Warmest Quarter [mm/quarter]
# Bio_19=Precipitation of Coldest Quarter [mm/quarter]


d.frac.AnnTemp <- varpart.results.b %>% 
  left_join(meanAnnTemp.biomes) %>%
  filter(Fraction == "[Space]")

c.frac.AnnTemp <- varpart.results.b %>% 
  left_join(meanAnnTemp.biomes) %>%
  filter(Fraction == "[LGM]")

a.frac.AnnTemp <- varpart.results.b %>% 
  left_join(meanAnnTemp.biomes) %>%
  filter(Fraction == "[Current Climate]")

b.frac.AnnTemp <- varpart.results.b %>% 
  left_join(meanAnnTemp.biomes) %>%
  filter(Fraction == "[Mid-Holocene]")

efgl.frac.AnnTemp <- varpart.results.b %>% 
  left_join(meanAnnTemp.biomes) %>%
  filter(Fraction == "[Contemporary:Paleoclimate]")

efgl.frac.AnnTemp

ggplot(efgl.frac.AnnTemp, aes(x=cur_bio_1, y=Explained)) + 
  geom_point()+
  geom_smooth(method = lm)

plot(d.frac.AnnTemp$cur_bio_1, d.frac.AnnTemp$Explained)
plot(c.frac.AnnTemp$cur_bio_1, c.frac.AnnTemp$Explained)

plot(d.frac.AnnTemp$cur_bio_1, d.frac.AnnTemp$Explained)
plot(c.frac.AnnTemp$cur_bio_1, c.frac.AnnTemp$Explained)



