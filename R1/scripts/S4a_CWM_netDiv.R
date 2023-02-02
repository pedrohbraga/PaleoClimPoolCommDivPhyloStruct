#################################################################################
### Code to obtain local diversification rates for each community based on    ### 
### rates estimated for each tip in the phylogeny                             ###
#                                                                               #
#                                                                               #
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2021-05-26"                                                     #
#                                                                               # 
#################################################################################

# Check lengths
length(tip.rates$lambda.avg); length(tip.rates$mu.avg)

#### Calculating CWM for Net Diversification Rates ####

### CWM for Speciation Rates

lambda.CWM.Chiroptera.Comm <- CWM_Std_TW(Trait = tip.rates$lambda.avg, 
                                         Distrib = Chiroptera.FaurSven.comm)

colnames(lambda.CWM.Chiroptera.Comm) <- paste("lambda", 
                                              colnames(lambda.CWM.Chiroptera.Comm),
                                              "rs_1", 
                                              sep = "_")

mu.CWM.Chiroptera.Comm <- CWM_Std_TW(Trait = tip.rates$mu.avg, 
                                     Distrib = Chiroptera.FaurSven.comm)

colnames(mu.CWM.Chiroptera.Comm) <- paste("mu", 
                                          colnames(mu.CWM.Chiroptera.Comm), 
                                          "rs_1",
                                          sep = "_")

## Subtract mu from lambda to obtain per community net diversification rates
# (note that this is different then subtracting mu from lambda directly at the
# tip-rates)

netDiv.CWM.Chiroptera.Comm <- lambda.CWM.Chiroptera.Comm - mu.CWM.Chiroptera.Comm

colnames(netDiv.CWM.Chiroptera.Comm) <- c("netDiv_CWM_rs_1", 
                                          "netDiv_CWM_std_tw_rs_1", 
                                          "netDiv_CWM_std_w_rs_1" )

netDiv.CWM.Chiroptera <- data.frame(lambda.CWM.Chiroptera.Comm,
                                    mu.CWM.Chiroptera.Comm, 
                                    netDiv.CWM.Chiroptera.Comm) %>% 
  rownames_to_column("ID")

netDiv.CWM.Chiroptera$ID <- as.numeric(netDiv.CWM.Chiroptera$ID)

# Joining netDiv.CWM.Chiroptera to the MPD.MNTD data

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff %>%
  left_join(netDiv.CWM.Chiroptera, 
            by = c("ID" = "ID")
  ) # %>%
# left_join(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.rarefaction.relative %>%
#             dplyr::select(ID_SamplingPool, 
#                           nri.rarefac.mean, 
#                           nti.rarefac.mean, 
#                           nri.sec, 
#                           nti.sec), 
#           by = c("ID_SamplingPool" = "ID_SamplingPool")
# ) 

# Not used ####

# ### ClaDS ###
# 
# #### Calculating CWM for Net Diversification Rates ####
# 
# ### CWM for Speciation Rates
# 
# lambda.CWM.Chiroptera.Comm.ClaDS <- CWM_Std_TW(Trait = tip.rates.lambda.ClaDS, 
#                                                Distrib = Chiroptera.FaurSven.comm)
# 
# colnames(lambda.CWM.Chiroptera.Comm.ClaDS) <- paste("lambda", 
#                                                     colnames(lambda.CWM.Chiroptera.Comm.ClaDS), 
#                                                     "ClaDS",  
#                                                     sep = "_")
# 
# lambda.CWM.Chiroptera.Comm.ClaDS <- lambda.CWM.Chiroptera.Comm.ClaDS  %>% 
#   rownames_to_column("ID")
# 
# lambda.CWM.Chiroptera.Comm.ClaDS$ID <- as.numeric(lambda.CWM.Chiroptera.Comm.ClaDS$ID)
# 
# MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div <- MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
#   left_join(lambda.CWM.Chiroptera.Comm.ClaDS, 
#             by = c("ID" = "ID"))


