#############################################################################################
#### Code to prepare reproducible HTML tables with the results obtained in this project  ####
#                                                                                           #
# Author: Pedro Henrique Pereira Braga                                                      #
# Last Update: "2019-06-20 17:00:30 EST"
#
#                                                                                           #
#############################################################################################

library(knitr)
library(kableExtra)

# Merge both NRI and NTI data

NRI.NTI.Realm <- filter(MPD.LatLong.Env.AllScales,
                        is.na(ID_Realm) == FALSE) %>%
  right_join(filter(MNTD.LatLong.Env.AllScales,
                    is.na(ID_Realm) == FALSE), 
             by = c("ID", "SamplingPool")) %>%
  group_by(ID_Realm.x, SamplingPool) %>%
  select(ID_Realm.x, NTI, NRI)


# Summarise the NRI and NTI data for the globe

summary.NRI.NTI.Global <- NRI.NTI.Realm %>%
  filter(SamplingPool == "Global sampling") %>%
  ungroup() %>%
  summarise(mean_NTI = mean(NTI, na.rm = TRUE), 
            sd_NTI = sd(NTI, na.rm = TRUE),
            mean_NRI = mean(NRI, na.rm = TRUE), 
            sd_NRI = sd(NRI, na.rm = TRUE)) %>%
  add_column(ID_Realm.x = "Global", 
             SamplingPool = "Global", 
             .before = TRUE)

# Summarise NRI and NTI data per Realm

summary.NRI.NTI.Realm <- NRI.NTI.Realm %>%
  summarise(mean_NTI = mean(NTI, na.rm = TRUE), 
            sd_NTI = sd(NTI, na.rm = TRUE),
            mean_NRI = mean(NRI, na.rm = TRUE), 
            sd_NRI = sd(NRI, na.rm = TRUE))

# Row-bind both global and per-realm summaries

summary.NRI.NTI <- summary.NRI.NTI.Global %>%
  bind_rows(summary.NRI.NTI.Realm) %>%
  as.data.frame() %>%
  mutate(mean_NTI.sd_NTI = paste(round(mean_NTI, 2), " (±", round(sd_NTI, 2),")", sep = ""),
         mean_NRI.sd_NRI = paste(round(mean_NRI, 2), " (±", round(sd_NRI, 2),")", sep = "")) %>%
  select(-mean_NTI, -sd_NTI,
         -mean_NRI, -sd_NRI)
      
# Create a table using knitr::kable

kable(summary.NRI.NTI, 
      row.names = FALSE,
      col.names = c("Spatial extent",
                    "Sampling Pool",
                    "NTI_mean (±NTI_sd)",
                    "NRI_mean (±NRI_sd)"),
      align=c("l", "l", "c", "c")) %>%
  kable_styling(full_width = F)  %>%
  collapse_rows(row_group_label_position = 'stack', 
                valign = "top") %>%
  save_kable(file = "tableS1.summary.NRI.NTI.html", 
             self_contained = T)

# write.table(summary.NRI.NTI)


