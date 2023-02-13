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
library(dplyr)

TableS1 <- MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
  select(ID, SamplingPool, nri, nti, ID_Realm) %>%
  drop_na() %>%
  group_by(ID_Realm, SamplingPool)
  
summary.NRI.NTI <- bind_rows(TableS1 %>%
            filter(SamplingPool == "Global sampling") %>%
            ungroup() %>%
            summarise(mean_NTI = mean(nti, na.rm = TRUE), 
                      sd_NTI = sd(nti, na.rm = TRUE),
                      mean_NRI = mean(nri, na.rm = TRUE), 
                      sd_NRI = sd(nri, na.rm = TRUE)) %>%
            add_column(ID_Realm = "Global", 
                       SamplingPool = "Global", 
                       .before = TRUE),
          TableS1 %>%
            summarise(mean_NTI = mean(nti, na.rm = TRUE), 
                      sd_NTI = sd(nti, na.rm = TRUE),
                      mean_NRI = mean(nri, na.rm = TRUE), 
                      sd_NRI = sd(nri, na.rm = TRUE))) %>%
  as.data.frame() %>%
  mutate(mean_NTI.sd_NTI = paste(round(mean_NTI, 2), " (±", round(sd_NTI, 2),")", sep = ""),
         mean_NRI.sd_NRI = paste(round(mean_NRI, 2), " (±", round(sd_NRI, 2),")", sep = "")) %>%
  select(-mean_NTI, -sd_NTI,
         -mean_NRI, -sd_NRI)

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


