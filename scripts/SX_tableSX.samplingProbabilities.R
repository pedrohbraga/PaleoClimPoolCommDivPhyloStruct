#############################################################################################
#### Code to prepare reproducible HTML tables with the results obtained in this project  ####
#                                                                                           #
# Author: Pedro Henrique Pereira Braga                                                      #
# Last Update: "2020-03-20 22:00:30 EST"
#
#############################################################################################

ASM_speciesRichnessGenus_tib <- ASM_BatsData %>% 
    group_by(Genus) %>% 
    summarize(count_ASM=n()) 

FaurSven_speciesRichnessGenus_tib <- sppTable %>%
  group_by(V2) %>% 
  summarize(count_FaurSven=n()) 

samplingProbsBats_tib <- FaurSven_speciesRichnessGenus_tib %>%
  left_join(ASM_speciesRichnessGenus_tib,
            by = c("V2" = "Genus")) %>%
  mutate(SamplingFraction = count_FaurSven/count_ASM) %>%
  mutate(SamplingFraction = ifelse(SamplingFraction > 1, 1, SamplingFraction)) %>%
  as.data.frame()

## Generating a table 

kable(samplingProbsBats_tib,
      row.names = FALSE,
      col.names = c("Clade (Genus)",
                    "Species kept from Faurby & Svenning (2015)",
                    "Species within ASM",
                    "Sampling fraction"),
      align=c("l", 
              "c",
              "c",
              "l")) %>%
  kable_styling(full_width = F)  %>%
#  collapse_rows(row_group_label_position = 'stack', 
#                valign = "top") %>%
  save_kable(file = "tableS1.rich.samplingProbsBats.html", 
             self_contained = T)


# Create a table using knitr::kable

kable(samplingProbsBats[-1, -1] %>%
        mutate(prob = round(as.numeric(prob), 4)*100) %>%
        mutate(prob = paste(prob, "%", sep = " ")), 
      row.names = FALSE,
      col.names = c("Clade (Genus)",
                    "Sampling Fraction"),
      align=c("l", "c")) %>%
  kable_styling(full_width = F)  %>%
  collapse_rows(row_group_label_position = 'stack', 
                valign = "top") %>%
  save_kable(file = "tableS1.samplingProbsBats.html", 
             self_contained = T)

