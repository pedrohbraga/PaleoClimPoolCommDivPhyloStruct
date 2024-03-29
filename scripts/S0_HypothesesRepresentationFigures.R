######################################################################################
### Code to create figures to represent the hypotheses and predictions for this    ###
##  study                                                                           ##
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2019-01-28"                                                          #
#                                                                                    # 
######################################################################################

# Alternative H1 - Strong dispersal barriers between regions lead to within-realm diversification --------

Hypothesis.1.alt <- data.frame(SamplingPool = factor(rep(c("Broad scale",
                                                           "Regional scale",
                                                           "Biome scale",
                                                           "Ecoregion scale"), 
                                                         each = 2000), 
                                                     levels = c("Broad scale",
                                                                "Regional scale",
                                                                "Biome scale",
                                                                "Ecoregion scale")
                                                     ), 
                               PhyloStructure = c(rnorm(2000, 
                                                        sd = 5, 
                                                        mean = 18),
                                                  rnorm(2000, 
                                                        sd = 1.8, 
                                                        mean = 9),
                                                  rnorm(2000, 
                                                        sd = 1.5, 
                                                        mean = 5),
                                                  rnorm(2000, 
                                                        sd = 1, 
                                                        mean = 2)))

(Hypothesis.1.alt.plot <- ggplot(Hypothesis.1.alt, 
                                 aes(x=SamplingPool, y=PhyloStructure)) + 
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    ylim(1, 20) +
    #    geom_hline(yintercept = 0, alpha = 0.4) +
    labs(x = "⬅ Sampling Pool Restriction ➡",
         y = "Phylogenetic Relatedness") +
    theme(strip.background=element_rect(fill = "white",
                                        linetype = NULL,
                                        color = "white"),
          strip.text=element_text(color = "black",
                                  face = "bold",
                                  size = 15),
          #axis.text.x = element_text(size = 15, 
          #                           face="bold"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 16, 
                                    face="bold"),
          axis.title.x = element_text(size = 16, 
                                      face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
    )
)


ggsave(filename = "figures/Hypothesis.1.alt.png", 
       dpi = 300, 
       width = 4.5, height = 5, 
       units = "in")

# Null H1 - Strong dispersal barriers between regions lead to within-realm diversification --------

Hypothesis.1.null <- data.frame(SamplingPool = factor(rep(c("Broad scale",
                                                            "Regional scale",
                                                            "Biome scale",
                                                            "Ecoregion scale"), 
                                                          each = 2000), levels=c("Broad scale",
                                                                                 "Regional scale",
                                                                                 "Biome scale",
                                                                                 "Ecoregion scale")), 
                                PhyloStructure = c(rnorm(2000, 
                                                         sd = 3, 
                                                         mean = 2.5),
                                                   rnorm(2000, 
                                                         sd = 1.8, 
                                                         mean = 2.4),
                                                   rnorm(2000, 
                                                         sd = 1.5, 
                                                         mean = 2.3),
                                                   rnorm(2000, 
                                                         sd = 1, 
                                                         mean = 2)))

(Hypothesis.1.null.plot <- ggplot(Hypothesis.1.null, 
                                  aes(x=SamplingPool, y=PhyloStructure)) + 
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    ylim(0, 20) +
    #    geom_hline(yintercept = 0, alpha = 0.4) +
    labs(x = "⬅ Sampling Pool Restriction ➡",
         y = "Phylogenetic Relatedness") +
    theme(strip.background=element_rect(fill = "white",
                                        linetype = NULL,
                                        color = "white"),
          strip.text=element_text(color = "black",
                                  face = "bold",
                                  size = 15),
          #axis.text.x = element_text(size = 15, 
          #                           face="bold"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 16, 
                                    face="bold"),
          axis.title.x = element_text(size = 16, 
                                      face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
    )
)


ggsave(filename = "figures/Hypothesis.1.null.plot.png", 
       dpi = 300, 
       width = 4.5, height = 5, 
       units = "in")


# Alternative H2 - Lineages inhabiting more stable climate over extended geological time should show higher net diversification rates through time --------

Hypothesis.2.alt <- data.frame(PhyloStructure = 1:2000)

Hypothesis.2.alt$ClimaticStability.alt = 2 + 5*Hypothesis.2.alt$PhyloStructure + rnorm(2000, sd = 20)

Hypothesis.2.alt$ClimaticStability.exp = Hypothesis.2.alt$PhyloStructure^2

ggplot(data = Hypothesis.2.alt, aes(x = ClimaticStability.exp, y = PhyloStructure)) +
  stat_smooth(method = "glm", 
              se = FALSE,
              colour = "black",
              method.args = list(family = gaussian(link = "log"))) +
  # theme_minimal() +
  theme_classic() +
  #    geom_hline(yintercept = 0, alpha = 0.4) +
  labs(x = "Climatic Stability",
       y = "Phylogenetic Relatedness") +
  theme(strip.background = element_rect(fill = "white",
                                        linetype = NULL,
                                        color = "white"),
        strip.text = element_text(color = "black",
                                  face = "bold",
                                  size = 15),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16, 
                                  face="bold"),
        axis.title.x = element_text(size = 16, 
                                    face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
  )

ggsave(filename = "figures/Hypothesis.2.alt.png", 
       dpi = 300, 
       width = 5, height = 5, 
       units = "in")

# Alternative H3 - Increased in situ diversification rates creates clusters of closely-related species, leading to high levels of phylogenetic clustering --------

Hypothesis.3.alt <- data.frame(PhyloStructure = 1:2000)

Hypothesis.3.alt$ClimaticStability.alt = 2 + 5*Hypothesis.3.alt$PhyloStructure + rnorm(2000, sd = 20)

Hypothesis.3.alt$ClimaticStability.exp = Hypothesis.3.alt$PhyloStructure^2

ggplot(data = Hypothesis.3.alt, aes(x = ClimaticStability.exp, y = PhyloStructure)) +
  stat_smooth(method = "glm", se = FALSE,
              colour = "black",
              method.args = list(family = gaussian(link = "log")))+
  # theme_minimal() +
  theme_classic() +
  #    geom_hline(yintercept = 0, alpha = 0.4) +
  labs(x = "Local Diversification Rates",
       y = "Phylogenetic Relatedness") +
  theme(strip.background = element_rect(fill = "white",
                                        linetype = NULL,
                                        color = "white"),
        strip.text = element_text(color = "black",
                                  face = "bold",
                                  size = 15),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16, 
                                  face="bold"),
        axis.title.x = element_text(size = 16, 
                                    face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
  )

ggsave(filename = "figures/Hypothesis.3.alt.png", 
       dpi = 300, 
       width = 5, height = 5, 
       units = "in")
