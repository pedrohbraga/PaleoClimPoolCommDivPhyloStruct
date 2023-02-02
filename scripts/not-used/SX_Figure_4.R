(NRI.poly.4.netDiv_CWM_std_tw.rq.plot <- ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
         filter(SamplingPool == "Global sampling")%>%
         drop_na(ID_Realm), 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI)) +
  stat_quantile(quantiles = q10, 
                formula = y ~ poly(x, 4), 
                colour = "black"
  ) +
  xlim(min(prebinning_test$netDiv_CWM_std_tw),
       0.35) +
  # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic(base_size = 17) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1))
  )


ggsave(filename = "figures/Fig.X.rq.poly.4.NRI.netDiv_CWM_std_tw_clean.png", 
       dpi = 300, 
       width = 8.5 , height = 8.5, 
       units = "in")


(NTI.poly.4.netDiv_CWM_std_tw.rq.plot <- ggplot(data = CWM.Div.MPD.MNTD.Chiroptera.Comm %>%
                                                 filter(SamplingPool == "Global sampling")%>%
                                                 drop_na(ID_Realm), 
                                               aes(x = netDiv_CWM_std_tw, 
                                                   y = NTI)) +
  stat_quantile(quantiles = q10, 
                formula = y ~ poly(x, 4), 
                colour = "black"
  ) +
  xlim(min(prebinning_test$netDiv_CWM_std_tw),
       0.35) +
  # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(x = c(expression("Net Diversification Rate"[CWM[STD[tw]]]))) +
  theme_classic(base_size = 17) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.box = "horizontal") + 
  guides(colour = guide_legend(nrow = 1)
         )
  )


ggarrange(NRI.poly.4.netDiv_CWM_std_tw.rq.plot + theme(axis.title.x = element_blank() #,
                                                       # axis.text = element_text(size = 12),
                                                       # axis.title.y = element_text(size = 12)
                                                       ),
          NTI.poly.4.netDiv_CWM_std_tw.rq.plot + theme(# axis.text = element_text(size = 12),
                                                       #axis.title.y = element_text(size = 12)
            ),
          ncol = 1, nrow = 2, hjust = 0.25,
          align = "hv",
          common.legend = TRUE)


ggsave(filename = "NRI.NTI.poly.4.netDiv_CWM_std_tw.rq.plot.png", 
       dpi = 300,
       width = 6*1,
       height = 6*2, 
       units = "in")

