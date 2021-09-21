library(quantreg)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggarrange)


# Obtain predicted values for the same x

pred.rq.poly.4.NRI.netDiv_CWM_std_tw <- as.data.frame(predict(rq.poly.4.NRI.netDiv_CWM_std_tw, 
           newdata = prebinning_test %>% select(netDiv_CWM_std_tw)))

pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw <- prebinning_test$netDiv_CWM_std_tw

# Reshape predictions to long format
pred.rq.poly.4.NRI.netDiv_CWM_std_tw<- melt(pred.rq.poly.4.NRI.netDiv_CWM_std_tw, 
                                            id.var = "netDiv_CWM_std_tw",
                                            variable.name = "tau",
                                            value.name = "NRI")

pred.rq.poly.4.NRI.netDiv_CWM_std_tw$tau <- as.numeric(gsub("tau= (.*)", "\\1", 
                                                            pred.rq.poly.4.NRI.netDiv_CWM_std_tw$tau))
pred.rq.poly.4.NRI.netDiv_CWM_std_tw$label <- NA
pred.rq.poly.4.NRI.netDiv_CWM_std_tw$label[which(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw == max(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw))] <- pred.rq.poly.4.NRI.netDiv_CWM_std_tw$tau[which(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw == max(pred.rq.poly.4.NRI.netDiv_CWM_std_tw$netDiv_CWM_std_tw))]

(NRI.poly.4.netDiv_CWM_std_tw.rq.plot <- ggplot(data = pred.rq.poly.4.NRI.netDiv_CWM_std_tw, 
       aes(x = netDiv_CWM_std_tw, 
           y = NRI,
           colour = tau, group = factor(tau)))  +   
  geom_line(size = 1) + 
  geom_text_repel(
    data = pred.rq.poly.4.NRI.netDiv_CWM_std_tw %>%
      drop_na() %>%
      group_by(tau) %>%
      summarise(tau = mean(tau),
                NRI = mean(NRI),
                netDiv_CWM_std_tw = mean(netDiv_CWM_std_tw),
                label = mean(label)
      ),
    aes(label = paste(tau*100,"%")),
    size = 4.5,
    nudge_x = 0.05,
    segment.color = "grey",
    min.segment.length = 0.1,
    colour = "black"
  )  +
  geom_hline(yintercept = 0, alpha = 0.25) +  
  scale_x_continuous( limits = c(c(-0.025, 0.05) + range(prebinning_test$netDiv_CWM_std_tw))) +
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

ggsave(filename = "figures/Fig.X.rq.poly.4.NRI.netDiv_CWM_std_tw_clean.png", 
       dpi = 300, 
       width = 8.5 , height = 8.5, 
       units = "in")

# Obtain predicted values for the same x

pred.rq.poly.4.NTI.netDiv_CWM_std_tw <- as.data.frame(predict(rq.poly.4.NTI.netDiv_CWM_std_tw, 
                                                              newdata = prebinning_test %>% select(netDiv_CWM_std_tw)))

pred.rq.poly.4.NTI.netDiv_CWM_std_tw$netDiv_CWM_std_tw <- prebinning_test$netDiv_CWM_std_tw

# Reshape predictions to long format
pred.rq.poly.4.NTI.netDiv_CWM_std_tw<- melt(pred.rq.poly.4.NTI.netDiv_CWM_std_tw, 
                                            id.var = "netDiv_CWM_std_tw",
                                            variable.name = "tau",
                                            value.name = "NTI")

pred.rq.poly.4.NTI.netDiv_CWM_std_tw$tau <- as.numeric(gsub("tau= (.*)", "\\1", 
                                                            pred.rq.poly.4.NTI.netDiv_CWM_std_tw$tau))
pred.rq.poly.4.NTI.netDiv_CWM_std_tw$label <- NA
pred.rq.poly.4.NTI.netDiv_CWM_std_tw$label[which(pred.rq.poly.4.NTI.netDiv_CWM_std_tw$netDiv_CWM_std_tw == max(pred.rq.poly.4.NTI.netDiv_CWM_std_tw$netDiv_CWM_std_tw))] <- pred.rq.poly.4.NTI.netDiv_CWM_std_tw$tau[which(pred.rq.poly.4.NTI.netDiv_CWM_std_tw$netDiv_CWM_std_tw == max(pred.rq.poly.4.NTI.netDiv_CWM_std_tw$netDiv_CWM_std_tw))]

(NTI.poly.4.netDiv_CWM_std_tw.rq.plot <- ggplot(data = pred.rq.poly.4.NTI.netDiv_CWM_std_tw, 
       aes(x = netDiv_CWM_std_tw, 
           y = NTI,
           colour = tau, group = factor(tau)))  +   
  geom_line(size = 1) + 
  geom_text_repel(
    data = pred.rq.poly.4.NTI.netDiv_CWM_std_tw %>%
      drop_na() %>%
      group_by(tau) %>%
      summarise(tau = mean(tau),
                NTI = mean(NTI),
                netDiv_CWM_std_tw = mean(netDiv_CWM_std_tw),
                label = mean(label)
      ),
    aes(label = paste(tau*100,"%")),
    size = 4.5,
    nudge_x = 3,
    segment.color = "grey",
    min.segment.length = 0.1,
    colour = "black"
  )  +
    geom_hline(yintercept = 0, alpha = 0.25) +  
  scale_x_continuous( limits = c(c(-0.025, 0.05) + range(prebinning_test$netDiv_CWM_std_tw))) +
  scico::scale_colour_scico(palette = "acton", direction = -1) + 
  # facet_grid(. ~ ID_Realm,
  #  scales = "free") +
  labs(y = c(expression("NTI"["Global"])),
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

ggsave(filename = "figures/Fig.X.rq.poly.4.NTI.netDiv_CWM_std_tw_clean.png", 
      dpi = 300, 
      width = 8.5 , height = 8.5, 
      units = "in")


ggarrange(NRI.poly.4.netDiv_CWM_std_tw.rq.plot + theme(axis.title.x = element_blank() #,
                                                       # axis.text = element_text(size = 12),
                                                       # axis.title.y = element_text(size = 12)
                                                       ),
          NTI.poly.4.netDiv_CWM_std_tw.rq.plot + theme(# axis.text = element_text(size = 12),
                                                       #axis.title.y = element_text(size = 12)
            ),
          labels = c("A", "B"),
          ncol = 1, nrow = 2, # hjust = 0.25,
          align = "v",
          heights = c(0.9, 1),
          common.legend = FALSE)

ggsave(filename = "NRI.NTI.poly.4.netDiv_CWM_std_tw.rq.plot.png", 
       dpi = 200,
       width = 8*1,
       height = 8*2, 
       units = "in")

