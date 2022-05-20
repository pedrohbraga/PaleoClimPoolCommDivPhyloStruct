# Unconditional quantile summary statistics

quantile.XY <- function(dataset = dataset,
                        X, Y, 
                        n.classes = 100, 
                        lab_x = "X", 
                        lab_y = "Y") {
  

  X = dataset %>% as.data.frame() %>% pull(X)
  Y = dataset %>% as.data.frame() %>% pull(Y)
  
  # divide X into n.classes quantiles and then calculate the variance of Y in
  # each interval.
  
  quantile.values <- cumsum(rep(1 / (n.classes + 1), n.classes))
  X.quant <- quantile(X, quantile.values)
  matrix.summaries <- matrix(0, n.classes, 13)
  colnames(matrix.summaries) <- c("median", "mean", "variance", 
                                  "sd", "se", "IQR", 
                                  "min", "max",
                                  "min.X", "max.X", "mean.X", "sd.X", "se.X")
  for (i in 1:n.classes) {
    if (i == 1) {
      matrix.summaries[i, "median"] <- median(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "mean"] <- mean(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "variance"] <- var(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "sd"] <- sd(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "se"] <- plotrix::std.error(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "IQR"] <- IQR(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "min"] <- min(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "max"] <- max(Y[which(X <= X.quant[i])])
      matrix.summaries[i, "min.X"] <- min(X[which((X <= X.quant[i]))])
      matrix.summaries[i, "max.X"] <- max(X[which((X <= X.quant[i]))])
      matrix.summaries[i, "mean.X"] <- mean(X[which((X <= X.quant[i]))])
      matrix.summaries[i, "sd.X"] <- sd(X[which((X <= X.quant[i]))])
      matrix.summaries[i, "se.X"] <- plotrix::std.error(X[which((X <= X.quant[i]))])
    }
    if (i == n.classes) {
      matrix.summaries[i, "median"] <- median(Y[which(X > X.quant[i])])
      matrix.summaries[i, "mean"] <- mean(Y[which(X > X.quant[i])])
      matrix.summaries[i, "variance"] <- var(Y[which(X > X.quant[i])])
      matrix.summaries[i, "sd"] <- sd(Y[which(X > X.quant[i])])
      matrix.summaries[i, "se"] <- plotrix::std.error(Y[which(X > X.quant[i])])
      matrix.summaries[i, "IQR"] <- IQR(Y[which(X > X.quant[i])])
      matrix.summaries[i, "min"] <- min(Y[which(X > X.quant[i])])
      matrix.summaries[i, "max"] <- max(Y[which(X > X.quant[i])])
      matrix.summaries[i, "min.X"] <- min(X[which((X > X.quant[i]))])
      matrix.summaries[i, "max.X"] <- max(X[which((X > X.quant[i]))])
      matrix.summaries[i, "mean.X"] <- mean(X[which((X > X.quant[i]))])
      matrix.summaries[i, "sd.X"] <- sd(X[which((X > X.quant[i]))])
      matrix.summaries[i, "se.X"] <- plotrix::std.error(X[which((X > X.quant[i]))])
    }
    if ((i > 1) & (i < n.classes)) {
      matrix.summaries[i, "median"] <- median(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "mean"] <- mean(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "variance"] <- var(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "sd"] <- sd(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "se"] <- plotrix::std.error(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "IQR"] <- IQR(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))], na.rm = TRUE)
      matrix.summaries[i, "min"] <- min(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))], na.rm = TRUE)
      matrix.summaries[i, "max"] <- max(Y[which((X > X.quant[i]) & (X <= X.quant[i + 1]))], na.rm = TRUE)
      matrix.summaries[i, "min.X"] <- min(X[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "max.X"] <- max(X[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "mean.X"] <- mean(X[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "sd.X"] <- sd(X[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
      matrix.summaries[i, "se.X"] <- plotrix::std.error(X[which((X > X.quant[i]) & (X <= X.quant[i + 1]))])
    }
  }
  
  # data.df <- data.frame(X.quantiles = X.quant,mean.Y=matrix.summaries[,"mean"],sd.Y=matrix.summaries[,"sd"])
  
  # plot.list <- list()
  # plot.list[[1]] <- ggplot(data.df, aes(x=X.quantiles, y=mean.Y)) +
  #   geom_point() + theme_classic() +
  #   theme(axis.text.y = element_text(size = 14)) +
  #   ylab("NRI (mean per quantile of X)") +
  #   theme(axis.title.y.left = element_text(size = 14)) +
  #   theme(axis.title.x = element_blank())
  
  # plot.list[[2]] <- ggplot(data.df, aes(x=X.quantiles, y=sd.Y)) +
  #   geom_point() + theme_classic() +
  #   theme(axis.text.y = element_text(size = 14)) +
  #   ylab("NRI (std. deviations per quantile of X)") +
  #   theme(axis.title.y.left = element_text(size = 14)) +
  #   theme(axis.title.x = element_blank())
  
  # panel.biotic <- ggarrange(plotlist=plot.list,ncol=2,nrow=1)
  # panel.biotic <- annotate_figure(panel.biotic,bottom=text_grob("Historical Change in Temperature", face = "bold", color = "black", size=14,vjust=0.5,hjust=0.45))+
  #   theme(plot.margin=grid::unit(c(1,1,1,1), "cm"))
  # return(panel.biotic)
  
  data.df <- data.frame(X.quantiles = X.quant, 
                        mean.Y = matrix.summaries[, "mean"],
                        sd.Y = matrix.summaries[, "sd"],
                        se.Y = matrix.summaries[, "se"],
                        IQR.Y = matrix.summaries[, "IQR"],
                        min.Y = matrix.summaries[, "min"],
                        max.Y = matrix.summaries[, "max"],
                        min.X = matrix.summaries[, "min.X"],
                        max.X = matrix.summaries[, "max.X"],
                        mean.X = matrix.summaries[, "mean.X"],
                        sd.X = matrix.summaries[, "sd.X"],
                        se.X = matrix.summaries[, "se.X"])
  
  # Represent with ggplot() and geom_point()
  
  (fig.subset.mean.var_y.var_x.quantile <- ggplot(
    data.df,
    aes(
      x = mean.X,
      y = mean.Y
    )
  ) +
      # geom_errorbar(aes(ymin = mean.Y - 1.96*se.Y,
      #                   ymax = mean.Y + 1.96*se.Y),
      #               cex = 0.1) +
      # geom_errorbarh(aes(xmin = min.X,
      #                    xmax = max.X),
      #                cex = 1) +
      # scico::scale_fill_scico_d(
      #   palette = "acton",
      #   direction = -1,
      #   drop = FALSE
      # ) +
    geom_point(cex = 2.1) +
      labs(
        x = lab_x,
        y = lab_y
      ) +
      geom_hline(yintercept = 0, alpha = 0.25) +
      geom_vline(xintercept = 0, alpha = 0.25) +
      scale_y_continuous(
        breaks = pretty(c(data.df$mean.Y, 
                          0), n = 6),
        limits = c(-0.05, 0.05) + range(pretty(c(data.df$mean.Y, 0), n = 7)),
        #   expand = expansion(mult = c(0, 0)),
        position = "left"
      ) +
      scale_x_continuous(
        breaks = pretty(c(X, 0), n = 5),
        limits = c(-0.05, 0.05) + range(pretty(c(data.df$mean.X, 0), 
                                               n = 6)),
        #        expand = expansion(mult = c(0, 0)),
      ) +
      theme_classic(base_size = 17 * 1.6) +
      theme(
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.box = "horizontal",
        # axis.title.x = element_blank(),
        axis.text.x = element_text(size = 17 * 1.6),
        axis.text.y = element_text(size = 17 * 1.6),
        axis.title.y = element_text(
          size = 17 * 2,
          face = "bold"
        ),
        axis.title.x = element_text(
          size = 17 * 1.4,
          face = "bold"
        ),
        plot.margin = unit(c(0.5, 0.6, 0, 0.6), "cm")
      ) +
      guides(colour = guide_legend(nrow = 1))
  )
  return(fig.subset.mean.var_y.var_x.quantile)
}

# p-value computation for bootstrap estimations using confidence interval inversion

boot.pvalue.CI_inversion <- function(boot.vector, theta_null = 0) {
  
  # as described in Hall (1992)
  n.boot <- length(boot.vector)
  pval_precision <- 1 / n.boot
  alpha_seq <- seq(1e-16, 1 - 1e-16, 
                   pval_precision)
  CIs <- matrix(0, length(alpha_seq), 5)
  conf <- 1 - alpha_seq
  
  # this need some optimization to make it faster when n.boot is large
  for (i in 1:length(conf)) {
    CIs[i, 1] <- conf[i]
    CIs[i, 2:3] <- c(alpha_seq[i] / 2, 1 - alpha_seq[i] / 2)
    CIs[i, 4:5] <- quantile(boot.vector, 
                            c(alpha_seq[i] / 2, 
                              1 - alpha_seq[i] / 2)
    ) # CI[1,5] should be NA but it works fine
  }
  
  bounds <- CIs[, 4:5] # bounds of CI
  p.value <- alpha_seq[which.min(theta_null >= bounds[, 1] & theta_null <= bounds[, 2])]
  
  # perhaps consider bias-corrected and accelerated (BCa) bootstrap interval
  return(p.value)
}


# Bootstrapping approach for logistic generalized linear models on unconditional
# quantiles of the response variable

logistic.Phylo.Env <- function(dataset = dataset,
                               community.phylo.variable, 
                               X,
                               quantiles = c(0.7, 0.9), 
                               n.boot = 200,
                               boot.size = 500) {
  X = dataset[, X]
  
  community.phylo.variable = dataset[, community.phylo.variable]
  
  n.quantiles <- length(quantiles)
  n <- nrow(X)
  
  phylo.condition <- matrix(0, n, 1)
  phylo.quant <- quantile(community.phylo.variable, quantiles)
  boot.list <- list()
  model.rob.list <- list()
  
  for (j in 1:n.quantiles) {
    phylo.condition[which(community.phylo.variable <= phylo.quant[j])] <- 0 # no need to set this (all values already set as zero), it's just for completion
    phylo.condition[which(community.phylo.variable > phylo.quant[j])] <- 1
    data.tmp <- data.frame(phylo.condition, scale(X))
    model.rob.list[[j]] <- glmRob(phylo.condition ~ ., 
                                  family = binomial(link = "logit"), 
                                  data = data.tmp)
    boot.coefs <- matrix(NA, 
                         nrow = n.boot, ncol = length(coef(model.rob.list[[j]])), 
                         dimnames = list(rep = seq(n.boot), 
                                         coef = names(coef(model.rob.list[[j]]))))
    for (i in 1:n.boot) {
      print(i)
      boot.data <- data.tmp[sample(nrow(data.tmp), 
                                   size = boot.size,
                                   replace = TRUE), ]
      boot.fit <- update(model.rob.list[[j]], 
                         data = boot.data) # refit the robust model with resampled data
      boot.coefs[i, ] <- coefficients(boot.fit)
    }
    boot.list[[j]] <- boot.coefs
  }
  
  p.values <- matrix(0, 
                     n.quantiles, 
                     ncol(boot.coefs))
  colnames(p.values) <- colnames(boot.coefs)
  rownames(p.values) <- paste0("q",
                               quantiles)
  for (j in 1:n.quantiles) {
    for (i in 1:ncol(boot.coefs)) {
      p.values[j, i] <- boot.pvalue.CI_inversion(boot.list[[j]][, i], 
                                                 theta_null = 0)
    }
  }
  result <- list(model.rob = model.rob.list, 
                 boot.coefs = boot.list, 
                 p.values = p.values)
  return(result)
}