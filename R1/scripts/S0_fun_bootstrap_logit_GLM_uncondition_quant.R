######################################################################################
##  Functions to describe and inferentially test whether how changes in predictive ###
##  variables affect community phylogenetic structure                              ###

# Code Author: Pedro Peres-Neto and Pedro Henrique Pereira Braga

# Code Last Updated: 2022-11-22

# We used two complementary approaches (one descriptive and the other inferential) to
# assess the predictions that phylogenetically clustered communities are more frequent
# in historically climatically stable regions (H2) and that increased in situ
# diversification rates generate regional clusters of closely related species (H3).

# To describe how phylogenetic relatedness of bat communities (for both NRI and NTI)
# changed as a function of historical change in temperature, historical change in
# precipitation and in situ net diversification rates, we plotted the mean
# phylogenetic relatedness of bat communities (for both NRI and NTI) across each
# percentile (100 quantiles) of the predictor variables of interest (i.e., historical
# change in temperature, historical change in precipitation and in situ net
# diversification rates; see Figures 3 and 4). This representation allowed us to
# describe how phylogenetic structure varies as a response to the predictors of
# interest.

# The function quantile.XY() below takes a dataset (data.frame) containing response
# and predictive variables (which are specified as names within the X.var and Y.var
# arguments), and a column named SamplingPool specifying the geographical scale used
# to compute Y.var (here, community phylogenetic structure). It will then represent
# the mean of Y.var across the n.classes (numerical) conditional quantiles, for each
# X.var, at the geographical scale you define (SamplingPool.i).

# The output of quantile.XY() returns a data.frame containing multiple summary
# statistics, which can be represented in a figure.

# To test hypotheses H2 and H3 inferentially, we explicitly tested how changes in
# historical climatic stability and in situ diversification rates independently
# increased (or decreased) the likelihood of a community being composed of highly
# phylogenetically related species. For this, we used two upper conditional
# proportional percentiles of community phylogenetic relatedness, the 90th and the
# 75th percentiles, as thresholds to consider whether a community was highly
# phylogenetically related. These values were chosen somewhat arbitrarily but focusing
# on these percentiles allowed us to focus on the effects of historical climatic
# stability and local net diversification affecting long-term persistence of species,
# and thus influencing highly phylogenetically related communities. We started by
# attributing a value of one for a community if its phylogenetic structure (either NRI
# or NTI) was greater than its upper conditional quantile (either the 90th or the 75th
# percentile) across all communities. If smaller, a value of zero was assigned
# instead. Based on this classification (of whether the community was highly
# phylogenetically related or not), we then applied conditionally unbiased bounded
# influence robust logistic regressions (i.e., robust to reduce the influence of
# potential outliers; see Kunsch et al., 1989) in which the response variable was the
# vector of binary outcomes (ones and zeros) representing relatively high or low
# phylogenetic structure (separately for NRI and NTI) and the predictors were z-score
# standardized values (i.e., to have mean zero and variance one) of historical change
# in temperature, historical change in precipitation and in situ net diversification
# rates. As such, we were able to estimate the relative importance of each predictor.
# To estimate confidence intervals for each predictor, we used a bootstrap approach
# based on 1,000 resamples of 2,500 random communities each. These logistic
# regressions allowed us to describe how changes in historical climatic stability and
# in situ diversification rates independently affected the log-odds (the logistic
# response) of a community being highly phylogenetically related.

# This second approach is done by the logistic.Phylo.Env() function (dependent on the
# boot.pvalue.CI_inversion() function, specified below).

# logistic.Phylo.Env() takes the same dataset (data.frame) containing response and
# predictive variables. In this case, you must have already filtered the data for the
# desired sampling pool.

# The argument community.phylo.variable is the name of the response variable ("nri" or
# "nti") within the dataset. 
# X are is a vector containing the names of the predictive variable;
# quantiles is a vector containing the quantiles of interest. In this study, we used
# 0.75 and 0.90.
# n.boot is the number of bootstraps;
# and, boot.size is the size of the bootstrap.

# The output is contains model, bootstrap coefficients, and p-values distributed in a
# list for each quantile used in the computation.

######################################################################################

# Unconditional quantile summary statistics

quantile.XY <- function(dataset = dataset,
                        X.var, 
                        Y.var,
                        SamplingPool.i,
                        n.classes = 100 #, 
                        # lab_x = "X", 
                        # lab_y = "Y"
                        ) {
  
  
  X = dataset %>%
    filter(SamplingPool == SamplingPool.i) %>%
    drop_na(Y.var) %>%
    drop_na(X.var)  %>% 
    as.data.frame() %>% 
    pull(X.var)
  
  Y = dataset %>%
    filter(SamplingPool == SamplingPool.i) %>%
    drop_na(Y.var) %>%
    drop_na(X.var) %>% 
    as.data.frame() %>% 
    pull(Y.var)
  
  # divide X into n.classes quantiles and then calculate the variance of Y in
  # each interval.
  
  quantile.values <- cumsum(rep(1 / (n.classes + 1), n.classes))
  
  # quantile.values <- 1:100/100
  
  X.quant <- stats::quantile(X, quantile.values, type = 7) # ; length(unique(X.quant))
  
  if(length(unique(X.quant)) != n.classes){
    X.quant <- X.quant[!duplicated(X.quant)]

    print("Quantiles cover non-unique values. Keeping the first duplicated occurrence and removing the following duplicates.")
  }
  
    matrix.summaries <- matrix(0, length(X.quant), 13)
    
    colnames(matrix.summaries) <- c("median", "mean", "variance", 
                                  "sd", "se", "IQR", 
                                  "min", "max",
                                  "min.X", "max.X", "mean.X", "sd.X", "se.X")
  for (i in 1:length(X.quant)) {
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
    if (i == length(X.quant)) {
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
    if ((i > 1) & (i < length(X.quant))) {
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
  
  data.df <- data.frame(
    X.quantiles = X.quant, 
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
    se.X = matrix.summaries[, "se.X"],
    x.lab = X.var,
    y.lab = Y.var,
    SamplingPool = SamplingPool.i
  )
  
  data.df$pct <- row.names(data.df)
  
  # Represent with ggplot() and geom_point()
  # 
  #   (fig.subset.mean.var_y.var_x.quantile <- ggplot(
  #     data.df,
  #     aes(
  #       x = mean.X,
  #       y = mean.Y
  #     )
  #   ) +
  #       # geom_errorbar(aes(ymin = mean.Y - 1.96*se.Y,
  #       #                   ymax = mean.Y + 1.96*se.Y),
  #       #               cex = 0.1) +
  #       # geom_errorbarh(aes(xmin = min.X,
  #       #                    xmax = max.X),
  #       #                cex = 1) +
  #       # scico::scale_fill_scico_d(
  #       #   palette = "acton",
  #       #   direction = -1,
  #       #   drop = FALSE
  #       # ) +
  #     geom_point(cex = 2.1) +
  #       labs(
  #         x = lab_x,
  #         y = lab_y
  #       ) +
  #       geom_hline(yintercept = 0, alpha = 0.25) +
  #       geom_vline(xintercept = 0, alpha = 0.25) +
  #       scale_y_continuous(
  #         breaks = pretty(c(data.df$mean.Y,
  #                           0), n = 6),
  #         limits = c(-0.05, 0.05) + range(pretty(c(data.df$mean.Y, 0), n = 7)),
  #         #   expand = expansion(mult = c(0, 0)),
  #         position = "left"
  #       ) +
  #       scale_x_continuous(
  #         breaks = pretty(c(X, 0), n = 5),
  #         limits = c(-0.05, 0.05) + range(pretty(c(data.df$mean.X, 0),
  #                                                n = 6)),
  #         #        expand = expansion(mult = c(0, 0)),
  #       ) +
  #       theme_classic(base_size = 17 * 1.6) +
  #       theme(
  #         legend.position = "none",
  #         legend.direction = "horizontal",
  #         legend.key.size = unit(0.5, "cm"),
  #         legend.text = element_blank(),
  #         legend.title = element_blank(),
  #         legend.box = "horizontal",
  #         # axis.title.x = element_blank(),
  #         axis.text.x = element_text(size = 17 * 1.6),
  #         axis.text.y = element_text(size = 17 * 1.6),
  #         axis.title.y = element_text(
  #           size = 17 * 2,
  #           face = "bold"
  #         ),
  #         axis.title.x = element_text(
  #           size = 17 * 1.4,
  #           face = "bold"
  #         ),
  #         plot.margin = unit(c(0.5, 0.6, 0, 0.6), "cm")
  #       ) +
  #       guides(colour = guide_legend(nrow = 1))
  #   )
  
  return(
    #fig.subset.mean.var_y.var_x.quantile,
    data.df
    )
  
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


###

logistic.Phylo.Env.threhsold <- function(dataset = dataset,
                                         community.phylo.variable = "nri", 
                                         X,
                                         threshold = c(1.96), 
                                         n.boot = 200,
                                         boot.size = 500) {
  X = dataset[, X]
  
  community.phylo.variable = dataset[, community.phylo.variable]
  
  n.threshold <- length(threshold)
  
  n <- nrow(X)
  
  phylo.condition <- matrix(0, n, 1)
  phylo.threshold <- setNames(threshold, threshold)
  
  boot.list <- list()
  
  model.rob.list <- list()
  
  for (j in 1:n.threshold) {
    phylo.condition[which(community.phylo.variable <= phylo.threshold[j])] <- 0 # no need to set this (all values already set as zero), it's just for completion
    phylo.condition[which(community.phylo.variable > phylo.threshold[j])] <- 1
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
                     n.threshold, 
                     ncol(boot.coefs)
                     )
  
  colnames(p.values) <- colnames(boot.coefs)
  rownames(p.values) <- paste0("t",
                               threshold)
  for (j in 1:n.threshold) {
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

