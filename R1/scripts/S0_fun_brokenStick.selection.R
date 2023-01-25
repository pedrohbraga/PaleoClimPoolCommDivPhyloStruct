brokenStick.selection <- function(dudi.pca.object){
  
  ## First, choose how many axes are going to be kept
  ev <- dudi.pca.object$eig
  
  # Broken stick model, largely based on Pierre Legendre's class
  
  n <- length(dudi.pca.object$eig)
  bsm <- data.frame(j = seq(1:n), 
                    p = 0)
  bsm$p[1] <- 1/n
  
  for (i in 2:n){
    bsm$p[i] = bsm$p[i - 1] + (1/(n + 1 - i))
  }
  
  bsm$p <- 100*bsm$p/n
  
  ## Plot eigenvalues and % of variance for each axis
  
  plot.new()
  
  par(mfrow = c(2, 1))
  
  barplot(ev, 
          main = "Eigenvalues", 
          col = "bisque", las = 2)
  
  abline(h=mean(ev), 
         col="red")		# average eigenvalue
  
  legend("topright", 
         "Average Eigenvalue", 
         lwd = 1, 
         col = 2, 
         bty = "n")
  
  brokenStick.matrix <- t(cbind(100*ev/sum(ev), bsm$p[n:1]))
  brokenStick.kept <- length(which(brokenStick.matrix[1, ] > brokenStick.matrix[2, ]))
  
  barplot(brokenStick.matrix, 
          beside = TRUE,
          main = "% variance", 
          col = c("bisque", 2), 
          las = 2)
  
  legend("topright", 
         c("% Eigenvalue", 
           "Broken-stick model"), 
         pch= 15, 
         col= c("bisque", 2), 
         bty= "n")
  
  return(list(brokenStick.kept = brokenStick.kept,
              brokenStick.matrix = brokenStick.matrix))
}
