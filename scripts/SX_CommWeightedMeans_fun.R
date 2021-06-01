#######################################################################################################
###### Utility function to calculate community weighted means of trait data #####
#######################################################################################################

# Code Author: Pedro Henrique Pereira Braga
# Based on: 
# Last Update: 2020-03-10
#
# This script provides the calculation for community weighted means (CWM;
# Lavorel et al. 2008) of a trait for each ecological community. Three measures
# are calculated for each set of communities: CWM, which refers to the

# Arguments :
# "Trait" = Trait values for each OTU.
# "Distrib" = Presence-absence community data. Object must be from a class
# "data.frame", "matrix" or "vector"
#
#
# Returns : The function returns a data.frame with the community weighted means
#
# The row names of this data.frame correspond to the cell IDs of the original Distrib data.frame
#
#######################################################################################################


CWM_Std_TW <- function(Trait = Trait, Distrib = CommunityDataset){
  
  # Distrib = Comm
  # Trait = Amphibians.NicheSep[ , "Tol", drop = FALSE]
  
  Trait <- as.matrix(Trait)
  Distrib <- as.matrix(Distrib)
  
  Distrib <- Distrib[ , (colnames(Distrib) %in% row.names(Trait))]
  
  centering_mat <- function(X, w){ 
    X - rep(1, length(w)) %*% t(w) %*% X 
    }
  
  standardize_w <- function(X, w){
    ones <- rep(1, length(w))
    Xc <- X - ones %*% t(w) %*% X # equations 6 and 7
    Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w)) 
  } 
  
  # Check community matrix for communities with zero species and 
  # species with zero occurrences
  
  rows <- seq_len(nrow(Distrib))
  cols <- seq_len(ncol(Distrib))
  rni <- which(rowSums(Distrib) == 0)
  rni.nozeros.id <- which(rowSums(Distrib) != 0) # keep the ID of rows that have sum equal to zero
  
  repeat {
    if (length(rni)) {
      Distrib <- Distrib[-rni, , drop = FALSE]; rows <- rows[-rni]
      }
    ksi <- which(colSums(Distrib) == 0)
    if (length(ksi)) {
      Distrib <- Distrib[,-ksi, drop = FALSE]; cols <- cols[-ksi]
      }
    rni <- which(rowSums(Distrib) == 0)
    if (length(rni)==0 & length(ksi) == 0) { 
      break 
    }
  }
  
  # E <-E[rows, ,drop = FALSE]
  Trait <- Trait[cols, , drop = FALSE]
  
  # end check community matrix
  
  # Calculate proportional occurrences
  Distrib <- Distrib/sum(Distrib)
  
  # Estimate dimensions for the data
  
  S <- ncol(Distrib) # number of species
  n <- nrow(Distrib) # number of communities
  # p <- ncol(E) # number of environmental predictors
  q <- ncol(Trait) # number of traits
  
  # Setting up matrices
  Wn <- rowSums(Distrib) # Community Weights
  Ws <- colSums(Distrib) # Species Weights
  
  # cor matrices are trait by environment
  
  CWM <- Distrib %*% Trait / Wn  # traits weighted by the number of occurrences (or species abundances) in each community
  
  # CWM.cor <- cor(CWM,E)
  # SNC <- t(Distrib)%*%E/Ws  # weighted means wrt to E
  # SNC.cor <- cor(Trait,SNC)
  
  CWM_std_w  <- standardize_w(CWM, Wn) # CWM are weighted by the number of occurrences (or species abundances) in each community
  
  # Estd_w <- standardize_w(E,Wn)
  # wCWM.cor <- t(t(Estd_w)%*%(CWMstd_w*Wn))
  
  # SNCstd_w <- standardize_w(SNC,Ws)
  
  Tstd_w <-  standardize_w(Trait, Ws) # weighted standardized traits
  # wSNC.cor <- t(Tstd_w)%*%(SNCstd_w*Ws)
  
  # Fourth corner calculated as W_n weighted covariance between 
  # CWM and standardized Trait (trait)
  
  CWM_std_tw <- Distrib%*%Tstd_w/Wn # CWM calculated on weighted standardized traits
  
  # Variation within and among ommunities
  # Among communities 
  among.Comm.var <- sum(diag(t(CWM_std_tw)%*%(CWM_std_tw* Wn)))
  # Within communities 
  within.Comm.var <- 1 - among.Comm.var
  
  return(
    CommWeightMeans <- data.frame(row.names = rni.nozeros.id,
                                  CWM = CWM,
                                  CWM_std_tw = CWM_std_tw,
                                  CWM_std_w = CWM_std_w))
}