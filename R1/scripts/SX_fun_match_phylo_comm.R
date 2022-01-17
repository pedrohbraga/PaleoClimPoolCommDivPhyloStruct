#############################################################################################
#### Function to match (and subset) species between phylogenetic and community data sets ####
#### Modified from the picante package                                                   ####   
#############################################################################################

# The difference here is that it does not drop names.

match.phylo.comm. <- function(phy, comm, silent = TRUE) {
  if (!(is.data.frame(comm) | is.matrix(comm))) {
    stop("Community data should be a data.frame or matrix with samples in rows and taxa in columns")
  }
  res <- list()
  phytaxa <- phy$tip.label
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to match phylogeny and community data")
  }
  if (!all(commtaxa %in% phytaxa)) {
    if(silent == FALSE) {
      print("Dropping taxa from the community because they are not present in the phylogeny:")
      print(setdiff(commtaxa, phytaxa))
    }
    comm <- comm[, intersect(commtaxa, phytaxa), drop = FALSE]
    commtaxa <- colnames(comm)
  }
  if (any(!(phytaxa %in% commtaxa))) {
    if(silent == FALSE) {
      print("Dropping tips from the tree because they are not present in the community data:")
      print(setdiff(phytaxa, commtaxa))
    }
    res$phy <- prune.sample(comm, phy)
  }
  else {
    res$phy <- phy
  }
  res$comm <- comm[, res$phy$tip.label, drop = FALSE]
  return(res)
}

