#Example: if the minimum ontogenetic cutoff is 5 and the max is 100 (as in ontoCut = c(5,100))
#A series with a pith offset of 4, 100 measurements in length (measurements of ontogenetic age
#4 through 103), would have its first measurement and its final 3 measurements removed to produce
#a truncated series extending from ontogenetic age 5 to 100, 96 increments in length
#' truncate the series in the rwl file with min/max ontogenetic ages
#'
#' @param rwl as produced by dplR
#' @param po as produced by dplR
#' @param ontoCut c(x,y)
#'
#' @return rwl, po
#' @export
#'
ontoTrunc <- function(rwl, po, ontoCut){
  
  #Check that the sample IDs in the PO and RWL files match
  checks <- rep(NA, length(colnames(rwl)))
  for (i in 1:length(colnames(rwl))){
    checks[i] <- colnames(rwl)[i] == po[i,1]
    if (checks[i] == FALSE){
      stop(paste0("Sample IDs in rwl file do not match those in the PO file. First unmatched ID at: ", i))}
  }
  
  goHalf <- FALSE
  goFull <- FALSE
  #Check ontoCut input
  if (length(ontoCut) == 1 & is.numeric(ontoCut)){
    goHalf <- TRUE
  }else if (length(ontoCut) == 2 & is.numeric(ontoCut[1]) & is.numeric(ontoCut[2])){
    goFull <- TRUE
    goHalf <- TRUE
  }else{
    stop("ontoCut must be 1 or 2 integers as in ontoCut = 3 or ontoCut = c(3,80)")
  }
  
  if (goFull){
    for (i in 1:dim(po)[1]){
      seriesLen <- sum(!is.na(rwl[,i]))
      maxOnto <- po[i,2] + seriesLen - 1
      if (maxOnto > ontoCut[2]){
        diff1 <- maxOnto - ontoCut[2]
        cut <- 1:diff1
        cut <- cut + max(which(!is.na(rwl[,i]))) - diff1
        rwl[cut,i] <- NA
      }
    }
  }
  
  
  #Check po file for series whose first increment is less than the early age cutoff
  #Iterate over series and remove any info prior to the first ontogentic cutoff
  
  if (goHalf){
    for (i in 1:dim(po)[1]){
      if (po[i,2] < ontoCut[1]){
        diff1 <- ontoCut[1] - po[i,2]
        cut <- 1:diff1
        cut <- cut + min(which(!is.na(rwl[,i]))) - 1
        rwl[cut,i] <- NA
        po[i,2] <- ontoCut[1]
      }
    }
  }
  
  returns <- list("rwl" = rwl, "po" = po)
  return(returns)
  
  
}