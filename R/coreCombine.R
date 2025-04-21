#' Combine multiple ring width measurements from the same sample
#'
#' @param rwi
#'
#' @return rwi
#' @export
#'
combCores <- function(rwi, PO=NULL){

  if(!is.null(PO)){
    checks <- rep(NA, length(colnames(rwi)))
    for (i in 1:length(colnames(rwi))){
      checks[i] <- colnames(rwi)[i] == PO[i,1]
      if (checks[i] == FALSE){
        stop(paste0("Sample IDs in rwl file do not match those in the PO file. First unmatched ID at: ", i))
      }
    }
  }

  source("D:/Lab Backup/R/is_odd.R")
  newIDS <- data.frame(matrix(ncol = 2, nrow = dim(rwi)[2]))
  colnames(newIDS) <- c("tree", "core")
  #Get rwi stats
  for (i in 1:dim(rwi)[2]){
    #step through i=1 and i=2 setting the "tree" ID to 001 and the "core" ID to 1 and 2
    newIDS[i,"tree"] <- ceiling(i/2)
    newIDS[i,"core"] <- colnames(rwi)[i]

  }

  oneShell <- sapply(seq(1,ncol(rwi),2), function(i) {
    rowMeans(rwi[,c(i, i+1)], na.rm=T)
  })

  colnames(oneShell) <- unique(gsub('.{1}$', '', newIDS$core))


  if(!is.null(PO)){
    n <- length(PO[,2])
    group <- gl(n, 2, n)
    onePO <- aggregate(PO[,2] ~ group, FUN = min)

    PO <- data.frame("ID" = unique(gsub('.{1}$', '', newIDS$core)),
                     "PO" = onePO[,2])
  }

  returns <- list("SampleIDs" = newIDS,
                  "rwl" = oneShell,
                  "PO" = PO)

  return(returns)
}





