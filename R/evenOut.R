#' look for inconsistency in curve fit by age
#'
#' @param AgeAlin age aligned ring widths
#'
#' @return adjusted age aligned
#' @export
#' @description
#' There should probably be a citation here, but I didn't write it down when I first wrote this code
#' 
#'
evenOut <- function(AgeAlin){

  OntoMed <- apply(AgeAlin, 1, function(x) median(x, na.rm = TRUE))
  
  #Compare MSE of points below with points above
  RMSEabove <- 0
  RMSEbelow <- 0
  RMSEdataA <- naFrame(AgeAlin)
  RMSEdataB <- naFrame(AgeAlin)
  for (i in 1:dim(AgeAlin)[2]){
    commonRows <- !is.na(AgeAlin[,i]) & !is.na(OntoMed)
    series1 <- AgeAlin[,i]
    rows1 <- as.numeric(rownames(AgeAlin)[commonRows])
    for (j in rows1){
      if (series1[j] > OntoMed[j]){
        RMSEdataA[j,i] <- sqrt((series1[j] - OntoMed[j])^2)
        RMSEabove <- RMSEabove + RMSEdataA[j,i]
      }else if (series1[j] < OntoMed[j]){
        RMSEdataB[j,i] <- sqrt((OntoMed[j] - series1[j])^2)
        RMSEbelow <- RMSEbelow + RMSEdataB[j,i]
      }
    }
  }
  print(paste0("RMSE of points above the median: ", round(RMSEabove, 0)))
  print(paste0("RMSE of points below the median: ", round(RMSEbelow, 0)))
  
  #Look at RMSE above/below by ontogenetic age
  RMSEdataA <- RMSEdataA[apply(RMSEdataA, 1, function(x) sum(!is.na(x))) > 0,]
  RMSEdataB <- RMSEdataB[apply(RMSEdataB, 1, function(x) sum(!is.na(x))) > 0,]
  ontoRMSEa <- rowMeans(RMSEdataA, na.rm = TRUE)
  ontoRMSEb <- rowMeans(RMSEdataB, na.rm = TRUE)
  x1 <- as.numeric(rownames(RMSEdataB))
  plot(x1, ontoRMSEb, type = "l", main = "Average distance to Median by Ontogenetic Age", 
       sub = "Blue = Points above median, Black = Points below", xlab = "Ontogenetic Age",
       ylab = "Average offset from median (RMSE, microns)")
  lines(x1, ontoRMSEa, col = "blue")
  
  #Adjust points to create equal spread avove and below
  ontoAspline <- detrend.series(ontoRMSEa, return.info = TRUE, method = "Spline", nyrs = 21)
  ontoBspline <- detrend.series(ontoRMSEb, return.info = TRUE, method = "Spline", nyrs = 21)
  
  aboveLoc <- !is.na(RMSEdataA)
  newAgeAligned <- AgeAlin
  oldAgeAligned <- data.frame(AgeAlin)
  for (i in 1:dim(AgeAlin)[2]){
    nowLoc <- aboveLoc[,i]
    for (j in 1:length(nowLoc)){
      if (nowLoc[j]){
        newAgeAligned[j,i] <- AgeAlin[j,i] * ontoBspline$curves[[j]]/ontoAspline$curves[[j]]
      }
    }
  }

  return(newAgeAligned)
  
}