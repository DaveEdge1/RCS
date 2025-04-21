#' create a regional curve (average ring with by ontogenetic age)
#'
#' @param rwlFile as created by read.rwl in dplR
#' @param poFile two-column file created by dplR
#' @param truncRC limit curve calculation to regions of a designated minimum sample depth
#' @param aAligned TRUE/FALSE, is the rwl file already aligned by ontogenetic age?
#' @param evO adjust RC fit (see evenOut())
#' @param dMethod how should we fit the curve to the ontogenetic age values? methods from dplR::detrend()
#' @param spLen if dMethod is "Spline", what spline length (in number of years)
#' @param tvSpline TRUE/FALSE, use an age-varying spline (recommended, dMethod must be set to NULL)
#' @param tvRange min and max spline lengths in nyrs (if using tvSpline), eg. c(3,80)
#' @param tvStiff set the age at which the spline length will max out (if using tvSpline)
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import zoo
#'
#' @return curve, info, and plot
#' @export
#'
robustRC <- function(rwlFile = NULL, poFile = NULL, truncRC = 20, aAligned = FALSE,
                     evO = FALSE, dMethod = NULL, spLen = NULL, tvSpline = TRUE,
                     tvRange =c(3,80), tvStiff = NULL){

  #
  # rwlFile = read.rwl("TreeNobAllLumped10-7.csv")
  # poFile = read.csv("TN_POLumped_Oct_7_2020.csv")


  #Check that the sample IDs in the PO and RWL files match

  checks <- rep(NA, length(colnames(rwlFile)))
  for (i in 1:length(colnames(rwlFile))){
    checks[i] <- colnames(rwlFile)[i] == poFile[i,1]

    if (checks[i] == FALSE){
      stop(paste0("Sample IDs in rwl file do not match those in the PO file. First unmatched ID at: ", i))
    }
  }

  #Align the samples based on age
  if (aAligned){
    ageAligned <- rwlFile
  }else{
    ageAligned <- NULL
    for (i in 1:length(colnames(rwlFile))){
      naRemoved <- array(na.omit(rwlFile[,i]))
      naFill <- rep(NA, poFile[i,2])
      naFill2 <- rep(NA, (dim(rwlFile)[1] - length(naFill) - length(naRemoved)))
      newSeries <- c(naFill, naRemoved, naFill2)
      ageAligned <- cbind(ageAligned, newSeries)
    }
    colnames(ageAligned) <- colnames(rwlFile)
    ageAligned <- data.frame(ageAligned)
  }

  if(evO){
    ageAligned <- evenOut(AgeAlin = ageAligned)
  }


  #################################################################################
  #Sort by lifespan
  seriesLengths <- apply(ageAligned, 2, function(x) sum(!is.na(x)))
  seriesLongevity <- seriesLengths + poFile[,2]
  shortSeriesAligned <- ageAligned[,seriesLongevity<70]
  #################################################################################

  #Calculate average growth by biological age
  #avgCurve <- rowMeans(ageAligned, na.rm = TRUE)
  avgCurve <- apply(ageAligned, 1, tbrm)

  #Calculate the sample depth of the age-aligned data and format the data for plotting
  RCsampleDepth <- rep(NA, length(rownames(ageAligned)))
  for (i in 1:length(rownames(ageAligned))){
    RCsampleDepth[i] <- sum(!is.na(ageAligned[i,]))
  }
  rcsInfo <- data.frame(cbind(avgCurve, RCsampleDepth))
  rcsInfo <- subset(rcsInfo, rcsInfo[,2] >= truncRC)
  rcsInfo$BioAge <- as.numeric(rownames(rcsInfo))
  rownames(rcsInfo) <- 1:length(rcsInfo$BioAge)

  #Fit the RC with a spline
  if(!is.null(dMethod)){
    detrendInfo <- detrend.series(rcsInfo$avgCurve, method = dMethod, return.info = TRUE, nyrs = spLen)
    rcsInfo$RC <- detrendInfo$curves
  }else if(is.null(dMethod)){
    detrendInfo <- tvSpline(rcsInfo$avgCurve, SpLenRange = tvRange, goStiff = tvStiff)
    rcsInfo <- cbind.NA(rcsInfo, data.frame("RCS"=detrendInfo$spline))
  }

  #plot the age aligned series and overlay the RC
  plotData <- cbind.NA(ageAligned[as.numeric(rownames(ageAligned)) %in% rcsInfo[,3],], rcsInfo[,c(3,1,4)])
  plotData <- reshape2::melt(plotData, id.vars = "BioAge")

  p1 <- ggplot2::ggplot(data = plotData, mapping = aes(x=BioAge, y=value, color=variable, size=variable)) + geom_line() +
    theme(legend.position = "none") +
    scale_color_manual(values = c(rep("grey50",dim(ageAligned)[2]), "black","darkblue")) +
    scale_size_manual(values = c(rep(0.5,dim(ageAligned)[2]),1,1.5))


  print(p1)


  return(list("RC" = rcsInfo[,3:4], "rcInfo" = rcsInfo[,1:2], "AgeAligned" = ageAligned, "plot" = p1))
}
