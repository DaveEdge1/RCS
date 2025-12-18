#' RCS (regional curve standardization) detrending
#'
#' @param rwlFile as created by read.rwl
#' @param poFile two-column file created by dplR
#' @param ratios use ratios rather than differences (recommended)
#' @param truncRC limit curve calculation to regions of a designated minimum sample depth
#' @param rcIn two column regional curve as created from robustRC() or NULL to create new
#' @param ageMin max age of valid RC, or NULL if no truncation needed
#' @param ageMax min age of valid RC, or NULL if no truncation needed
#'
#' @import dplR
#'
#' @return rwi and assorted goodies
#' @export
#'
newRCS <- function(rwlFile = NULL, poFile = NULL, ratios = TRUE,
                   truncRC = 10, rcIn = NULL, ageMin = NULL, ageMax = NULL
                   ){
  #########################################################
  #########################################################
  #Step 1
  #Check that the sample IDs in the PO and RWL files match
  checks <- rep(NA, length(colnames(rwlFile)))
  for (i in 1:length(colnames(rwlFile))){
    checks[i] <- colnames(rwlFile)[i] == poFile[i,1]
    if (checks[i] == FALSE){
      stop(paste0("Sample IDs in rwl file do not match those in the PO file. First unmatched ID at: ", i))}
  }

  #########################################################
  #########################################################
  #Step 2
  #Add a small random numer to all values as unique identifiers
  addRand <- rnorm(length(rwlFile[!is.na(rwlFile)]), .0001, .0001)
  rwlFile[!is.na(rwlFile)] <- rwlFile[!is.na(rwlFile)] + addRand
  #########################################################
  #########################################################
  #Step 3
  #Build RC if user does not provide
  if(is.null(rcIn)){
    rcIn = robustRC(rwlFile = rwlFile, poFile = poFile)$RC
  }
  #########################################################
  #########################################################
  #Step 4
  #Truncate RC if user requests by "ageMix" and/or "ageMax"
  if (!is.null(ageMin)){
    rcIn <- rcIn[rcIn[,1] >= ageMin,]
  }

  if (!is.null(ageMax)){
    rcIn <- rcIn[rcIn[,1] <= ageMax,]
  }
  #########################################################
  #########################################################
  #Step 5
  #Align the samples based on age
  ageAligned <- c(1,2,3,4,5)

  for (i in 1:length(colnames(rwlFile))){
    naRemoved <- array(na.omit(rwlFile[,i]))
    naFill <- rep(NA, poFile[i,2])
    newSeries <- c(naFill, naRemoved)
    ageAligned <- cbind.NA(ageAligned, newSeries)
  }
  ageAligned[,1] <- NULL
  colnames(ageAligned) <- colnames(rwlFile)
  #########################################################
  #########################################################
  #Step 6
  #Truncate age-aligned data based on limits of RC
  ageLimits <- NULL
  aa <- NULL
  bb <- NULL
  warnTrig <- 0
   if (!is.null(rcIn)){
    rcIn <- apply(rcIn, 2, function(x) as.numeric(x))
    earlyNA <- 1:(min(rcIn[,1]) - 1)
    lateNA <- (max(rcIn[,1]) + 1):length(rownames(ageAligned))
    ageLimits <- c(min(rcIn[,1]), max(rcIn[,1]))
    ageAligned <- ageAligned[rcIn[,1],]
    if(suppressWarnings(as.integer(rownames(ageAligned)))[1] < min(rcIn[,1])){
      aa <- as.data.frame(matrix(nrow = length(earlyNA), ncol = length(colnames(ageAligned))))
      ageAligned <- rbind.data.frame(aa, ageAligned)
      warnTrig <- 1
    }
    if(max(suppressWarnings(as.integer(rownames(ageAligned))), na.rm = T) > max(rcIn[,1])){
      bb <- as.data.frame(matrix(nrow = length(lateNA), ncol = length(colnames(ageAligned))))
      ageAligned <- rbind.data.frame(ageAligned, bb)
      warnTrig <- 1
    }
    if (warnTrig == 1){
      warning("Data truncated based on age limits imposed from external RC input. See 'ageLimits'.")
    }
  }
  #########################################################
  #########################################################
  #Step 7
  #Calculate average growth by biological age
  avgCurve <- rowMeans(ageAligned, na.rm = TRUE)
  plot(x = suppressWarnings(as.integer(rownames(ageAligned))), y = avgCurve, type = "l", lwd = 2,
       main = "Average Growth by Ontogentetic Age vs. RC", ylab = "Growth (microns)", xlab = "Ontogenetic Age")
  lines(x = rcIn[,1], y = rcIn[,2], col = "blue")
  #########################################################
  #########################################################
  #Step 8
  #Calculate the sample depth of the age-aligned data and truncate the data based
  #on required sample depth - "truncRC"
  RCsampleDepth <- rep(NA, length(rownames(ageAligned)))
  for (i in 1:length(rownames(ageAligned))){
    RCsampleDepth[i] <- sum(!is.na(ageAligned[i,]))
  }
  rcsInfo <- cbind.data.frame(avgCurve, as.numeric(RCsampleDepth))
  rcsInfo <- rcsInfo[complete.cases(rcsInfo),]
  rcsInfo <- rcsInfo[rcsInfo[,1] >= truncRC,]
  rcsInfo$BioAge <- as.integer(rownames(rcsInfo))
  colnames(rcsInfo)[2] <- "SampDepth"

  ageAligned <- ageAligned[suppressWarnings(as.numeric(rownames(ageAligned))) >= min(rcsInfo$BioAge),]
  ageAligned <- ageAligned[suppressWarnings(as.numeric(rownames(ageAligned))) <= max(rcsInfo$BioAge),]
  #########################################################
  #########################################################
  #Step 9
  #Grab the RC for the period corresponding to the new truncated age-aligned data
  curveSub <- c(min(rcsInfo$BioAge), max(rcsInfo$BioAge))
  rcsInfo$RCurve <- rcIn[which(rcIn[,1] >= curveSub[1] & rcIn[,1] <= curveSub[2]),2]
  #########################################################
  #########################################################
  #Step 10
  #Find chronology year positions of age-aligned data before detrending
  seriesLengths <- apply(ageAligned, 2, function(x) sum(!is.na(x)))
  rwlPos <- data.frame(matrix(ncol = ncol(rwlFile), nrow = 2, data = NA))
  colnames(rwlPos) <- colnames(rwlFile)
  solveIt <- 0
  for (i in 1:length(seriesLengths)){
    firstNum <- min(which(!is.na(ageAligned[,i])))
    lastNum <- firstNum + seriesLengths[i] - 1
    minIndex <- which(ageAligned[firstNum,i] == rwlFile[,i])
    maxIndex <- which(ageAligned[lastNum,i] == rwlFile[,i])
    rwlPos[,i] <- c(minIndex,maxIndex)
  }
  #########################################################
  #########################################################
  #Step 11
  #Detrend age-aligned data
  ageAligned <- ageAligned[apply(ageAligned, 1, function(x) sum(!is.na(x))) != 0,]

  if (sum(na.omit(suppressWarnings(as.numeric(rownames(ageAligned)))) != rcsInfo$BioAge, na.rm = TRUE) != 0){
    stop("RC and age-aligned data do not match!")
  }
  #cat(rcsInfo$RCurve, "\n")
  newAgeAligned <-  cbind.data.frame(ageAligned, rcsInfo$RCurve)
  ageAlignedDet <- data.frame(matrix(ncol = ncol(ageAligned), nrow = nrow(ageAligned), data = NA))
  colnames(ageAlignedDet) <- colnames(ageAligned)
  rownames(ageAlignedDet) <- rownames(ageAligned)
  RCspot <- dim(newAgeAligned)[2]
  for (i in 1:dim(ageAligned)[2]){
    p1 <- min(which(!is.na(ageAligned[,i])))
    p2 <- p1 + seriesLengths[i] - 1
    #EDITED LINES***
    num <- newAgeAligned[p1:p2, i]
    den <- newAgeAligned[p1:p2, RCspot]
    if(ratios){
      ageAlignedDet[p1:p2,i] <- num / den
    } else {
      ageAlignedDet[p1:p2,i] <- num - den
    }
    #***
  }
  #########################################################
  #########################################################
  #Step 12
  #Place detrended data in chronology positions
  rwi <- data.frame(matrix(ncol = ncol(rwlFile), nrow = nrow(rwlFile), data = NA))
  colnames(rwi) <- colnames(rwlFile)
  rownames(rwi) <- rownames(rwlFile)
  for (i in 1:dim(rwi)[2]){
    getIndices <- rwlPos[1,i]:rwlPos[2,i]
    rwi[getIndices,i] <- na.omit(ageAlignedDet[i])
  }
  #########################################################
  #########################################################
  #Step 13
  #Plot the chronology
  #EDITED LINES***
  crn_obj <- dplR::chron(rwi)
  # Use generic plot on cron object
  plot(crn_obj,
       main = "Detrended Chronology",
       xlab = "Year",
       ylab = "Index")
  #***

  return(list("rwi" = rwi, "rcInfo" = rcsInfo, "ageAligned" = ageAligned, "ageAlignedDetrended" = ageAlignedDet,
              "ageLimits" = ageLimits))
}





