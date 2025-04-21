#' regional curve standardization with plotting and custom RC input
#'
#' @param rwl rwl file as from dplR
#' @param po po file as from dplR
#' @param rc regional curve formatted as the RC object returned from robustRC()
#'
#' @return series info, rwi, plot data, plot
#' @export
#'
#' @import dplR
#' @import ggplot2
#' @import RColorBrewer
#'
dRCS <- function(rwl, po, rc){


  naFrame <- function(df){
    df2 <- data.frame(matrix(data = NA, nrow = nrow(df), ncol = ncol(df)))
    colnames(df2) <- colnames(df)
    rownames(df2) <- rownames(df)
    return(df2)
  }

  #Check that the sample IDs in the PO and RWL files match
  checks <- rep(NA, length(colnames(rwl)))
  for (i in 1:length(colnames(rwl))){
    checks[i] <- colnames(rwl)[i] == po[i,1]
    if (checks[i] == FALSE){
      stop(paste0("Sample IDs in rwl file do not match those in the PO file. First unmatched ID at: ", i))}
  }

  #How many measurements were taken for each series
  seriesLengths <- apply(rwl, 2, function(x) sum(!is.na(x)))
  #Add the pith offset to get the sample longevity
  seriesLongevity <- seriesLengths + po[,2]

  #Truncate the series in the rwl based on 'ontoCut'

  #merge the RC and the given series into a single data frame based on the year,
  #aligning them ontogenetically

  #This sequence provides data frames for each series with Year, Ontogenetic Age,
  #raw values, regional curve values, and detrended series values
  seriesDat <- list()
  for (ab in 1:dim(rwl)[2]){
    seriesDat1 <- NULL
    seriesDat1 <- data.frame("Year" = as.numeric(rownames(rwl)[!is.na(rwl[,ab])]),
                                "OntoAge" = po[ab,2]:(seriesLengths[ab]+po[ab,2]-1),
                                "rawValues" = na.omit(rwl[,ab]))
    seriesDat1$RC <- rep(NA,dim(seriesDat1)[1])
    seriesDat1$RC[seriesDat1$OntoAge %in% rc$RC$BioAge] <- rc$RC$RC[rc$RC$BioAge %in% seriesDat1$OntoAge]
    seriesDat1$detValues <- seriesDat1$rawValues / seriesDat1$RC
    seriesDat[[ab]] <- seriesDat1
  }

  #Detrend the series and place them back into a data frame
  TNrwi <- naFrame(rwl)
  for (i in 1:dim(TNrwi)[2]){
    TNrwi[as.numeric(rownames(TNrwi)) %in% seriesDat[[i]][,1],i] <- seriesDat[[i]][,5]
  }

  #Get the order of series based on first measurement
  FirstRingOrder <- names(sort(apply(TNrwi,2,function(x) min(which((!is.na(x)))))))


  #Build a color palette for plotting
  pal <- colorRampPalette(c("red", "goldenrod", "darkgreen", "blue", "purple"))

  #Calculate the mean value chronology
  tnchron <- chron(TNrwi)
  plot(tnchron)
  #Smooth the chronology
  tnchron <- rollmean(tnchron[1,],11, na.pad = TRUE)

  print(paste0("Year: ", length(as.numeric(rownames(TNrwi)))))
  print(paste0("dim(rwi): ", dim(TNrwi[,FirstRingOrder])))
  print(paste0("dim(chron): ", dim(tnchron)))

  #Plot the rwi and overlay a smoothed mean chronology
  rwiDF <- cbind.data.frame("Year"=as.numeric(rownames(TNrwi)),
                            TNrwi[,FirstRingOrder],
                            "Chron" = tnchron)
  spagPlot <- melt(rwiDF, id.vars = "Year")
  p1 <- ggplot(spagPlot, aes(x=Year,y=value,color=variable, size=variable)) + geom_line() +
    theme(legend.position = "none") +
    scale_color_manual(values = c(pal(dim(TNrwi)[2]),"black")) +
    scale_size_manual(values = c(rep(0.5,dim(TNrwi)[2]),1.5))
  p1

  returns <- list("SeriesInfo" = seriesDat,
                  "rwi" = TNrwi,
                  "plotData" = spagPlot,
                  "plot" = p1)
  return(returns)
}





