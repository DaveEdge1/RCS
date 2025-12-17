#' regional curve standardization with plotting and custom RC input
#'
#' @param rwl rwl file as from dplR
#' @param po po file as from dplR
#' @param rc regional curve formatted as the RC object returned from robustRC()
#'
#' @return list containing series info, rwi, plot data, plot
#' @export
#'
#' @import dplR
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom tidyr pivot_longer
#' @importFrom zoo rollmean
dRCS <- function(rwl, po, rc){
  
  # helper: empty frame with same dims
  naFrame <- function(df){
    df2 <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
    colnames(df2) <- colnames(df)
    rownames(df2) <- rownames(df)
    df2
  }
  
  # Check sample ID consistency
  if(!all(colnames(rwl) == po[,1])){
    bad <- which(colnames(rwl) != po[,1])[1]
    stop("Sample IDs in rwl and PO do not match. First mismatch at column ", bad)
  }
  
  # Basic series stats
  seriesLengths <- apply(rwl, 2, function(x) sum(!is.na(x)))
  seriesLongevity <- seriesLengths + po[,2]
  
  # Build series-level data objects
  seriesDat <- vector("list", ncol(rwl))
  
  for(ab in seq_len(ncol(rwl))){
    # extract valid raw values
    rawVals <- rwl[!is.na(rwl[,ab]), ab]
    
    # Ontogenetic age sequence
    onto <- po[ab,2]:(po[ab,2] + length(rawVals) - 1)
    
    # Year vector
    years <- as.numeric(rownames(rwl)[!is.na(rwl[,ab])])
    
    df <- data.frame(
      Year      = years,
      OntoAge   = onto,
      RawValues = rawVals
    )
    
    # Match RC by ontogenetic age
    df$RC <- rc$RC$RC[ match(df$OntoAge, rc$RC$BioAge) ]
    
    # Detrended values
    df$detValues <- df$RawValues / df$RC
    
    seriesDat[[ab]] <- df
  }
  
  # Fill detrended series back into an rwi matrix
  TNrwi <- naFrame(rwl)
  for(i in seq_len(ncol(TNrwi))){
    idx <- match(seriesDat[[i]]$Year, as.numeric(rownames(TNrwi)))
    TNrwi[idx, i] <- seriesDat[[i]]$detValues
  }
  
  # Order by first ring year
  FirstRingOrder <- names(sort(apply(TNrwi, 2, function(x) min(which(!is.na(x))))))
  
  # Colors
  pal <- colorRampPalette(c("red", "goldenrod", "darkgreen", "blue", "purple"))
  
  # Build chronology
  ch <- dplR::chron(TNrwi, prefix = "")
  stdChron <- ch$std
  smChron <- zoo::rollmean(stdChron, 11, fill = NA)
  
  # Print for debugging
  print(paste0("Years in matrix: ", nrow(TNrwi)))
  print(paste0("dim(rwi): ", paste(dim(TNrwi[, FirstRingOrder]), collapse=" x ")))
  print(paste0("length(chron): ", length(smChron)))
  
  # Plotting DF
  rwiDF <- data.frame(
    Year = as.numeric(rownames(TNrwi)),
    TNrwi[, FirstRingOrder],
    Chron = smChron
  )
  
  spagPlot <- tidyr::pivot_longer(
    rwiDF, cols = -Year,
    names_to = "variable",
    values_to = "value"
  )
  
  # Plot
  p1 <- ggplot(spagPlot, aes(x=Year, y=value, color=variable, size=variable)) +
    geom_line() +
    theme(legend.position = "none") +
    scale_color_manual(values = c(pal(ncol(TNrwi)), "black")) +
    scale_size_manual(values = c(rep(0.5, ncol(TNrwi)), 1.5))
  
  list(
    SeriesInfo = seriesDat,
    rwi        = TNrwi,
    plotData   = spagPlot,
    plot       = p1
  )
}