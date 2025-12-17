#Age-varying spline

#' Age-varying spline
#'
#' @param series time series
#' @param SpLenRange min/max spline length in nyrs eg (3,80)
#' @param goStiff at what age should the spline length max out? or NULL
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplR detrend.series
#'
#' @return original series and spline
#' @export
#'
tvSpline <- function(series, SpLenRange = c(3,100), goStiff = NULL, plot=FALSE){

  if(is.null(goStiff)){
    goStiff <- length(series)
  }

  length(series)
  spLens <- seq(SpLenRange[1],SpLenRange[2],length.out=goStiff)

  keepVal <- rep(NA,length(spLens))
  for(i in 1:goStiff){
    #EDITED LINES***
    dCurve <- dplR::detrend.series(series, method = "Spline", nyrs = spLens[i], return.info = TRUE)
    keepVal[i] <- dCurve$curve[i]

  }
  if (plot) {
    plotData <- data.frame(Age = 1:length(series),
                           Original = series,
                           Spline = keepVal)
    plotDataMelt <- tidyr::pivot_longer(plotData,
                                        cols = c("Original", "Spline"),
                                        names_to = "variable",
                                        values_to = "value")

    p1 <- ggplot2::ggplot(plotDataMelt,
                          aes(x = Age, y = value, color = variable, size = variable)) +
      geom_line() +
      scale_color_manual(values = c("Original" = "grey50", "Spline" = "blue")) +
      scale_size_manual(values = c("Original" = 0.5, "Spline" = 1.5)) +
      labs(x = "Ontogenetic Age", y = "Growth (microns)", title = "Age-varying Spline") +
      theme_minimal() +
      theme(legend.title = element_blank())

    print(p1)
  }
  return(list(series = series, spline = keepVal))

  #***

}
