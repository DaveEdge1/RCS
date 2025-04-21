#Age-varying spline

#' Age-varying spline
#'
#' @param series time series
#' @param SpLenRange min/max spline length in nyrs eg (3,80)
#' @param goStiff at what age should the spline length max out? or NULL
#'
#' @importFrom dplR detrend.series
#'
#' @return original series and spline
#' @export
#'
tvSpline <- function(series, SpLenRange = c(3,100), goStiff = NULL){
  library(dplR)

  if(is.null(goStiff)){
    goStiff <- length(series)
  }

  length(series)
  spLens <- seq(SpLenRange[1],SpLenRange[2],length.out=goStiff)

  keepVal <- rep(NA,length(spLens))
  for(i in 1:goStiff){
    dCurve <- detrend.series(series, method = "Spline", nyrs = spLens[i], return.info = T, make.plot = F)
    keepVal[i] <- dCurve$curves[i]
  }

  par(mfrow=c(1,1))
  plot(series, type = "l")
  lines(keepVal, col="blue")

  returns <- list("series" = series,
                  "spline" = keepVal)

  return(returns)
}
