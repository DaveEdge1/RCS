#' build a blank data frame of the same size
#'
#' @param df 
#'
#' @return blank df
#' @export
#'
naFrame <- function(df){
  df2 <- data.frame(matrix(data = NA, nrow = nrow(df), ncol = ncol(df)))
  colnames(df2) <- colnames(df)
  rownames(df2) <- rownames(df)
  return(df2)
}