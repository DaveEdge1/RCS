#' Combine multiple ring width measurements from the same sample
#'
#' @param rwi ring width data as data.frame from dplR
#' @param po pith offset data with columns for sample ID and pith offset
#'
#' @return list with combined rwl and updated PO
#' @export
#'
combCores <- function(rwi, po) {

  if (ncol(po) != 2)
    stop("PO must have 2 columns: ID and PO")

  colnames(po) <- c("ID", "PO")

  # Remove exactly ONE trailing letter (upper or lowercase) to get core ID
  coreID <- sub("([A-Za-z])$", "", po$ID)

  uniqueCores <- unique(coreID)

  newRWL <- list()
  newPO  <- list()

  for (core in uniqueCores) {

    # indices of all series belonging to this core
    idx <- which(coreID == core)
    seriesNames <- po$ID[idx]
    nSeries <- length(idx)

    # Extract all series data for this core
    seriesData <- lapply(seriesNames, function(sname) {
      if (!(sname %in% colnames(rwi))) {
        stop(paste("Series", sname, "not found in RWI file."))
      }
      rwi[[sname]]
    })

    # Combine series by averaging (handles NA values properly)
    if (nSeries == 1) {
      # Single series - just use it directly
      combinedSeries <- seriesData[[1]]
    } else {
      # Multiple series - average them row-wise, ignoring NAs
      seriesMatrix <- do.call(cbind, seriesData)
      combinedSeries <- apply(seriesMatrix, 1, function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          return(mean(x, na.rm = TRUE))
        }
      })
    }

    # Store the combined series
    newRWL[[core]] <- combinedSeries

    # Use the minimum PO value for the combined series
    newPO[[core]] <- min(po$PO[idx])
  }

  newRWL <- as.data.frame(newRWL)

  # Preserve the year rownames from the original rwl object
  rownames(newRWL) <- rownames(rwi)

  newPO <- data.frame(
    ID = names(newPO),
    PO = as.numeric(newPO),
    stringsAsFactors = FALSE
  )

  return(list(rwl = newRWL, PO = newPO))
}
