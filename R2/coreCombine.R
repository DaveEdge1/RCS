#' Combine multiple ring width measurements from the same sample
#'
#' @param rwi
#'
#' @return rwi
#' @export
#'
combCores <- function(rwi, po) {

  if (ncol(po) != 2)
    stop("PO must have 2 columns: ID and PO")

  colnames(po) <- c("ID", "PO")

  # Remove exactly ONE trailing letter (upper or lowercase)
  coreID <- sub("([A-Za-z])$", "", po$ID)

  uniqueCores <- unique(coreID)

  newRWL <- list()
  newPO  <- list()

  # Labels for multi-series cores
  letterLabels <- c(LETTERS, letters)

  for (core in uniqueCores) {

    # indices of all series belonging to this core
    idx <- which(coreID == core)
    seriesNames <- po$ID[idx]
    nSeries <- length(idx)

    for (i in seq_len(nSeries)) {

      # Naming rules:
      label <- if (nSeries == 1) {
        core                     # single series → just TND123
      } else {
        paste0(core, "_", letterLabels[i])   # multi-series → TND123_A, TND123_B, …
      }

      # Extract series data
      if (!(seriesNames[i] %in% colnames(rwi))) {
        stop(paste("Series", seriesNames[i], "not found in RWI file."))
      }

      newRWL[[label]] <- rwi[[seriesNames[i]]]
      newPO[[label]]  <- po$PO[idx[i]]
    }
  }

  newRWL <- as.data.frame(newRWL)
  newPO <- data.frame(
    ID = names(newPO),
    PO = as.numeric(newPO),
    stringsAsFactors = FALSE
  )

  return(list(rwl = newRWL, PO = newPO))
}
