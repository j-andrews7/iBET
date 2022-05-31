#' @importFrom SummarizedExperiment rowData rowData<-
.swap_rownames <- function(object, swap.rownames = NULL) {

  if (!identical(swap.rownames, NULL) && is(object, "SummarizedExperiment")) {

    if (!swap.rownames %in% names(SummarizedExperiment::rowData(object))) {
      stop("'swap.rownames' is not a column of 'rowData(object)'")
    }
    rowData(object)$ORIGINAL_ROWS <- rownames(object)
    rownames(object) <- rowData(object)[,swap.rownames]
  }

  if (!identical(swap.rownames, NULL) && is(object, "data.frame")) {

    if (!swap.rownames %in% names(object)) {
      stop("'swap.rownames' is not a column of 'object'")
    }

    object$ORIGINAL_ROWS <- rownames(object)
    rownames(object) <- make.names(object[[swap.rownames]], unique = TRUE)
  }

  object
}
