#' Shift regions upstream or downstream
#'
#' @param gr regions as GenomicRanges
#' @param direction upstream or downstream
#' @param window size of shift
#' @export
#' 
shift <- function(gr, direction, window) {
  plus <- gr[strand(gr) == "+"]
  minus <- gr[strand(gr) == "-"]
  
  if(window < 0) 
    stop("window should not be negative!")
  
  if(direction == "upstream")
    window <- -window
  else if(direction != "downstream")
    stop("not a valid direction: ", direction)
  
  c(GenomicRanges::shift(plus, window), GenomicRanges::shift(minus, -window))
}

#' Extend regions upstream or downstream
#'
#' @param gr regions as GenomicRanges
#' @param direction upstream or downstream
#' @param window size of extension
#' @export
#' 
extend <- function(gr, direction, window) {
  plus <- gr[strand(gr) == "+"]
  minus <- gr[strand(gr) == "-"]
  
  if(window < 0) 
    stop("window should not be negative!")
  
  if(direction == "upstream") {
    start(plus) <- start(plus) - window
    end(minus) <- end(minus) + window
  } else if(direction == "downstream") {
    end(plus) <- end(plus) + window
    start(minus) <- start(minus) - window
  } else {
    stop("not a valid direction: ", direction)
  }
  
  c(plus, minus)
  
}

#' Inject a complete sample from shifted ods into original ods
#' 
#' @param ods OutriderDataSet containing real counts
#' @param shifted RangedSummarizedExperiment on shifted features
#' @param sampleID ID of sample (patient)
#' @return simulated OutriderDataSet containing all simulated outliers for 1 sample
#' @export
#' 
injectSample <- function(ods, shifted, sampleID) {
  
  simulated <- copyODS(ods)
  
  f <- which(names(shifted) %in% names(simulated))
  s <- which(colnames(simulated) == sampleID)
  SummarizedExperiment::assay(simulated)[, s] <- SummarizedExperiment::assay(shifted)[f, s]
  
  # remove zero counts
  simulated <- simulated[rowMeans2(counts(simulated)) > 0]
  OUTRIDER:::checkCountRequirements(simulated)
  
  # add injection coordinates
  simulated@metadata$simulated <- data.table(row = 1:nrow(simulated), col = s)
  
  simulated
}

#' Take random features and replace these in simulated
#' 
#' @param ods OutriderDataSet containing real counts
#' @param shifted RangedSummarizedExperiment on shifted features
#' @param n number of features to randomly choose
#' @param prob probability of an outlier
#' @param seed for sampling the positions of injection
#' @return simulated OutriderDataSet containing outliers at random positions and INDICES of injections
#' @export
#' 
injectFeatures <- function(ods, shifted, prob = 1E-4, seed = NULL) {
  
  # match dimensions
  shifted_matched <- shifted[rownames(ods)]
  
  shiftedCounts <- SummarizedExperiment::assay(shifted_matched)

  # get random positions
  if (!is.null(seed)) set.seed(seed)
  inject_vector <- sample(c(T,F),
                          nrow(ods)*ncol(ods),
                          prob = c(prob, 1-prob), 
                          replace = T)
  inject_matrix <- matrix(inject_vector, nrow = nrow(ods))
  
  # inject counts
  simulated <- copyODS(ods)
  simul_counts <- counts(simulated)
  simul_counts[inject_matrix] <- shiftedCounts[inject_matrix]
  
  # remove zero rows
  zeroRows <- rowMeans2(simul_counts) == 0
  simul_counts <- simul_counts[!zeroRows,]
  inject_matrix <- inject_matrix[!zeroRows,]
  simulated <- simulated[!zeroRows]
  
  # save into ODS
  counts(simulated) <- simul_counts
  simulated@metadata$simulated <- as.data.table(which(inject_matrix, arr.ind = T))
  OUTRIDER:::checkCountRequirements(simulated)
  
  simulated
}

#' Remove simulated 0-count rows
#'
#' @param ods simulated OutriderDataSet containing injections
#' @deprecated
#' 
removeZeroCounts <- function(ods) {
  
  zeroRows <- rowMeans2(counts(ods)) == 0
  
  if (sum(zeroRows) > 0) {
    # remove from index
    injections <- ods@metadata$simulated
    ods@metadata$simulated <- injections[!row %in% which(zeroRows)]
    # remove rows from count matrix
    ods <- ods[!zeroRows]
  }
  
  ods
}

#' Simulated Outliers
#' 
#' Conveniently retrieve the feature ID and sample IDs of the injected outliers
#' 
#' @param ods OutriderDataSet after injection i.e. containing the metadata entry "injections"
#' @param indices if TRUE return the indices, else translate to feature and sample IDs
#' @return return the data.table of injections
#' @export
#' 
simulatedOutliers <- function(ods, indices = F) {
  
  injections <- ods@metadata$simulated
  if (indices) injections
  else data.table(feature_id = as.integer(rownames(ods)[injections$row]),
                  sampleID = colnames(ods)[injections$col])
}
