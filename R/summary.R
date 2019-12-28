#' Select columns and write to tsv file
#' 
#' @param results results table from pipeline output
#' @param filename specify location
#' @param type output type, either feature_wise or gene_wise
#' @export
#' 
writeResults <- function(results, filename, type = "gene_wise") {
  
  if (type == "gene_wise") {
    genes <- results[, 
                     .(gene_name, gene_id, gene_type, transcript_id, transcript_type, 
                       sampleID, EXOME_ID, FIBROBLAST_ID,
                       outrider, CANDIDATE_GENE, KNOWN_MUTATION,
                       pheno_v, biochem_v, CLINICAL_SYMPTOMS,
                       COMMENT, DISEASE, MITOCARTA, HANS_CLASS,
                       zScore_genewise, zscore_dist = max(zscore_dist)
                     ), by = gene_id]
    genes <- unique(genes)
    write.table(genes, filename, sep = "\t", quote = F, row.names = F)
  } else if (type == "feature_wise") {
    write.table(
      results[, .(gene_name, gene_id, transcript_id, type, feature_id, sampleID,
                  pValue, padjust, zScore, zScore_genewise, zscore_dist, outrider,
                  MITOCARTA, HANS_CLASS)],
      filename, sep = "\t", quote = F, row.names = F
    )
  }
  else stop(type, " is not a valid output type!")
  
}

#' Get all OUTRIDER results
#'
#' Efficiently retrieve only certain assays for all features
#' @param ods OUTRIDER data set
#' @param assays to subset which values are to be extracted. Vector can either be empty or contain only names of existing assays in the ods (names(assays(ods)))
#' @param meta metadata columns to include
#' @param id.as.int wether gene name should be coerced to integer (useful for feature_ids)
#' @export
#' 
get_all_results <- function(ods, assays = c("pValue", "zScore"), meta = c("feature_id"), id.as.int = F) {
  
  sampleIDs <- colnames(ods)
  features <- rownames(ods)
  if (id.as.int) features <- as.integer(rownames(ods))
  results <- data.table(feature_id = rep(features, times = length(samples)),
                        sampleID = rep(sampleIDs, each = length(features)))
  for (column in assays) {
    results[, tmp := c(assay(ods, column))]
    setnames(results, old = "tmp", new = column)
  }
  
  # add meta columns
  meta <- unique(c(meta, "feature_id"))
  if (length(meta) > 1)
    results <- as.data.table(merge(results, mcols(ods)[, meta]))
  
  results
}


#'
#'convenience function that adds summary statistics to given counts as data.table
#' @param se summarized experiment
#' @return a data.table containing the summarised values and annotations
#' @export
#'
summarise_counts <- function(se) {
  
  se_dt <- assay(se) %>% as.data.table
  sample_ids <- colnames(se)
  
  se_dt[, mean_counts := rowMeans(.SD[, sample_ids, with = F])]
  se_dt[, sum_counts := rowSums(.SD[, sample_ids, with = F])]
  se_dt[, counted := .SD[, mean_counts > 0]]
    
  cbind(mcols(se) %>% as.data.table, se_dt)
}

#'
#' for checking the coverage
#' @param exons_dt data.table containing the counts and summary values as returned by summarise_counts()
#' @export
#'
coverage_ratio <- function(exons_dt) {
  total <- table(exons_dt$counted)
  not_counted <- total[which(row.names(total) == F)]
  counted <- total[which(row.names(total) == T)]
  
  as.numeric(counted/(counted + not_counted))
}
