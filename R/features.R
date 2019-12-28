#' Retrieve info as data.table
#' 
#' Get data.table of metadata from RangedSummarizedExperiment or OutriderDataSet object
#' 
#' @param gr GRanges object or any object containing it (i.e. RangedSummarizedExperiment, OutriderDataSet), must support mcols function
#' @param columns a character vector for selecting only a part of the mcols. If NULL then all the columns will be included
#' @param summary include summary per feature over all samples, if TRUE
#' @param rownumbers if TRUE, the row numbers will be included as columns "row"
#' @export
#'
featureInfo <- function(gr, columns = NULL, summary = F, rownumbers = F) {
  
  # select columns
  if (is.null(columns))
    feature_dt <- as.data.table(mcols(gr))
  else {
    columns <- unique(c(columns, c("feature_id", "feature", "location")))
    feature_dt <- as.data.table(mcols(gr)[, columns])
  }
  
  # add type for easy subsetting
  feature_dt[, type := paste(feature, location, sep = "_")]
  feature_dt[, type := factor(type, unique(type), order = T)]
  
  # add summary info
  if (summary) {
    feature_dt[, meanCounts := rowMeans(assay(gr))]
  }
  
  # add row numbers
  if (rownumbers)
    feature_dt[, row := 1:.N]
  
  feature_dt
}

#' maps the exons to the transcripts
#' 
#'@param txdb annotation database
#'@return
#'
#'@importFrom dplyr %>%
#'@import GenomicFeatures
#'@export
#'
exon_to_tx <- function(txdb) {
  
  exons_tx <- exonsBy(txdb, by = "tx") %>% unlist
  tx_mapping <- mcols(transcripts(txdb)) %>% as.data.table
  
  tx_ids <- as.integer(names(exons_tx))
  tx_names <- tx_mapping[tx_ids, tx_name]
  
  exons_tx_dt <- exons_tx %>% as.data.table
  exons_tx_dt[, transcript_idx := tx_ids]
  exons_tx_dt[, transcript_name := tx_names]
  
  exons_tx_dt
}

#' Extract regions around TSS and PAS and returns a GRanges object
#' 
#' Convenience function for TSS and PAS regions of same length
#' 
#' @param exons data.table containing at least the exon coordinates (for GRanges), "transcript_id" and "exon_id" for each exon. Must be strand-specific.
#' @param window number of bp that the regions should have
#' @export
#' 
windowRegions <- function(exons, window = 500) {
  
  tss <- tssRegions(exons, window)
  pas <- pasRegions(exons, window)
  
  feature_dt <- rbind(tss, pas)
  # feature_dt <- merge(exons[, -c("start", "end")], feature_dt, by = c("transcript_id", "exon_id")) # add annotation info
  feature_dt[, feature_id := 1:.N]
  
  regions <- makeGRangesFromDataFrame(feature_dt, keep.extra.columns = T)
  names(regions) <- regions$feature_id
  
  regions
}

#' Transcription start site regions
#'
#' Get regions of specified length around the TSS
#' @param exons data.table containing at least the exon coordinates (for GRanges), "transcript_id" and "exon_id" for each exon. Must be strand-specific.
#' @param window number of bp for regions
#' @return data.table with coordinates and some metadata on the new regions
#' exon_id and transcript_id are necessary for joining the features with the rest of the exon annotation.
#' @export
#'
tssRegions <- function(exons, window) {
  
  ## strand: +
  setorder(exons, start)
  plus <- exons[strand == "+", .SD[1], by = "transcript_id"] # first exons
  plus[, site := .SD[, start]] # TSS based on first exons
  
  US <- plus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                 start = .SD$site - window, end = .SD$site - 1,
                 feature = "TSS", location = "upstream", site)]
  DS <- plus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                 start = .SD$site, end = .SD$site + window - 1,
                 feature = "TSS", location = "downstream", site)]
  plus <- rbind(US, DS)
  
  ## strand: -
  setorder(exons, -end)
  minus <- exons[strand == "-", .SD[1], by = "transcript_id"]
  minus[, site := .SD[, end]]
  
  US <- minus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                  start = .SD$site + 1, end = .SD$site + window,
                  feature = "TSS", location = "upstream", site)]
  DS <- minus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                  start = .SD$site - window + 1, end = .SD$site,
                  feature = "TSS", location = "downstream", site)]
  minus <- rbind(US, DS)
  
  rbind(plus, minus)
}


#' Polyadenylation site regions
#' 
#' Get regions of specified length around the PAS
#' @param exons data.table containing at least the exon coordinates (for GRanges), "transcript_id" and "exon_id" for each exon. Must be strand-specific.
#' @param window number of bp for regions
#' @return data.table with coordinates and some metadata on the new regions
#' exon_id and transcript_id are necessary for joining the features with the rest of the exon annotation.
#' @export
#'
pasRegions <- function(exons, window) {
  
  ## strand: +
  setorder(exons, -start)
  plus <- exons[strand == "+", .SD[1], by = "transcript_id"] # last exons
  plus[, site := .SD[, end]] # PAS
  
  US <- plus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                 start = .SD$site - window + 1, end = .SD$site,
                 feature = "PAS", location = "upstream", site)]
  DS <- plus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                 start = .SD$site + 1, end = .SD$site + window,
                 feature = "PAS", location = "downstream", site)]
  plus <- rbind(US, DS)
  
  ## strand: -
  setorder(exons, end)
  minus <- exons[strand == "-", .SD[1], by = "transcript_id"]
  minus[, site := .SD[, start]]
  
  US <- minus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                  start = .SD$site, end = .SD$site + window - 1,
                  feature = "PAS", location = "upstream", site)]
  DS <- minus[, .(exon_id, transcript_id, transcriptID, gene_name, seqnames, strand,
                  start = .SD$site - window, end = .SD$site - 1,
                  feature = "PAS", location = "downstream", site)]
  minus <- rbind(US, DS)
  
  rbind(plus, minus)
}


#'
#' simple function that takes the first hit after sorting
#' @param dt containing exon info
#' @param by sorting criterium (transcript or gene)
#' @export
#'
get_first_exons <- function(dt, by) {
  setorder(dt, start)
  plus <- dt[strand == "+", .SD[1], by = by]
  setorder(dt, -end)
  minus <- dt[strand == "-", .SD[1], by = by]
  rbind(plus, minus)
}

#'
#' simple function that takes the first hit after sorting
#' @param dt containing exon info
#' @param by sorting criterium (transcript or gene)
#' @export
#'
get_last_exons <- function(dt, by) {
  # dt <- dt[seqnames != "MT" & seqnames != "chrM"] # remove mitochondrial DNA
  setorder(dt, -start)
  plus <- dt[strand == "+", .SD[1], by = by]
  setorder(dt, end)
  minus <- dt[strand == "-", .SD[1], by = by]
  rbind(plus, minus)
}