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

#' Get exons from txdb
#'
#' @param txdb
#' @param filter_length filter transcripts
#'               'tx_length': by length of transcript
#'               'gene_length': by length of gene 
#' @param filter_exon keep transcript with longest exons
#' @return exons in a data.table convertible to GRanges object, contains
#'         transcript and gene mapping
#' @import data.table
#' @importFrom GenomicFeatures exons
#' @export
#'
getExons <- function(txdb, filter_length = "tx_length", filter_exon = TRUE) {
    
    longest_tx <- getLongestTranscript(txdb, filter_length, filter_exon)
    
    # exon mapping
    message("extract exons into data.table")
    range_columns <- c("seqnames", "start", "end", "strand")
    id_cols       <- c("exon_id", "exon_name")
    list_cols     <- c("tx_id", "tx_name", "gene_id")
    
    exons <- exons(txdb, columns = c(id_cols, list_cols),
                   filter=list(tx_id=longest_tx))
    
    exon_dt <- unnest_dt(as.data.table(exons), 
                         id_columns = c(range_columns, id_cols),
                         list_columns = list_cols)
    exon_dt
}

#' 
#' @import data.table
#' @importFrom GenomicFeatures transcriptsBy transcriptLengths
#' 
getLongestTranscript <- function(txdb, filter_length = "tx_length", filter_exon = TRUE) {
    
    message("get transcript by gene")
    tx_by_gene <- as.data.table(transcriptsBy(txdb, by = "gene"))
    
    if (filter_length == "tx_length") {
        # longest by tx length
        message("select longest transcripts by ", filter_length)
        tx_lengths <- as.data.table(transcriptLengths(txdb))
        longest_tx <- tx_lengths[, .SD[which_max(tx_len)], by = gene_id]
      
    } else if (filter_length == "gene_length") {
        # longest by width = longest by position in genome
        message("select longest transcripts by ", filter_length)
        max_width <- lapply(width(tx_by_gene), which_max)
        longest_gene <- as.data.table(tx_by_gene[max_width])
        
        tx_lengths <- as.data.table(transcriptLengths(txdb))
        longest_tx <- tx_lengths[tx_id %in% longest_gene$tx_id]
    } else {
        longest_tx <- tx_by_gene
    }
    
    if (isTRUE(filter_exon)) {
        # longest by number of exons
        message("select longest transcripts by number of exons")
        longest_tx <- longest_tx[, .SD[which_max(nexon)], by = gene_id]
    }
    
    longest_tx$tx_id
}

#' Extract regions around TSS and PAS and returns a GRanges object
#' 
#' Convenience function for TSS and PAS regions of same length
#' 
#' @param exon_dt data.table convertible into GRanges, and contain a transcript
#'                column that is specified in tx_id. Data must be strand-specific.
#' @param window number of bp that the regions should have
#' @export
#' 
getDassieRegions <- function(exon_dt, tssWindow = 500, pasWindow = 1000, 
                             tx_id="tx_id", remove_redundant = TRUE) {
    
    tss <- tssRegions(exon_dt, window=tssWindow, tx_id=tx_id)
    pas <- pasRegions(exon_dt, window=pasWindow, tx_id=tx_id)
    
    feature_dt <- rbind(tss, pas)
    feature_dt[, feature_id := 1:.N]
    
    # Remove redundant regions
    if (isTRUE(remove_redundant)) {
      message("removing duplicated features")
      feature_dt <- feature_dt[, .SD[1], 
                               by = c("seqnames", "start", "end", 
                                      "strand", "feature", "location")]
    }
    
    # Convert to GenomicRanges and save
    regions <- makeGRangesFromDataFrame(feature_dt, keep.extra.columns = T)
    names(regions) <- regions$feature_id
    
    regions
}

#' Transcription start site regions
#'
#' Get regions of specified length around the TSS
#' @param exons_dt data.table of GRanges object with columns "tx_id" for each
#'                 exon for merging by transcripts. Ranges must be strand-specific.
#' @param window number of bp for regions
#' @return GRanges object with upstream and downstream TSS regions
#' @import data.table
#' @export
#'
tssRegions <- function(exon_dt, window, tx_id="tx_id") {
  
  #TODO: deal with non-stranded data
  
  tssPlus <- function() {
    # first exons
    setorder(exon_dt, start)
    plus <- exon_dt[strand == "+", .SD[1], by = tx_id]
    # TSS
    plus[, site := start]
    
    # ranges upstream of site
    US <- copy(plus)
    US[, start := site - window]
    US[, end := site - 1]
    US[, location := "upstream"]
    
    # ranges downstream of site
    DS <- copy(plus)
    DS[, start := site]
    DS[, end := site + window - 1]
    DS[, location := "downstream"]
    
    rbind(US, DS)
    
  }
  
  tssMinus <- function() {
      # first exons
      setorder(exon_dt, -end)
      minus <- exon_dt[strand == "-", .SD[1], by = tx_id]
      
      # TSS
      minus[, site := end]
      
      # ranges upstream of site
      US <- copy(minus)
      US[, start := site + 1]
      US[, end := site + window]
      US[, location := "upstream"]
      
      # ranges downstream of site
      DS <- copy(minus)
      DS[, start := site - window + 1]
      DS[, end := site]
      DS[, location := "downstream"]
      
      rbind(US, DS)
  }
  
  tssFeatures <- rbind(tssPlus(), tssMinus())
  tssFeatures[, feature := "TSS"]
  
  tssFeatures
  
}


#' Polyadenylation site regions
#' 
#' Get regions of specified length around the PAS
#' @param exons_dt data.table of GRanges object with columns "tx_id" for each
#'                 exon for merging by transcripts. Ranges must be strand-specific.
#' @param window number of bp for regions
#' @return data.table with coordinates and some metadata on the new regions
#' exon_id and transcript_id are necessary for joining the features with the rest of the exon annotation.
#' @import data.table
#' @export
#'
pasRegions <- function(exon_dt, window, tx_id="tx_id") {
  
  pasPlus <- function() {
      # last exons
      setorder(exon_dt, -start)
      plus <- exon_dt[strand == "+", .SD[1], by = tx_id]
      
      # PAS
      plus[, site := end]
      
      # ranges upstream of site
      US <- copy(plus)
      US[, start := site - window + 1]
      US[, end := site]
      US[, location := "upstream"]
      
      # ranges downstream of site
      DS <- copy(plus)
      DS[, start := site + 1]
      DS[, end := site + window]
      DS[, location := "downstream"]
      
      rbind(US, DS)
  }
  
  pasMinus <- function() {
      # last exons
      setorder(exon_dt, end)
      minus <- exon_dt[strand == "-", .SD[1], by = tx_id]
      
      # PAS
      minus[, site := start]
      
      # ranges upstream of site
      US <- copy(minus)
      US[, start := site]
      US[, end := site + window - 1]
      US[, location := "upstream"]
      
      # ranges downstream of site
      DS <- copy(minus)
      DS[, start := site - window]
      DS[, end := site - 1]
      DS[, location := "downstream"]
      
      rbind(US, DS)
  }
  
  pasFeatures <- rbind(pasPlus(), pasMinus())
  pasFeatures[, feature := "PAS"]
  
  pasFeatures
  
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