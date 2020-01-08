#' Flatten data.tables with list entries
#' 
#' @param dt any object that can be coerced to data.table
#' @param id_columns all columns that do not contain any lists as entries
#' @param list_columns all columns that have list entries
#' @param debug print unlisting statement to be evaluated for each list_column
#' @importFrom rlang sym syms expr
#' @importFrom purrr reduce
#' @importFrom dplyr left_join
#' @import data.table
#' @export
#' 
unnest_dt <- function(dt, id_columns, list_columns, debug=FALSE) {
  
  dt <- as.data.table(dt)
  
  # unnest function for single dt
  unnest_dt_single <- function(tbl, col) {
    
    tbl <- as.data.table(tbl)
    col <- sym(col)
    clnms <- syms(setdiff(colnames(tbl), as.character(col)))
    #tbl <- as.data.table(tbl)
    
    cmd <- expr(tbl[, as.character(unlist(!!col)), by = list(!!!clnms)])
    if (isTRUE(debug)) print(cmd)
    tbl <- eval(cmd)
    colnames(tbl) <- c(as.character(clnms), as.character(col))
    
    tbl
  }
  
  # flatten all list columns independently
  flattened <- lapply(list_columns, function(col) {
    dt_sub <- dt[, .SD, .SDcols=c(id_columns,col)]
    unnest_dt_single(dt_sub, col)
  })
  
  # merge all flattened data.tables
  merged <- data.table(reduce(flattened, left_join, by = id_columns))
  merged
}

#' Summarise gene type
#' 
#' Add a column named general_gene_type to the data.table containing generalised gene type names 
#' @param dt data.table with at least a column containing specific gene types
#' @param column specifies the column of the specific gene type
#' @param mapping mapping of specific to general gene type as data.table or list of vectors
#' if kept NULL, then a basic mapping of pc (protein_coding) vs. npc is created
#' @param general_column name of general gene type column to add
#' @export
#'
generalise_gene_type <- function(dt, mapping = NULL, column = "gene_type", protein_coding = "protein_coding", general_column = "general_gene_type") {
  
  if (is.null(mapping)) {
    specific <- unique(dt[[column]])
    mapping <- data.table(specific, ifelse(specific == protein_coding, "coding", "non-coding"))
  } else if (class(mapping) == "list") {
    mapping <- rbindlist(mapply(data.table, mapping, names(mapping), SIMPLIFY = F))
  } else if (!class(mapping) %in% c("data.table")) {
    stop("invalid mapping table of type ", class(mapping))
  }
  setnames(mapping, c(column, general_column))
  
  merge(dt, mapping, by = column)
  
}

#' Plot multiple results with a loop
#' 
#' @param ods Outrider dataset
#' @param geneName name of gene to be aggregated over
#' @param feature_dt dt containing feature information, must be subsettable by gene_name
#' @param plot a plot function with params (ods, geneID)
#' @param export
#' 
multiplotLoop <- function(ods, ods_gene, geneName, feature_dt = NULL, plot = plotQQ, ...){
  
  if (is.null(feature_dt))
    feature_dt <- featureInfo(ods, columns = c("gene_name"))
  features <- feature_dt[gene_name == geneName, feature_id]
  show(feature_dt[feature_id %in% features])
  
  par(mfrow = c(2,2))
  for (id in features) {
    plot(ods, as.character(id), main = feature_dt[feature_id == id, type], ...)
  }
  
}

#' Make a clean copy of the original ODS
#' 
#' Removes old assays and keeps only original counts
#' Removes any meta data
#' 
#' @param ods OUTRIDER Dataset or RangedSummarizedExperiment
#' @param keep.meta inf TRUE keep mcols, otherwise remove (default FALSE)
#' @return copy of ods with 2 assays "counts" and "trueCounts"
#' @export
#' 
copyODS <- function(ods, keep.meta = T) {
  
  require(SummarizedExperiment)
  copy <- ods
  counts <- assays(copy)$counts
  
  assays(copy) <- list(counts = counts)
  assays(copy)$trueCounts <- counts
  
  if (!keep.meta) mcols(copy) <- NULL
  
  copy
}

#' converts GRanges object to bed file
#' 
#' @param gr GRanges object
#' @param filename for output bed file
#' @param invert strand if TRUE
#' @export
#'
gr_to_bed <- function(gr, filename, invertStrand = F) {
  
  if(invertStrand)
    strand(gr) <- ifelse(strand(gr) == '+', '-', '+')
  
  mcols(gr)$type <- paste(gr$feature, gr$location, sep = "_")
  types <- sort(unique(gr$type), decreasing = T)
  color_map <- data.table(type = types, 
                          col = apply(col2rgb(rainbow(length(types))),
                                      2, paste, collapse = ",")
  )
  setkey(color_map, type)
  
  dt <- data.table(
    chrom = as.character(seqnames(gr)), 
    chromStart = start(gr), 
    chromEnd = end(gr), 
    # name = paste(gr$transcript_id, gr$type),
    name = gr$transcript_id,
    score = gr$feature_id, 
    strand = as.character(strand(gr))
    ,
    thickStart = start(gr),
    thickEnd = start(gr),
    itemRgb = color_map[gr$type, col]
  )
  
  setorderv(dt, c("chrom", "chromStart"), order = c(1,1))
  write.table(dt, file=filename, quote=F, sep="\t", row.names=F, col.names=F)
  dt
}

## small helpers

#'
#' @return logical vector of all positions that have maximum value
#' @export
#' 
which_max <- function(x) which(x == max(x))
