#' Visualise coverage for multiple genes using average counts
#' @param cvg_list a named list of coverage Rle's
#' @param regions GRanges object containing the feature regions, must contain a column "site"
#' @param strand either "+" or "-"
#' @param feature type of feature: "TSS" or "PAS"
#' @param window number of positions upstream and downstream of feature site
#' @param log log representation of counts
#' @export
#' @deprecated
#'
plotCoverageHeatmap <- function(cvg_list, regions, strand = "+", feature = "TSS", window = 500, log = T, parallel = F) {
  
  cvg <- coverageMatrix(cvg_list, regions, strand, feature, window, parallel)
  drawCoverageHeatmap(cvg, strand = strand, feature = feature, window = window, log = log)
  
}

#' Coverage Matrix
#' 
#' Generates a matrix with average per-base counts around TSS or PAS of given transcripts
#' @param cvg_list list of coverage Rles
#' @param regions data.table containing the feature regions, must contain a column "site"
#' @param strand either "+" or "-"
#' @param feature type of feature: "TSS" or "PAS"
#' @param window number of positions upstream and downstream of feature site
#' @export
#'
coverageMatrix <- function(cvg_list, regions, strand, feature, window, parallel = F) {
  
  sites <- regions[feature == feature & strand == strand, .(seqnames, feature, gene_name)]
  
  if (parallel) {
    require(BiocParallel)
    total_cvg <- bplapply(1:nrow(sites), function(i)
      colMeans2(coverageMatrixPerSite(cvg_list, regions, gene = sites$gene_name[i], feature = sites$feature[i], window)))
  } else {
    total_cvg <- lapply(1:nrow(sites), function(i)
      colMeans2(coverageMatrixPerSite(cvg_list, regions, gene = sites$gene_name[i], feature = sites$feature[i], window)))
  }
  do.call(rbind, total_cvg) # as matrix
}

#' Coverage Heatmap per site
#' 
#' Visualise coverage across all samples for a given gene name
#' @param cvg_matrix matrix to be plotted (created by coverMatrixPerSite)
#' @param regions data.table containing the feature regions, must contain a column "site"
#' @param gene gene name
#' @param strand either "+" or "-"
#' @param feature type of feature: "TSS" or "PAS"
#' @param window number of positions upstream and downstream of feature site
#' @param log log representation of counts
#' @export
#'
plotPerSiteCoverageHeatmap <- function(cvg_matrix, regions, gene, feature = "TSS", window = 500, log = T) {
  
  ft <- feature
  site <- unique(regions[gene_name == gene & feature == ft, .(seqnames, strand, site)])
  if(nrow(site) != 1) 
    stop(paste("slice should only have 1 row but has", nrow(site)))
  
  drawCoverageHeatmap(cvg_matrix, site$strand, feature, window, gene, log)
  
}

#' Coverage Matrix per Site
#' 
#' Generates a matrix with per-base counts around a specific TSS or PAS for all samples
#' @param cvg_list a named list of coverage Rle's
#' @param regions data.table containing the feature regions, must contain a column "site"
#' @param gene name of gene
#' @param strand either "+" or "-"
#' @param feature type of feature: "TSS" or "PAS"
#' @param window number of positions upstream and downstream of feature site
#' @export
#'
coverageMatrixPerSite <- function(cvg_list, regions, gene, feature = "TSS", window = 500) {
  ft <- feature
  site <- unique(regions[gene_name == gene & feature == ft, .(seqnames, strand, site)])
  if(nrow(site) != 1) 
    stop(paste("slice should only have 1 row but has", nrow(site), "for gene", gene))
  
  start <- site$site - window
  end <- site$site + window
  t(sapply(cvg_list, function(cvg) as.vector(cvg[[site$seqnames]][start:end])))
}

#' Draw a heatmap
#' 
#' Wraps the actual plotting, only dependent on matrix and meta info
#' @param cvg coverage matrix
#' @param strand either "+" or "-"
#' @param feature either "TSS" or "PAS"
#' @param window positions upstream/downstream
#' @param transcript ID of the transcript to be displayed in the plot for single-Transcript plots
#' @param log TRUE or FALSE for log10-transformation of counts + 1
#' @return a heatmap depending on the input (per site if transcript is specified or for multiple transcripts).
#' @export
#' 
drawCoverageHeatmap <- function(cvg_matrix, strand, feature, window, gene = NULL, log) {
  
  require(ComplexHeatmap)
  if (sum(cvg_matrix) == 0) stop("no counts in matrix, so plot not drawn")
  
  n <- max(cvg_matrix)
  n_mid <- n/2
  # use log-scale
  if (log) {
    cvg_matrix <- log10(cvg_matrix+1)
    n <- log10(n)
    n_mid <- log10(n_mid)
  }
  count_label <- ifelse(log, "log10(counts + 1)", "read counts")
  heat_colors <- circlize::colorRamp2(c(0, n/2, n), c("red", "yellow", "white"))
  
  # annotation
  if (strand == "+") {
    loc_df <- data.frame(location = rep(c("upstream", "downstream"), c(window, window + 1)))
  } else if (strand == "-") {
    loc_df <- data.frame(location = rep(c("downstream", "upstream"), c(window, window + 1)))
  } else {
    stop(paste0("invalid strand: ", strand))
  }
  loc_colors <- list(location = c("upstream" = "chartreuse3", "downstream" = "darkorchid"))
  ha <- HeatmapAnnotation(df = loc_df, col = loc_colors)
  
  # adjust header according to viz type
  title <- ifelse(is.null(gene), yes = paste(feature, "on", strand, "strand."),
                  no = paste(gene, feature, "on", strand, "strand."))
  title <- paste(title, "window:", window)
  
  ht <- Heatmap(cvg_matrix, column_title = title,
                name = count_label, 
                col = heat.colors(20),
                cluster_rows = T, cluster_columns = F, 
                top_annotation = ha,
                heatmap_legend_param = list(legend_direction = "horizontal",
                                            title_position = "lefttop"),
                row_names_gp = gpar(fontsize = 8),
                show_row_dend = F
                )
                          
  draw(ht, heatmap_legend_side = "bottom")
  
}
