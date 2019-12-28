#' Expression Ranks with marked sample
#' 
#' @param ods OUTRIDER data set
#' @param gene_id ID of ods features
#' @param padj cutoff for marking most significant
#' @param normalized if T then normalized counts will be used
#' @param sampleIDs sample IDs to be additionally marked, if not already marked signif
#' @export
#' 
plotExpressionRank2 <- function(ods, id, sampleIDs = c(), padj = 0.05, normalized = T, main = NULL, log = T) {
  
  require(ggplot2)
  require(cowplot)
  dt <- get_all_results(ods, assays = c("pValue", "padjust", "counts", "normalizationFactors"), id.as.int = F)
  dt <- dt[feature_id == id]
  dt[, marked := "gray40"]
  dt[, marked := ifelse(sampleID %in% sampleIDs, "chartreuse", marked)]
  dt[, marked := ifelse(padjust <= padj, "firebrick", marked)]
  
  if(isTRUE(normalized)) 
    dt[, counts := counts / normalizationFactors * mean(normalizationFactors)]
  dt[order(counts), rank := 1:.N]
  
  if (is.null(main))
    main <- paste("Expression ranks for", id)
  main <- paste(main, ifelse(normalized, "(normalized)", "(raw counts)"))
  
  p <- ggplot(dt, aes(x = rank, y = counts, col = marked)) + 
    geom_point(size = 2) + 
    ggtitle(main) +
    scale_color_identity() +
    background_grid(major = "xy")
  
  if (log) {
    p <- p + scale_y_log10()
  }
  
  p
}

#' QQ plots with marked samples
#' 
#' @param ods OUTRIDER data set
#' @param gene_id ID of ods features
#' @param padj cutoff for marking most significant
#' @param sampleIDs sample IDs to be additionally marked, if not already marked signif
#' @export
#' 
plotQQ2 <- function(ods, id, sampleIDs = c(), padj = 0.05, main = NULL, ci = 0.95) {
  
  require(ggplot2)
  dt <- get_all_results(ods[id], assays = c("pValue", "padjust"), id.as.int = F)
  aberrantEvent <- aberrant(ods[id], padjCutoff = padj, by = "sample")
  dt[, col := "gray40"]
  dt[, col := ifelse(sampleID %in% sampleIDs, "chartreuse", col)]
  dt[, col := ifelse(aberrantEvent, "firebrick", col)]

  # observed
  dt[, obs := -log10(pValue)]
  dt[is.na(obs) | is.infinite(obs), obs:=1]
  dt[obs == 0, obs := min(dt[obs != 0, obs]) * 1e-2]
  setorder(dt, -obs)
  
  # expected
  dt[, exp := -log10(ppoints(.N))]
  
  # confidence band https://gist.github.com/slowkow/9041570
  dt[, clower := -log10(qbeta((1 - ci)/2, 1:.N, .N - 1:.N + 1))]
  dt[, cupper := -log10(qbeta((1 + ci)/2, 1:.N, .N - 1:.N + 1))]
  
  if (is.null(main))
    main <- paste("QQ Plot for", id)
  
  ggplot(dt, aes(x = exp, y = obs, col = col)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.8, col = "firebrick") +
    geom_line(aes(exp, cupper), linetype = 2, col = "gray") +
    geom_line(aes(exp, clower), linetype = 2, col = "gray") +
    geom_point(size = 2) +
    scale_color_identity() +
    labs(
      title = main,
      x = "-log10(expected p-value)",
      y = "-log10(observed p-value)"
    )
}

#' gene_type must be included in counts data table
#' 
#' @param counts must be a data.table with genetic features in rows and samples in columns
#'     must also contain at least gene type
#' @param select gene types
#' @param stratify for stratifying according to gene type
#' @param title displayed in plot
#' @param filename if not specified, the plot will not be saved
#'  
plot_exon_coverage <- function(counts, select = NULL, stratify = F, title = "Exon coverage", filename = NA, ...) {
  
  if (is.null(select)) select <- unique(counts$gene_type) # all gene types
  
  if (stratify) { 
    by_type <- as.data.table(table(counts[gene_type %in% select, c("gene_type", "counted"), with = F]))
    ggplot(by_type, aes(gene_type, N, fill = counted)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(
        aes(label = N, y = N + 100),
        position = position_dodge(1),
        vjust = 0
      ) +
      labs(x = "Gene type", y = "Number of exons", title = paste0(title, " by gene type")) + 
      theme_bw() + scale_fill_ptol(name = "", labels = c("0 reads", "> 0 reads"))
  }
  else {
    total_counted <- as.data.table(table(counts[gene_type %in% select, counted]))
    
    ggplot(total_counted[c(2,1),], aes(V1, N, fill = V1)) +
      geom_bar(stat = "identity") +
      geom_text(
        aes(label = N, y = N + 100),
        position = position_dodge(1),
        vjust = 0
      ) +
      labs(x = "", y = "Number of exons", title = title) + 
      theme_bw() +
      scale_fill_ptol(name = "counted")
  }
  
  if(!is.na(filename)) 
    ggsave(filename, ...)
  else 
    show(last_plot())
  
}

#'
#'
#' plots number of exons according to gene type
#' @param exons must be a data.table with genetic features in rows and samples in columns
#'     must also contain at least gene type
#' @param select gene types
#' @param title displayed in plot
#' @param filename if not specified, the plot will not be saved
#'
plot_exon_distribution <- function(exons, select = NULL, title = "Number of exons by gene type", filename = NA, ...) {
  
  if (is.null(select)) 
    select <- unique(exons[, gene_type]) # all gene types
  
  exon_count_by_type <- exons[gene_type %in% select, .N, by = c("gene_type", "gene_name")]
  
  ggplot(exon_count_by_type, aes(gene_type, N)) + 
    geom_boxplot() + 
    scale_y_log10(breaks = c(0,1,2, 5,10, 20, 50, 100, max(exon_count_by_type[,N])) ) +
    labs(x = "Gene type", y = "Number of exons", title = title) + 
    theme_bw()
  
  if(!is.na(filename)) 
    ggsave(filename, ...)
  else 
    show(last_plot())
}

#'
#' Boxplots the distribution of mean counts per gene compared between gene types
#' @param exons must be a data.table with genetic features in rows and samples in columns
#'     must also contain at least gene type
#' @param select gene types
#' @param title displayed in plot
#' @param filename if not specified, the plot will not be saved
#' 
plot_mean_counts <- function(exons, select = NULL, title = "Mean exon counts", filename = NA, ...) {
  
  if (is.null(select)) 
    select <- unique(exons[, gene_type]) # all gene types
  
  counts_by_type <- exons[gene_type %in% select, c("gene_type", "mean_counts")]
  
  ggplot(counts_by_type, aes(gene_type, mean_counts)) +
    geom_boxplot() + 
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, max(counts_by_type[, mean_counts]))) +
    labs(x = "Gene type", y = "Mean read count", title = title) + 
    theme_bw() + 
    scale_fill_ptol() + 
    stat_summary(fun.data = function(x) c(y = median(x), label = length(x)), geom = "text")
    
    if(!is.na(filename)) 
      ggsave(filename, ...)
  else 
    show(last_plot())
  
}