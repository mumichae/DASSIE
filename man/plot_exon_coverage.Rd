% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_exon_coverage}
\alias{plot_exon_coverage}
\title{gene_type must be included in counts data table}
\usage{
plot_exon_coverage(
  counts,
  select = NULL,
  stratify = F,
  title = "Exon coverage",
  filename = NA,
  ...
)
}
\arguments{
\item{counts}{must be a data.table with genetic features in rows and samples in columns
must also contain at least gene type}

\item{select}{gene types}

\item{stratify}{for stratifying according to gene type}

\item{title}{displayed in plot}

\item{filename}{if not specified, the plot will not be saved}
}
\description{
gene_type must be included in counts data table
}
