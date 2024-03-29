% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{get_all_results}
\alias{get_all_results}
\title{Get all OUTRIDER results}
\usage{
get_all_results(
  ods,
  assays = c("pValue", "zScore"),
  meta = c("feature_id"),
  id.as.int = F
)
}
\arguments{
\item{ods}{OUTRIDER data set}

\item{assays}{to subset which values are to be extracted. Vector can either be empty or contain only names of existing assays in the ods (names(assays(ods)))}

\item{meta}{metadata columns to include}

\item{id.as.int}{wether gene name should be coerced to integer (useful for feature_ids)}
}
\description{
Efficiently retrieve only certain assays for all features
}
