% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{essential_lfc_qc}
\alias{essential_lfc_qc}
\title{Computes essential gene recovery AUC.}
\usage{
essential_lfc_qc(
  df,
  screens,
  negative_controls = NULL,
  append_to_negatives = FALSE
)
}
\arguments{
\item{df}{LFC dataframe.}

\item{screens}{List of screens generated with \code{add_screens}.}

\item{negative_controls}{List of negative control genes to append to default list of
non-essential genes (default NULL).}

\item{append_to_negatives}{Whether to append the negative controls specified in
negative_controls to the list of gold-standard human non-essential genes if TRUE or to
replace the list of non-essential genes if FALSE (default FALSE).}
}
\value{
Returns a dataframe with three columns for replicate name, essential AUC relative 
  to all other genes, and essential AUC relative to a specified set of non-essentials.
}
\description{
Computes area under the curve for ROC curves that measure how well each technical replicate
recovers signal for essential-targeting guides. Only computes AUC for guides that target 
essential genes twice, guides that target two different essential genes, or guides that 
target an essential gene and an intergenic region.
}
